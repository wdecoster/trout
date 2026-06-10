use std::collections::HashMap;
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::sync::Arc;

pub struct Allele {
    /// Sample name. Every allele read from one VCF shares the same name, so we store it as a
    /// reference-counted `Arc<str>` and clone the pointer (not the bytes) per allele. With ~49M
    /// alleles across a cohort this turns a per-allele heap string into one allocation per file.
    pub sample: Arc<str>,
    /// Allele sequence as `Box<str>` rather than `String`: we never grow it after parsing, so the
    /// extra capacity word `String` carries is pure overhead at this multiplicity.
    pub seq: Box<str>,
}

/// All alleles seen at a locus plus the REF sequence from the VCF.
/// REF is the same across all input VCFs at the same coordinates (genome reference),
/// so the first writer wins on merge.
pub struct LocusData {
    pub ref_seq: Box<str>,
    pub alleles: Vec<Allele>,
}

/// One finalised locus: its `(chrom, pos, end)` key paired with the gathered alleles.
pub type LocusEntry = ((Box<str>, u32, u32), LocusData);

/// A single parsed VCF data line: the alleles one sample contributes at one locus. Held as the
/// "head" of a stream during the k-way merge.
struct Record {
    /// Rank of `chrom` in genome/contig order (from the `##contig` header lines, or assigned on
    /// first sight). Drives the merge ordering so output comes out in reference order.
    rank: u32,
    chrom: Box<str>,
    pos: u32,
    end: u32,
    ref_seq: Box<str>,
    /// The chosen allele sequences for this sample at this locus, after `--min-support` filtering.
    /// Never empty: lines whose alleles are all filtered out are skipped during parsing.
    seqs: Vec<Box<str>>,
}

/// One open VCF, decompressed on the fly, with its next unconsumed record peeked as `head`.
struct Stream {
    reader: Box<dyn BufRead>,
    sample: Arc<str>,
    head: Option<Record>,
}

/// Streaming k-way merge over a cohort of coordinate-sorted VCFs.
///
/// All files are opened at once and advanced in lockstep by genomic coordinate. `next_batch`
/// gathers the alleles for the next run of loci across every file and returns them in reference
/// order, so peak memory is bounded by one batch of loci rather than the whole cohort. This is
/// what lets the tool scale to large repeat catalogs: cost is O(batch × samples), not
/// O(loci × samples).
///
/// File-handle note: every input VCF is held open for the duration. With a few hundred samples
/// this is fine under the default `ulimit -n` (1024); for larger cohorts raise it (`ulimit -n`)
/// — `open` surfaces a clear error if the limit is hit.
pub struct VcfMerger {
    streams: Vec<Stream>,
    /// chrom -> rank in genome order. Seeded from `##contig` headers; grows if a record names a
    /// contig that no header declared (assigned the next rank on first sight).
    contig_rank: HashMap<Box<str>, u32>,
    next_rank: u32,
    min_support: Option<u32>,
    /// Total alleles dropped by `--min-support` so far. Read after the merge is drained.
    pub n_dropped: usize,
}

impl VcfMerger {
    /// Open every VCF, parse its header (contig order + sample name), and prime the first record
    /// of each. Files that fail to open or have no usable header are skipped with a warning.
    pub fn open(
        paths: &[PathBuf],
        min_support: Option<u32>,
    ) -> Result<VcfMerger, Box<dyn std::error::Error>> {
        // Genome/contig order, accumulated across files. The first file establishes the order;
        // contigs only later files declare are appended.
        let mut contig_rank: HashMap<Box<str>, u32> = HashMap::new();
        let mut next_rank: u32 = 0;

        let mut streams: Vec<Stream> = Vec::with_capacity(paths.len());
        for path in paths {
            match open_stream(path, &mut contig_rank, &mut next_rank) {
                Ok(stream) => streams.push(stream),
                Err(e) => {
                    // EMFILE/ENFILE means we ran out of file descriptors — almost certainly the
                    // ulimit rather than a bad file, so say so explicitly.
                    if e.raw_os_error() == Some(24) || e.raw_os_error() == Some(23) {
                        return Err(format!(
                            "ran out of open file handles after {opened} of {total} VCFs ({e}). \
                             trout holds every input open at once and already raised its soft limit \
                             to the hard cap at startup, so the hard limit itself is too low for \
                             {total} files. Ask your admin to raise the hard limit (e.g. \
                             /etc/security/limits.conf `nofile`), or run trout on fewer files at a \
                             time and combine the per-locus outputs.",
                            opened = streams.len(),
                            total = paths.len(),
                        )
                        .into());
                    }
                    eprintln!("Warning: failed to open {}: {}", path.display(), e);
                }
            }
        }

        let mut merger = VcfMerger {
            streams,
            contig_rank,
            next_rank,
            min_support,
            n_dropped: 0,
        };

        // Prime each stream's first record.
        let VcfMerger {
            streams,
            contig_rank,
            next_rank,
            min_support,
            n_dropped,
        } = &mut merger;
        for s in streams.iter_mut() {
            match next_record(&mut s.reader, *min_support, contig_rank, next_rank) {
                Ok((head, dropped)) => {
                    *n_dropped += dropped;
                    s.head = head;
                }
                Err(e) => {
                    eprintln!("Warning: error reading {}: {}", s.sample, e);
                    s.head = None;
                }
            }
        }

        Ok(merger)
    }

    /// Sample names of all successfully opened input VCFs (one per file). This is the set of
    /// identifiers `--samples` is matched against.
    pub fn sample_names(&self) -> Vec<Arc<str>> {
        self.streams.iter().map(|s| Arc::clone(&s.sample)).collect()
    }

    /// Gather up to `batch_size` loci in genome order, draining the merge as it goes. Returns an
    /// empty vec once every file is exhausted. Each entry is `((chrom, pos, end), LocusData)`.
    pub fn next_batch(&mut self, batch_size: usize) -> Vec<LocusEntry> {
        let VcfMerger {
            streams,
            contig_rank,
            next_rank,
            min_support,
            n_dropped,
        } = self;

        let mut out: Vec<LocusEntry> = Vec::new();

        while out.len() < batch_size {
            // Smallest head across all streams = the next locus to finalise.
            let mut min_key: Option<(u32, u32, u32, &str)> = None;
            for s in streams.iter() {
                if let Some(h) = &s.head {
                    let cand = (h.rank, h.pos, h.end, h.chrom.as_ref());
                    if min_key.is_none_or(|m| cmp_key(cand, m) == std::cmp::Ordering::Less) {
                        min_key = Some(cand);
                    }
                }
            }
            let Some((_, mpos, mend, mchrom)) = min_key else {
                break; // all streams drained
            };
            let mchrom: Box<str> = Box::from(mchrom);

            // Gather every stream sitting on this exact locus, advancing each as we consume it.
            let mut ref_seq: Option<Box<str>> = None;
            let mut alleles: Vec<Allele> = Vec::new();
            for s in streams.iter_mut() {
                while s
                    .head
                    .as_ref()
                    .is_some_and(|h| h.pos == mpos && h.end == mend && h.chrom == mchrom)
                {
                    let rec = s.head.take().unwrap();
                    if ref_seq.is_none() {
                        ref_seq = Some(rec.ref_seq);
                    }
                    for seq in rec.seqs {
                        alleles.push(Allele {
                            sample: Arc::clone(&s.sample),
                            seq,
                        });
                    }
                    match next_record(&mut s.reader, *min_support, contig_rank, next_rank) {
                        Ok((head, dropped)) => {
                            *n_dropped += dropped;
                            s.head = head;
                        }
                        Err(e) => {
                            eprintln!("Warning: error reading {}: {}", s.sample, e);
                            s.head = None;
                        }
                    }
                }
            }

            out.push((
                (mchrom, mpos, mend),
                LocusData {
                    ref_seq: ref_seq.unwrap_or_default(),
                    alleles,
                },
            ));
        }

        out
    }
}

/// Total order over `(rank, pos, end, chrom)`. The chrom string only breaks ties between contigs
/// that share a rank (i.e. both assigned the fallback rank), keeping the merge deterministic.
fn cmp_key(a: (u32, u32, u32, &str), b: (u32, u32, u32, &str)) -> std::cmp::Ordering {
    a.0.cmp(&b.0)
        .then_with(|| a.1.cmp(&b.1))
        .then_with(|| a.2.cmp(&b.2))
        .then_with(|| cmp_chrom(a.3, b.3))
}

/// Natural chromosome order: numeric contigs ascending, then non-numeric lexicographic. Only used
/// as a tie-break for contigs not declared in any `##contig` header.
fn cmp_chrom(a: &str, b: &str) -> std::cmp::Ordering {
    let a = a.strip_prefix("chr").unwrap_or(a);
    let b = b.strip_prefix("chr").unwrap_or(b);
    match (a.parse::<u64>(), b.parse::<u64>()) {
        (Ok(na), Ok(nb)) => na.cmp(&nb),
        (Ok(_), Err(_)) => std::cmp::Ordering::Less,
        (Err(_), Ok(_)) => std::cmp::Ordering::Greater,
        (Err(_), Err(_)) => a.cmp(b),
    }
}

fn reader(path: &Path) -> Result<Box<dyn BufRead>, std::io::Error> {
    let file = std::fs::File::open(path)?;
    if path.to_string_lossy().ends_with(".gz") {
        Ok(Box::new(std::io::BufReader::new(
            flate2::read::MultiGzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(std::io::BufReader::new(file)))
    }
}

/// Open one VCF and consume its header: record `##contig` order into `contig_rank` and pull the
/// sample name from the `#CHROM` line. The returned reader is positioned at the first data line.
fn open_stream(
    path: &Path,
    contig_rank: &mut HashMap<Box<str>, u32>,
    next_rank: &mut u32,
) -> Result<Stream, std::io::Error> {
    let mut reader = reader(path)?;
    let mut sample: Arc<str> = Arc::from("unknown");
    let mut line = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break; // empty / header-only file
        }
        if line.starts_with("##") {
            if let Some(id) = parse_contig_id(&line) {
                contig_rank.entry(Box::from(id)).or_insert_with(|| {
                    let r = *next_rank;
                    *next_rank += 1;
                    r
                });
            }
            continue;
        }
        if line.starts_with('#') {
            // #CHROM ... FORMAT <SAMPLE>: the sample name is the last column.
            let name = line.trim_end().rsplit('\t').next().unwrap_or("unknown");
            sample = Arc::from(name);
            break; // header done; reader now at the first data line
        }
        // A data line with no #CHROM header is malformed VCF; stop scanning.
        break;
    }

    Ok(Stream {
        reader,
        sample,
        head: None,
    })
}

/// Read forward until the next valid data record (skipping genotype-missing and all-filtered
/// lines) or EOF. Returns the record plus the number of alleles dropped by `--min-support` on the
/// way (dropped counts from skipped lines are folded into the next returned record).
fn next_record(
    reader: &mut Box<dyn BufRead>,
    min_support: Option<u32>,
    contig_rank: &mut HashMap<Box<str>, u32>,
    next_rank: &mut u32,
) -> Result<(Option<Record>, usize), std::io::Error> {
    let mut line = String::new();
    let mut dropped_total = 0usize;

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            return Ok((None, dropped_total));
        }

        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        if fields.len() < 10 {
            continue;
        }

        let chrom = fields[0];
        let pos: u32 = match fields[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let ref_seq = fields[3];
        let alt_field = fields[4];
        let info = fields[7];
        let format = fields[8];
        let sample_field = fields[9];

        let gt_idx = format.split(':').position(|f| f == "GT").unwrap_or(0);
        let gt_field = sample_field.split(':').nth(gt_idx).unwrap_or(".");
        if gt_field.contains('.') {
            continue;
        }
        let allele_indices = parse_gt(gt_field);

        let alt_seqs: Vec<&str> = if alt_field == "." {
            Vec::new()
        } else {
            alt_field.split(',').collect()
        };

        // STRdust SUP is per-allele in GT order; parse lazily only when filtering is on.
        let sup_values: Vec<u32> = if min_support.is_some() {
            parse_sup(format, sample_field)
        } else {
            Vec::new()
        };

        let mut seqs: Vec<Box<str>> = Vec::with_capacity(allele_indices.len());
        for (i, idx) in allele_indices.iter().enumerate() {
            if let Some(min) = min_support {
                // Missing/unparseable SUP for this allele counts as 0 support → dropped. This is
                // intentional: if the user asked for a quality filter and the file can't supply
                // the metric, we err on the side of caution rather than silently passing.
                let sup = sup_values.get(i).copied().unwrap_or(0);
                if sup < min {
                    dropped_total += 1;
                    continue;
                }
            }
            let seq: Box<str> = if *idx == 0 {
                Box::from(ref_seq)
            } else {
                match alt_seqs.get(idx - 1) {
                    Some(&s) => Box::from(s),
                    None => continue,
                }
            };
            seqs.push(seq);
        }

        if seqs.is_empty() {
            // Whole line filtered out (or no resolvable alleles) — keep scanning.
            continue;
        }

        let end = parse_end(info, pos, ref_seq.len() as u32);
        let rank = *contig_rank.entry(Box::from(chrom)).or_insert_with(|| {
            let r = *next_rank;
            *next_rank += 1;
            r
        });

        return Ok((
            Some(Record {
                rank,
                chrom: Box::from(chrom),
                pos,
                end,
                ref_seq: Box::from(ref_seq),
                seqs,
            }),
            dropped_total,
        ));
    }
}

fn parse_contig_id(line: &str) -> Option<&str> {
    // ##contig=<ID=chr1,length=...>
    let rest = line.strip_prefix("##contig=<")?;
    let id_start = rest.find("ID=")? + 3;
    let after = &rest[id_start..];
    let end = after.find([',', '>']).unwrap_or(after.len());
    Some(after[..end].trim())
}

fn parse_sup(format: &str, sample_field: &str) -> Vec<u32> {
    let sup_idx = match format.split(':').position(|f| f == "SUP") {
        Some(i) => i,
        None => return Vec::new(),
    };
    let sup_field = match sample_field.split(':').nth(sup_idx) {
        Some(s) => s,
        None => return Vec::new(),
    };
    // Negative or non-numeric SUP entries clamp to 0 (the strictest interpretation — they
    // will be dropped under any positive --min-support).
    sup_field
        .split(',')
        .map(|v| v.parse::<u32>().unwrap_or(0))
        .collect()
}

fn parse_end(info: &str, pos: u32, ref_len: u32) -> u32 {
    for field in info.split(';') {
        if let Some(val) = field.strip_prefix("END=")
            && let Ok(end) = val.parse()
        {
            return end;
        }
    }
    pos + ref_len
}

fn parse_gt(gt: &str) -> Vec<usize> {
    gt.split(['|', '/'])
        .filter_map(|a| a.trim().parse().ok())
        .collect()
}
