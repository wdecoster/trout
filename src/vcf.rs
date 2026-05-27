use rayon::prelude::*;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::io::BufRead;
use std::path::{Path, PathBuf};

pub struct Allele {
    pub sample: String,
    pub seq: String,
}

/// All alleles seen at a locus plus the REF sequence from the VCF.
/// REF is the same across all input VCFs at the same coordinates (genome reference),
/// so the first writer wins on merge.
pub struct LocusData {
    pub ref_seq: String,
    pub alleles: Vec<Allele>,
}

/// (chrom, pos, end) -> locus data
pub type LocusMap = HashMap<(String, u32, u32), LocusData>;

pub fn read_vcfs(paths: &[PathBuf]) -> LocusMap {
    paths
        .par_iter()
        .filter_map(|path| match read_vcf(path) {
            Ok(m) => Some(m),
            Err(e) => {
                eprintln!("Warning: failed to read {}: {}", path.display(), e);
                None
            }
        })
        .reduce(LocusMap::new, |mut acc, local| {
            for (key, data) in local {
                match acc.entry(key) {
                    Entry::Occupied(mut e) => e.get_mut().alleles.extend(data.alleles),
                    Entry::Vacant(e) => {
                        e.insert(data);
                    }
                }
            }
            acc
        })
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

fn read_vcf(path: &Path) -> Result<LocusMap, Box<dyn std::error::Error + Send + Sync>> {
    let rdr = reader(path)?;
    let mut sample_name = String::new();
    let mut locus_map = LocusMap::new();

    for line_result in rdr.lines() {
        let line = line_result?;

        if line.starts_with("##") {
            continue;
        }

        if line.starts_with('#') {
            let fields: Vec<&str> = line.split('\t').collect();
            sample_name = fields.last().unwrap_or(&"unknown").to_string();
            log::debug!("Reading sample: {}", sample_name);
            continue;
        }

        if sample_name.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            continue;
        }

        let chrom = fields[0].to_string();
        let pos: u32 = fields[1].parse()?;
        let ref_seq = fields[3];
        let alt_field = fields[4];
        let info = fields[7];
        let format = fields[8];
        let sample_field = fields[9];

        let alt_seqs: Vec<&str> = if alt_field == "." {
            Vec::new()
        } else {
            alt_field.split(',').collect()
        };

        let end = parse_end(info, pos, ref_seq.len() as u32);

        let gt_idx = format.split(':').position(|f| f == "GT").unwrap_or(0);
        let gt_field = sample_field.split(':').nth(gt_idx).unwrap_or(".");

        if gt_field.contains('.') {
            continue;
        }

        let allele_indices = parse_gt(gt_field);
        let key = (chrom, pos, end);

        let entry = locus_map.entry(key).or_insert_with(|| LocusData {
            ref_seq: ref_seq.to_string(),
            alleles: Vec::new(),
        });

        for idx in allele_indices {
            let seq = if idx == 0 {
                ref_seq.to_string()
            } else {
                match alt_seqs.get(idx - 1) {
                    Some(&s) => s.to_string(),
                    None => continue,
                }
            };

            entry.alleles.push(Allele {
                sample: sample_name.clone(),
                seq,
            });
        }
    }

    Ok(locus_map)
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
