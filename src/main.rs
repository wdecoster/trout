use clap::Parser;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

mod features;
mod outlier;
mod plot;
mod vcf;

#[derive(Parser)]
#[command(
    name = "trout",
    version,
    about = "Tandem Repeat OUTlier identification based on sequence composition and length",
    long_about = "trout reads STRdust VCF files from a cohort and detects outlier samples \
                  at each repeat locus using multidimensional DBSCAN over repeat length \
                  and k-mer composition."
)]
struct Args {
    /// Input VCF files (plain or gzip-compressed)
    #[arg(required = true)]
    vcfs: Vec<PathBuf>,

    #[arg(
        short = 'k',
        long = "kmer",
        default_value = "2",
        hide_default_value = true,
        help = "K-mer length for sequence composition features, or 'auto' to detect per locus \
                from sequence periodicity (REF and median-length allele must agree) [default: 2]"
    )]
    kmer: String,

    #[arg(
        long = "eps",
        default_value = "0.2",
        hide_default_value = true,
        help = "DBSCAN neighborhood radius in normalized [0,1] feature space. \
                Typical range: 0.1–0.5; smaller values flag more outliers [default: 0.2]"
    )]
    eps: f64,

    #[arg(
        long = "min-samples",
        help = "DBSCAN minimum cluster size [default: log2(n_samples)+1]"
    )]
    min_samples: Option<usize>,

    #[arg(
        long = "min-locus-samples",
        default_value = "5",
        hide_default_value = true,
        help = "Skip loci with fewer than this many samples [default: 5]"
    )]
    min_locus_samples: usize,

    #[arg(
        long = "length-weight",
        default_value = "1.0",
        hide_default_value = true,
        help = "Scale factor for the normalized length feature. Higher values make pure \
                length outliers easier to detect; 2^k gives length equal Euclidean mass \
                to the combined kmer space but can be too aggressive with the default eps. \
                Tune together with --eps [default: 1.0]"
    )]
    length_weight: f64,

    #[arg(
        long = "features-out",
        value_name = "FILE",
        help = "Write a full feature matrix (length + k-mer frequencies per sample per locus) \
                to FILE for scatter plotting"
    )]
    features_out: Option<PathBuf>,

    #[arg(
        long = "plot",
        value_name = "FILE",
        help = "Write an interactive SVG scatter plot (length vs. top composition axis) \
                for each outlier locus to FILE. Hover over points to see sample names."
    )]
    plot: Option<PathBuf>,

    #[arg(
        long = "min-axis-dev",
        value_name = "THRESH",
        help = "Suppress outlier calls where the top axis deviation from the cluster mean \
                is below THRESH. Values are in the normalized [0,1] feature space: 0.05 \
                corresponds to a 5 pp difference in a k-mer frequency or 5%% of the maximum \
                allele length. Also applied per-outlier when generating scatter plots: \
                individual outliers whose top-axis deviation is below THRESH are shown as \
                context (orange) in other plots rather than generating their own axis group. \
                Useful range: 0.05–0.15. Off by default."
    )]
    min_axis_dev: Option<f64>,

    #[arg(
        long = "min-length",
        value_name = "N",
        help = "Only report outliers where at least one flagged allele is ≥ N bp. \
                Useful to suppress noise from short normal-variation alleles. Off by default."
    )]
    min_length: Option<usize>,

    #[arg(
        long = "min-fold-length",
        value_name = "FOLD",
        help = "Only report outliers whose flagged allele length differs from the locus cluster \
                mean by at least FOLD-fold in either direction (e.g. 1.5 keeps alleles ≥1.5x \
                longer, or ≤0.67x shorter, than the cluster mean). Suppresses outliers with only \
                a modest length change. Combine with --expansions-only to keep large expansions \
                only. Off by default."
    )]
    min_fold_length: Option<f64>,

    #[arg(
        long = "expansions-only",
        help = "Only report length outliers that are longer than the cohort (expansions). \
                A flagged allele whose dominant deviation is on the length axis but which is \
                shorter than the cluster (a contraction) is suppressed. Composition (k-mer) \
                outliers and length expansions are unaffected. Off by default."
    )]
    expansions_only: bool,

    #[arg(
        long = "jitter",
        help = "Add a small deterministic jitter to scatter-plot points so samples sharing the \
                same length and composition (which otherwise stack into a single marker) fan out \
                into a visible cloud. Purely cosmetic — does not affect outlier calling. Off by \
                default."
    )]
    jitter: bool,

    #[arg(
        long = "summary",
        value_name = "FILE",
        help = "Write a per-sample QC TSV to FILE: one row per sample with the number of loci \
                where the sample contributed data, the number of loci where it was a DBSCAN \
                noise point, and the resulting outlier rate. Sorted by outlier count desc. \
                Samples flagged at many loci are typically QC issues (low coverage, \
                contamination) rather than biologically interesting. The QC counts ignore the \
                `--samples` filter, so controls are included."
    )]
    summary: Option<PathBuf>,

    #[arg(
        long = "min-support",
        value_name = "N",
        help = "Drop alleles whose STRdust SUP (read support) is below N. Off by default — set \
                to suppress outlier calls driven by low-support assemblies. SUP is per-allele \
                in the VCF FORMAT (one value per called GT allele); alleles whose SUP is \
                missing or unparseable are treated as 0 support and dropped."
    )]
    min_support: Option<u32>,

    #[arg(
        long = "samples",
        value_name = "SAMPLES",
        help = "Samples of interest: a sample name, a comma-separated list, or a path to a \
                file with one sample name per line. Only these samples are reported as outliers; \
                all other samples act as controls and contribute to DBSCAN and scatter plots \
                but are never reported."
    )]
    samples: Option<String>,

    #[arg(
        long = "threads",
        value_name = "N",
        help = "Number of parallel threads for per-locus processing \
                [default: all available CPUs]"
    )]
    threads: Option<usize>,

    #[arg(
        long = "repeat",
        value_name = "FILE",
        help = "TSV file with columns (chrom, start, end, name, k). Matching loci are \
                renamed to <name> in output and scatter-plot titles, and processed with \
                the given k instead of --kmer."
    )]
    repeat: Option<PathBuf>,
}

struct RepeatEntry {
    name: String,
    k: usize,
}

type RepeatMap = HashMap<(String, u32, u32), RepeatEntry>;

/// How the user specified `-k`. `Auto` means detect per locus from sequence periodicity.
enum KSpec {
    Fixed(usize),
    Auto,
}

struct OutlierCall {
    sample: String,
    /// Raw allele length in bp of the sample's most-deviant flagged allele (the one this row
    /// represents). For a hom outlier both alleles are identical; for a compound het the longer/
    /// more-deviant of the two flagged alleles is reported and `zygosity` records that it is het.
    allele_length: usize,
    /// Feature axis with the largest deviation from the cluster mean for the represented allele.
    top_axis: String,
    /// Magnitude of that deviation in the normalized [0,1] feature space (length deviation
    /// is divided by length_weight first so it's comparable to k-mer deviations).
    deviation: f64,
    /// Genotype zygosity of the sample at this locus: "hom" (alleles identical), "het" (alleles
    /// differ), or "hemi" (a single allele was called, e.g. haploid chrX/Y).
    zygosity: &'static str,
}

struct LocusResult {
    chrom: String,
    pos: u32,
    end: u32,
    label: Option<String>,
    /// One row per outlier allele (post `--samples` filter), sorted by deviation desc.
    outlier_calls: Vec<OutlierCall>,
    n_outliers: usize,
    /// Each row: (sample, raw_length, feature_vector, is_outlier, k_used, k_source)
    feat_rows: Vec<(String, usize, Vec<f64>, bool, usize, features::KSource)>,
    plot_data: Vec<plot::LocusPlotData>,
    /// k actually used for this locus (for the --kmer auto summary).
    locus_k: usize,
    /// Whether k was detected, fell back to default, or specified by the user.
    k_source: features::KSource,
    /// Sample names present at this locus, deduplicated (used by --summary).
    unique_samples: Vec<String>,
    /// Sample names flagged as DBSCAN noise at this locus, pre `--samples` filter (--summary).
    /// `outlier_samples` above is post-filter and used for the main outliers TSV.
    noise_samples: Vec<String>,
}

fn main() {
    env_logger::init();
    let args = Args::parse();

    // The k-way merge holds every input VCF open at once. Raise our own open-file soft limit to
    // the hard limit so cohorts larger than the default 1024 ulimit just work without the user
    // having to run `ulimit -n` first. Best-effort: any failure leaves the limit untouched.
    raise_fd_limit();

    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to build thread pool");
    }

    let n_vcfs = args.vcfs.len();
    let min_samples = args
        .min_samples
        .unwrap_or_else(|| (n_vcfs as f64).log2() as usize + 1);

    let k_spec = parse_k_spec(&args.kmer);

    let repeat_map: Option<RepeatMap> = args.repeat.as_deref().map(parse_repeat_file);

    // Collect all unique k values that may be needed.
    //  - Fixed: just that k (plus any per-locus k from --repeat)
    //  - Auto: precompute all of AUTO_K_MIN..=AUTO_K_MAX since we don't know upfront which a
    //    locus will pick; the cost is at most ~5 small tables.
    let all_ks: HashSet<usize> = {
        let mut ks = HashSet::new();
        match k_spec {
            KSpec::Fixed(k) => {
                ks.insert(k);
            }
            KSpec::Auto => {
                for k in features::AUTO_K_MIN..=features::AUTO_K_MAX {
                    ks.insert(k);
                }
            }
        }
        if let Some(ref rm) = repeat_map {
            for e in rm.values() {
                ks.insert(e.k);
            }
        }
        ks
    };
    // multi_k drives the features-out column layout (key=value vs one column per k-mer).
    // Auto always implies multi_k since different loci can land on different k.
    let multi_k = matches!(k_spec, KSpec::Auto) || all_ks.len() > 1;

    // Precompute KmerTable and feature names for every k that will be needed.
    let kmer_tables: HashMap<usize, features::KmerTable> = all_ks
        .iter()
        .map(|&k| (k, features::KmerTable::new(k)))
        .collect();
    let kmer_names: HashMap<usize, Vec<String>> = kmer_tables
        .iter()
        .map(|(&k, t)| (k, features::feature_names(t)))
        .collect();
    // Names for the global k — used for the features-out header in single-k mode.
    let global_names_owned = match k_spec {
        KSpec::Fixed(k) => kmer_names.get(&k).cloned(),
        KSpec::Auto => None,
    };

    let k_display = match k_spec {
        KSpec::Fixed(k) => k.to_string(),
        KSpec::Auto => "auto".to_string(),
    };
    eprintln!(
        "Loading {} VCFs (k={}, eps={}, min_samples={}, length_weight={})",
        n_vcfs, k_display, args.eps, min_samples, args.length_weight
    );

    let mut merger = match vcf::VcfMerger::open(&args.vcfs, args.min_support) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    };

    // Optional feature matrix writer — opened before the parallel section.
    let mut feat_writer: Option<BufWriter<Box<dyn Write>>> = args.features_out.as_ref().map(|p| {
        let f: Box<dyn Write> =
            Box::new(std::fs::File::create(p).expect("Cannot create features-out file"));
        let mut w = BufWriter::new(f);
        if repeat_map.is_some() {
            write!(w, "chrom\tstart\tend\tname\tsample\tis_outlier\tlength").unwrap();
        } else {
            write!(w, "chrom\tstart\tend\tsample\tis_outlier\tlength").unwrap();
        }
        if multi_k {
            // Per-locus k can differ, so use a single key=value column instead of one column
            // per canonical k-mer (which would require a union of columns and waste space).
            // k_source distinguishes detected periods from default fallbacks vs user overrides.
            write!(w, "\tk\tk_source\tkmer_freqs").unwrap();
        } else if let Some(ref names) = global_names_owned {
            for name in names.iter().skip(1) {
                write!(w, "\t{}", name).unwrap();
            }
        }
        writeln!(w).unwrap();
        w
    });
    let need_feat_rows = feat_writer.is_some();
    let collect_plots = args.plot.is_some();
    let collect_summary = args.summary.is_some();

    // Copy scalars so the parallel closure doesn't need to borrow Args.
    let eps = args.eps;
    let length_weight = args.length_weight;
    let min_axis_dev = args.min_axis_dev;
    let min_locus_samples = args.min_locus_samples;
    let min_length = args.min_length;
    let min_fold_length = args.min_fold_length;
    let expansions_only = args.expansions_only;
    let jitter = args.jitter;
    let samples_of_interest: Option<HashSet<String>> = args.samples.as_deref().map(parse_samples);

    // Validate --samples against the cohort. A name-format mismatch (e.g. the file lists bare
    // sample IDs but the VCFs carry prefixed/suffixed names) silently yields zero outliers, so
    // report how many matched and stop if any requested sample is absent from the dataset.
    if let Some(ref wanted) = samples_of_interest {
        let cohort: HashSet<String> = merger
            .sample_names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        let mut missing: Vec<&str> = Vec::new();
        for s in wanted {
            if !cohort.contains(s.as_str()) {
                missing.push(s.as_str());
            }
        }
        missing.sort_unstable();
        eprintln!(
            "--samples: {} requested, {} found in cohort of {} samples",
            wanted.len(),
            wanted.len() - missing.len(),
            cohort.len()
        );
        if !missing.is_empty() {
            let preview: Vec<&str> = missing.iter().take(10).copied().collect();
            eprintln!(
                "Error: {} of {} --samples names are not present in the dataset (showing {}): {}{}",
                missing.len(),
                wanted.len(),
                preview.len(),
                preview.join(", "),
                if missing.len() > preview.len() {
                    ", ..."
                } else {
                    ""
                }
            );
            let mut examples: Vec<&str> = cohort.iter().map(|s| s.as_str()).collect();
            examples.sort_unstable();
            examples.truncate(3);
            eprintln!(
                "Cohort sample names look like: {}. Ensure --samples uses the same identifiers.",
                examples.join(", ")
            );
            std::process::exit(1);
        }
    }

    let has_repeat = repeat_map.is_some();

    if has_repeat {
        println!("chrom\tstart\tend\tname\tsample\tallele_length\ttop_axis\tdeviation\tzygosity");
    } else {
        println!("chrom\tstart\tend\tsample\tallele_length\ttop_axis\tdeviation\tzygosity");
    }

    // Each locus is fully independent: feature computation, DBSCAN, and plot-data generation are
    // all embarrassingly parallel. The merger yields loci already in genome order, so we process
    // them in bounded batches: rayon parallelises within a batch (and preserves order), while the
    // batch boundary caps peak memory at one batch of alleles instead of the whole cohort.
    // Count loci dropped for having too few samples, so a too-high --min-locus-samples or a
    // cohort with no shared loci (e.g. disjoint coordinates) surfaces as a number rather than a
    // silent zero-outlier result.
    let skipped_low_sample = std::sync::atomic::AtomicUsize::new(0);

    let process_locus = |((chrom, pos, end), locus_data): vcf::LocusEntry| -> Option<LocusResult> {
        let chrom = String::from(chrom);
        let vcf::LocusData { ref_seq, alleles } = locus_data;

        let unique_samples: HashSet<&str> = alleles.iter().map(|a| a.sample.as_ref()).collect();
        if unique_samples.len() < min_locus_samples {
            skipped_low_sample.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            return None;
        }

        let max_len = alleles.iter().map(|a| a.seq.len()).max().unwrap_or(1);
        if max_len == 0 {
            return None;
        }

        // Per-locus k precedence:
        //   1. --repeat row's explicit k (most authoritative — user-curated)
        //   2. auto detection from REF + median-length allele (-k auto)
        //   3. global -k value
        let repeat_entry = repeat_map
            .as_ref()
            .and_then(|rm| rm.get(&(chrom.clone(), pos, end)));
        let (locus_k, k_source) = if let Some(e) = repeat_entry {
            (e.k, features::KSource::User)
        } else {
            match k_spec {
                KSpec::Fixed(k) => (k, features::KSource::User),
                KSpec::Auto => {
                    let median = median_length_allele(&alleles);
                    features::detect_locus_k(&ref_seq, median)
                }
            }
        };
        let locus_table = kmer_tables.get(&locus_k).unwrap();
        let locus_names = kmer_names.get(&locus_k).unwrap();
        let label = repeat_entry.map(|e| e.name.clone());

        let points: Vec<Vec<f64>> = alleles
            .iter()
            .map(|a| {
                features::build_feature_vector(&a.seq, locus_table, locus_k, max_len, length_weight)
            })
            .collect();

        let mut is_noise = outlier::find_outliers(&points, eps, min_samples);

        // --min-length: clear the noise flag for alleles shorter than the cutoff.
        if let Some(min_len) = min_length {
            for (noise, allele) in is_noise.iter_mut().zip(alleles.iter()) {
                if *noise && allele.seq.len() < min_len {
                    *noise = false;
                }
            }
        }

        let cluster_mean = cluster_mean(&points, &is_noise);

        // --expansions-only: a length outlier shorter than the cohort cluster is a contraction,
        // not an expansion. Clear the noise flag for flagged alleles whose dominant deviation is
        // on the length axis and whose length is below the cluster mean, so they drop out of
        // calls, counts, and plots alike. Composition (k-mer) outliers are unaffected.
        if expansions_only {
            for (noise, point) in is_noise.iter_mut().zip(points.iter()) {
                if *noise
                    && point[0] < cluster_mean[0]
                    && per_point_top_axis(point, &cluster_mean, locus_names, length_weight).0
                        == "length"
                {
                    *noise = false;
                }
            }
        }

        // --min-fold-length: clear the noise flag for flagged alleles whose length is within
        // min_fold-fold of the cluster mean length (in either direction), keeping only outliers
        // with a substantial relative length change. Cluster mean is the raw-bp mean over the
        // non-noise alleles (not the normalized length feature).
        if let Some(min_fold) = min_fold_length {
            let (sum, n) = alleles.iter().zip(is_noise.iter()).fold(
                (0usize, 0usize),
                |(sum, n), (a, &noise)| {
                    if noise {
                        (sum, n)
                    } else {
                        (sum + a.seq.len(), n + 1)
                    }
                },
            );
            if n > 0 {
                let cluster_len_mean = sum as f64 / n as f64;
                for (noise, allele) in is_noise.iter_mut().zip(alleles.iter()) {
                    if *noise {
                        let l = allele.seq.len() as f64;
                        let fold = (l / cluster_len_mean).max(cluster_len_mean / l);
                        if fold < min_fold {
                            *noise = false;
                        }
                    }
                }
            }
        }

        let locus_top_axes = top_axes(
            &points,
            &is_noise,
            &cluster_mean,
            locus_names,
            2,
            length_weight,
        );

        if let Some(thresh) = min_axis_dev
            && locus_top_axes
                .first()
                .map(|(_, dev)| *dev < thresh)
                .unwrap_or(true)
        {
            is_noise.iter_mut().for_each(|n| *n = false);
        }

        let _ = locus_top_axes; // locus-level axes are only used for the suppression above

        // One outlier record per *sample* (not per allele): a homozygous expansion would
        // otherwise emit two byte-identical rows. We group each sample's alleles at the locus,
        // derive zygosity from the full genotype, and report the sample's most-deviant flagged
        // allele. noise_samples (for --summary) is the set of flagged samples and ignores the
        // `--samples` filter so QC reflects the whole cohort.
        struct SampleAgg<'a> {
            /// All of the sample's allele sequences at this locus (flagged or not) — used to
            /// classify zygosity.
            seqs: Vec<&'a str>,
            /// Flagged (DBSCAN-noise) alleles: (length, feature point) for the represented-allele
            /// pick. Empty if the sample has no outlier allele here.
            flagged: Vec<(usize, &'a Vec<f64>)>,
        }
        let mut by_sample: HashMap<&str, SampleAgg> = HashMap::new();
        for ((allele, &noise), point) in alleles.iter().zip(is_noise.iter()).zip(points.iter()) {
            let agg = by_sample
                .entry(allele.sample.as_ref())
                .or_insert_with(|| SampleAgg {
                    seqs: Vec::new(),
                    flagged: Vec::new(),
                });
            agg.seqs.push(&allele.seq);
            if noise {
                agg.flagged.push((allele.seq.len(), point));
            }
        }

        let mut noise_samples: Vec<String> = Vec::new();
        let mut outlier_calls: Vec<OutlierCall> = Vec::new();
        for (sample, agg) in &by_sample {
            if agg.flagged.is_empty() {
                continue;
            }
            if collect_summary {
                noise_samples.push(sample.to_string());
            }
            let in_list = samples_of_interest
                .as_ref()
                .map(|s| s.contains(*sample))
                .unwrap_or(true);
            if !in_list {
                continue;
            }
            // Represent the sample by its most-deviant flagged allele.
            let mut best: Option<(usize, String, f64)> = None;
            for (len, point) in &agg.flagged {
                let (top_axis, deviation) =
                    per_point_top_axis(point, &cluster_mean, locus_names, length_weight);
                if best.as_ref().is_none_or(|(_, _, d)| deviation > *d) {
                    best = Some((*len, top_axis, deviation));
                }
            }
            let (allele_length, top_axis, deviation) = best.unwrap();
            outlier_calls.push(OutlierCall {
                sample: sample.to_string(),
                allele_length,
                top_axis,
                deviation,
                zygosity: zygosity(&agg.seqs),
            });
        }
        // Sort within locus by deviation desc so the most extreme calls appear first; tie
        // break on sample name for determinism (HashMap iteration order is unstable).
        outlier_calls.sort_by(|a, b| {
            b.deviation
                .partial_cmp(&a.deviation)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.sample.cmp(&b.sample))
        });
        let n_outliers = outlier_calls.len();

        let feat_rows = if need_feat_rows {
            let mut rows: Vec<(String, usize, Vec<f64>, bool, usize, features::KSource)> = alleles
                .iter()
                .zip(points.iter())
                .zip(is_noise.iter())
                .map(|((a, pt), &noise)| {
                    (
                        a.sample.to_string(),
                        a.seq.len(),
                        pt.clone(),
                        noise,
                        locus_k,
                        k_source,
                    )
                })
                .collect();
            rows.sort_by(|(a, _, _, _, _, _), (b, _, _, _, _, _)| a.cmp(b));
            rows
        } else {
            Vec::new()
        };

        let title_prefix = label
            .as_deref()
            .map(|s| s.to_string())
            .unwrap_or_else(|| format!("{}:{}-{}", chrom, pos, end));

        let plot_data = if collect_plots && n_outliers > 0 {
            build_plot_data(
                &alleles,
                &points,
                &is_noise,
                &cluster_mean,
                locus_names,
                &title_prefix,
                length_weight,
                min_axis_dev,
                &samples_of_interest,
                jitter,
            )
        } else {
            Vec::new()
        };

        // Only materialised for --summary; otherwise these per-locus sample lists would hold
        // ~one String per sample for every locus in the results vector for no reason.
        let unique_samples_vec: Vec<String> = if collect_summary {
            unique_samples.iter().map(|s| s.to_string()).collect()
        } else {
            Vec::new()
        };

        Some(LocusResult {
            chrom,
            pos,
            end,
            label,
            outlier_calls,
            n_outliers,
            feat_rows,
            plot_data,
            locus_k,
            k_source,
            unique_samples: unique_samples_vec,
            noise_samples,
        })
    };

    // Output is emitted batch by batch; genome order is preserved across batch boundaries because
    // the merger yields loci in order and rayon's collect keeps within-batch order.
    let mut total_outliers = 0usize;
    let mut plot_loci: Vec<plot::LocusPlotData> = Vec::new();
    // (k, source) -> count, so the auto-k summary can show e.g. k=2 detected separately from
    // k=2 fallback. In non-auto mode this is unused.
    let mut k_histogram: BTreeMap<(usize, features::KSource), usize> = BTreeMap::new();
    // sample -> (n_loci_with_data, n_loci_flagged_noise) for --summary.
    let mut sample_stats: HashMap<String, (usize, usize)> = HashMap::new();
    let mut n_loci = 0usize;

    // Loci per batch: bounds peak memory to ~BATCH_SIZE × samples alleles. Large enough to keep
    // all cores busy, small enough that one batch of alleles is cheap to hold.
    const BATCH_SIZE: usize = 4096;
    loop {
        let batch = match merger.next_batch(BATCH_SIZE) {
            Ok(b) => b,
            Err(e) => {
                eprintln!("Error: {}", e);
                std::process::exit(1);
            }
        };
        if batch.is_empty() {
            break;
        }
        n_loci += batch.len();
        let results: Vec<LocusResult> = batch.into_par_iter().filter_map(&process_locus).collect();

        for result in results {
            *k_histogram
                .entry((result.locus_k, result.k_source))
                .or_insert(0) += 1;

            if collect_summary {
                for s in &result.unique_samples {
                    sample_stats.entry(s.clone()).or_insert((0, 0)).0 += 1;
                }
                for s in &result.noise_samples {
                    sample_stats.entry(s.clone()).or_insert((0, 0)).1 += 1;
                }
            }

            if !result.outlier_calls.is_empty() {
                for call in &result.outlier_calls {
                    if has_repeat {
                        let name = result.label.as_deref().unwrap_or(".");
                        println!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}",
                            result.chrom,
                            result.pos,
                            result.end,
                            name,
                            call.sample,
                            call.allele_length,
                            call.top_axis,
                            call.deviation,
                            call.zygosity,
                        );
                    } else {
                        println!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}",
                            result.chrom,
                            result.pos,
                            result.end,
                            call.sample,
                            call.allele_length,
                            call.top_axis,
                            call.deviation,
                            call.zygosity,
                        );
                    }
                }
                total_outliers += result.n_outliers;
            }
            if let Some(ref mut w) = feat_writer {
                let name_col = result.label.as_deref().unwrap_or(".");
                for (sample, raw_len, point, noise, row_k, row_src) in &result.feat_rows {
                    if has_repeat {
                        write!(
                            w,
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            result.chrom,
                            result.pos,
                            result.end,
                            name_col,
                            sample,
                            *noise as u8,
                            raw_len
                        )
                        .unwrap();
                    } else {
                        write!(
                            w,
                            "{}\t{}\t{}\t{}\t{}\t{}",
                            result.chrom, result.pos, result.end, sample, *noise as u8, raw_len
                        )
                        .unwrap();
                    }
                    if multi_k {
                        // Key=value column: the k used at this locus plus the canonical k-mer
                        // frequencies. point[0] is the (weighted) length feature; skip it.
                        let row_names = kmer_names.get(row_k).unwrap();
                        let kv = point
                            .iter()
                            .zip(row_names.iter())
                            .skip(1)
                            .map(|(v, n)| format!("{}={:.3}", n, v))
                            .collect::<Vec<_>>()
                            .join(";");
                        write!(w, "\t{}\t{}\t{}", row_k, row_src.label(), kv).unwrap();
                    } else {
                        for v in point.iter().skip(1) {
                            write!(w, "\t{:.6}", v).unwrap();
                        }
                    }
                    writeln!(w).unwrap();
                }
            }

            plot_loci.extend(result.plot_data);
        }
    }

    if let Some(min) = args.min_support {
        eprintln!(
            "Dropped {} alleles with SUP < {} (--min-support)",
            merger.n_dropped, min
        );
    }
    eprintln!("Processed {} loci", n_loci);
    let skipped = skipped_low_sample.load(std::sync::atomic::Ordering::Relaxed);
    if skipped > 0 {
        eprintln!(
            "Skipped {} loci with fewer than {} samples (--min-locus-samples)",
            skipped, min_locus_samples
        );
    }

    if matches!(k_spec, KSpec::Auto) {
        let parts: Vec<String> = k_histogram
            .iter()
            .map(|((k, src), n)| format!("k={} ({}): {}", k, src.label(), n))
            .collect();
        eprintln!("Auto-k per locus: {}", parts.join(", "));
    }

    eprintln!("Done: {} outlier calls", total_outliers);

    if let Some(ref path) = args.summary {
        write_summary(path, &sample_stats);
    }

    if let Some(ref path) = args.plot {
        plot::render_scatter_plots(&plot_loci, path);
    }
}

fn write_summary(path: &Path, stats: &HashMap<String, (usize, usize)>) {
    let f = std::fs::File::create(path).expect("Cannot create summary file");
    let mut w = BufWriter::new(f);
    writeln!(w, "sample\tn_loci\tn_outlier\toutlier_rate").unwrap();
    let mut rows: Vec<(&String, usize, usize, f64)> = stats
        .iter()
        .map(|(s, (l, o))| {
            let rate = if *l > 0 { *o as f64 / *l as f64 } else { 0.0 };
            (s, *l, *o, rate)
        })
        .collect();
    // Sort by outlier count desc (the QC signal), then rate desc as tiebreaker so a sample
    // flagged 5/10 ranks above one flagged 5/100.
    rows.sort_by(|a, b| {
        b.2.cmp(&a.2)
            .then(b.3.partial_cmp(&a.3).unwrap_or(std::cmp::Ordering::Equal))
            .then(a.0.cmp(b.0))
    });
    for (sample, n_loci, n_outlier, rate) in rows {
        writeln!(w, "{}\t{}\t{}\t{:.4}", sample, n_loci, n_outlier, rate).unwrap();
    }
}

#[allow(clippy::too_many_arguments)]
fn build_plot_data(
    alleles: &[vcf::Allele],
    points: &[Vec<f64>],
    is_noise: &[bool],
    cluster_mean: &[f64],
    names: &[String],
    title_prefix: &str,
    length_weight: f64,
    min_axis_dev: Option<f64>,
    samples_of_interest: &Option<HashSet<String>>,
    jitter: bool,
) -> Vec<plot::LocusPlotData> {
    let mut pts: Vec<(&str, usize, &Vec<f64>, bool)> = alleles
        .iter()
        .zip(points.iter())
        .zip(is_noise.iter())
        .map(|((a, pt), &noise)| (a.sample.as_ref(), a.seq.len(), pt, noise))
        .collect();
    pts.sort_by_key(|(s, _, _, _)| *s);

    let mut unique_axes: Vec<String> = Vec::new();
    let outlier_top_axes: Vec<Option<String>> = pts
        .iter()
        .map(|(sample, _, point, is_out)| {
            if *is_out {
                // Controls (samples not in the --samples list) don't generate axis groups;
                // they appear as blue background in all plots.
                let is_listed = samples_of_interest
                    .as_ref()
                    .map(|s| s.contains(*sample))
                    .unwrap_or(true);
                if !is_listed {
                    return None;
                }
                let (ax, dev) = per_point_top_axis(point, cluster_mean, names, length_weight);
                // Per-outlier axis-dev filter: unattributed outliers go to blue, not orange.
                if min_axis_dev.map(|t| dev < t).unwrap_or(false) {
                    return None;
                }
                if !unique_axes.contains(&ax) {
                    unique_axes.push(ax.clone());
                }
                Some(ax)
            } else {
                None
            }
        })
        .collect();

    // Distinct samples contributing alleles at this locus, shown in each plot title so heavy
    // marker overlap (many samples on the same coordinate) isn't mistaken for missing data.
    let n_samples = alleles
        .iter()
        .map(|a| a.sample.as_ref())
        .collect::<HashSet<&str>>()
        .len();

    // Length-axis (x) span, used to scale jitter. Computed once since x is the same for every plot.
    let x_min = pts.iter().map(|(_, l, _, _)| *l).min().unwrap_or(0) as f64;
    let x_max = pts.iter().map(|(_, l, _, _)| *l).max().unwrap_or(0) as f64;

    let mut plot_data: Vec<plot::LocusPlotData> = Vec::new();

    for axis_name in &unique_axes {
        let y_feat_idx = if axis_name == "length" {
            let group_points: Vec<&Vec<f64>> = pts
                .iter()
                .zip(outlier_top_axes.iter())
                .filter(|((_, _, _, is_out), ax)| *is_out && ax.as_deref() == Some("length"))
                .map(|((_, _, pt, _), _)| *pt)
                .collect();
            // Exclude kmer axes that already have a dedicated plot — avoids
            // the same kmer name appearing as Y-axis in two different plots.
            let excluded: Vec<usize> = unique_axes
                .iter()
                .filter(|ax| *ax != "length")
                .filter_map(|ax| names.iter().position(|n| n == ax))
                .collect();
            best_kmer_y_for_length_group(&group_points, cluster_mean, &excluded)
        } else {
            names.iter().position(|n| n == axis_name).unwrap_or(1)
        };
        let y_label = names[y_feat_idx].clone();
        let title = format!("{} [{}] (n={})", title_prefix, axis_name, n_samples);

        let mut normal: Vec<(f64, f64, String)> = Vec::new();
        let mut outlier_pts: Vec<(f64, f64, String)> = Vec::new();
        let mut other_outlier_pts: Vec<(f64, f64, String)> = Vec::new();

        // Jitter amplitudes (a small fraction of each axis's spread) plus a tally of how many
        // points share each exact (length, feature-value) coordinate. Only points that collide
        // with another are jittered, so unique points — including a lone outlier — keep their true
        // position and --jitter is a no-op on plots with no overlap. x-jitter does most of the
        // de-overlapping (samples sharing a length fan out horizontally); y-jitter spreads a
        // same-composition row vertically.
        let (jx_amp, jy_amp, collisions) = if jitter {
            let mut y_lo = f64::INFINITY;
            let mut y_hi = f64::NEG_INFINITY;
            let mut counts: HashMap<(usize, u64), usize> = HashMap::new();
            for (_, raw_len, p, _) in &pts {
                let y = p[y_feat_idx];
                y_lo = y_lo.min(y);
                y_hi = y_hi.max(y);
                *counts.entry((*raw_len, y.to_bits())).or_insert(0) += 1;
            }
            ((x_max - x_min) * 0.015, (y_hi - y_lo) * 0.04, counts)
        } else {
            (0.0, 0.0, HashMap::new())
        };

        for ((sample, raw_len, point, is_out), top_ax) in pts.iter().zip(outlier_top_axes.iter()) {
            let y_raw = point[y_feat_idx];
            let overlapped = collisions
                .get(&(*raw_len, y_raw.to_bits()))
                .copied()
                .unwrap_or(0)
                > 1;
            let (ox, oy) = if jitter && overlapped {
                let (a, b) = jitter_offsets(sample, *raw_len);
                (a * jx_amp, b * jy_amp)
            } else {
                (0.0, 0.0)
            };
            let x = *raw_len as f64 + ox;
            let y_val = y_raw + oy;
            if *is_out && top_ax.as_deref() == Some(axis_name.as_str()) {
                // Listed outlier for this axis → red
                outlier_pts.push((x, y_val, sample.to_string()));
            } else if *is_out && top_ax.is_some() {
                // Listed outlier attributed to a different axis → orange (companion plot exists)
                other_outlier_pts.push((x, y_val, sample.to_string()));
            } else {
                // Normal allele, control outlier, or unattributed outlier → blue
                normal.push((x, y_val, sample.to_string()));
            }
        }

        plot_data.push(plot::LocusPlotData {
            title,
            y_label,
            normal,
            outliers: outlier_pts,
            other_outliers: other_outlier_pts,
        });
    }

    plot_data
}

/// Deterministic per-point jitter, each component in [-1, 1], seeded by sample name and allele
/// length so a given point always lands in the same place across runs (reproducible plots).
fn jitter_offsets(sample: &str, raw_len: usize) -> (f64, f64) {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    sample.hash(&mut h);
    raw_len.hash(&mut h);
    let v = h.finish();
    let a = (v & 0xffff_ffff) as f64 / u32::MAX as f64 * 2.0 - 1.0;
    let b = ((v >> 32) & 0xffff_ffff) as f64 / u32::MAX as f64 * 2.0 - 1.0;
    (a, b)
}

/// Mean feature vector of non-outlier alleles (the "cluster centre").
fn cluster_mean(points: &[Vec<f64>], is_noise: &[bool]) -> Vec<f64> {
    let n_features = points[0].len();
    let mut sum = vec![0.0f64; n_features];
    let mut count = 0usize;
    for (p, &noise) in points.iter().zip(is_noise.iter()) {
        if !noise {
            for (s, v) in sum.iter_mut().zip(p.iter()) {
                *s += v;
            }
            count += 1;
        }
    }
    if count == 0 {
        return sum;
    }
    sum.iter().map(|s| s / count as f64).collect()
}

/// Return the feature with the highest absolute deviation of one outlier point from the cluster
/// mean, together with that deviation. Length deviation is un-weighted for fair comparison.
fn per_point_top_axis(
    point: &[f64],
    cluster_mean: &[f64],
    names: &[String],
    length_weight: f64,
) -> (String, f64) {
    let (best_idx, best_dev) = point
        .iter()
        .zip(cluster_mean.iter())
        .enumerate()
        .map(|(i, (&pv, &mv))| {
            let delta = if i == 0 {
                (pv - mv).abs() / length_weight
            } else {
                (pv - mv).abs()
            };
            (i, delta)
        })
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or((1, 0.0));
    (names[best_idx].clone(), best_dev)
}

/// For a group of outlier points that share "length" as their top axis, return the kmer index
/// with the highest mean absolute deviation from the cluster mean. Used as the Y axis so the
/// plot shows composition context alongside the length separation on X.
/// `excluded` lists feature indices that already serve as the title axis of another plot —
/// those are skipped so the same kmer name doesn't appear as Y-axis in two different plots.
fn best_kmer_y_for_length_group(
    points: &[&Vec<f64>],
    cluster_mean: &[f64],
    excluded: &[usize],
) -> usize {
    if points.is_empty() {
        return 1;
    }
    let n_feat = points[0].len();
    let mut total = vec![0.0f64; n_feat];
    for p in points {
        for (i, (&pv, &mv)) in p.iter().zip(cluster_mean.iter()).enumerate().skip(1) {
            total[i] += (pv - mv).abs();
        }
    }
    total
        .iter()
        .enumerate()
        .skip(1)
        .filter(|(i, _)| !excluded.contains(i))
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(1)
}

/// Return the top-n features by mean absolute deviation of outlier alleles from the cluster
/// mean, as (name, deviation) pairs. The length feature deviation is divided by length_weight
/// first so axes are compared in the original (unweighted) scale.
fn top_axes(
    points: &[Vec<f64>],
    is_noise: &[bool],
    cluster_mean: &[f64],
    names: &[String],
    n: usize,
    length_weight: f64,
) -> Vec<(String, f64)> {
    let n_features = points[0].len();
    let mut total_dev = vec![0.0f64; n_features];
    let mut count = 0usize;

    for (p, &noise) in points.iter().zip(is_noise.iter()) {
        if noise {
            for (i, (dev, (&pv, &mv))) in total_dev
                .iter_mut()
                .zip(p.iter().zip(cluster_mean.iter()))
                .enumerate()
            {
                let delta = if i == 0 {
                    (pv - mv) / length_weight
                } else {
                    pv - mv
                };
                *dev += delta.abs();
            }
            count += 1;
        }
    }

    if count == 0 {
        return Vec::new();
    }

    let mut indexed: Vec<(usize, f64)> = total_dev
        .iter()
        .enumerate()
        .map(|(i, &d)| (i, d / count as f64))
        .collect();
    indexed.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    indexed
        .iter()
        .take(n)
        .map(|(i, dev)| (names[*i].clone(), *dev))
        .collect()
}

/// Parse `--samples`: if the argument is a readable file path, read one sample name per line;
/// otherwise treat it as a comma-separated list of names.
fn parse_samples(arg: &str) -> HashSet<String> {
    if let Ok(content) = std::fs::read_to_string(arg) {
        return content
            .lines()
            .map(|l| l.trim())
            .filter(|l| !l.is_empty())
            .map(|l| l.to_string())
            .collect();
    }
    arg.split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect()
}

/// Parse `--repeat`: TSV with columns chrom, start, end, name, k.
/// Lines starting with '#' and lines with fewer than 5 fields are skipped.
fn parse_repeat_file(path: &Path) -> RepeatMap {
    let mut map = RepeatMap::new();
    let content = match std::fs::read_to_string(path) {
        Ok(c) => c,
        Err(e) => {
            eprintln!(
                "Warning: could not read repeat file {}: {}",
                path.display(),
                e
            );
            return map;
        }
    };
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            eprintln!("Warning: skipping short line in repeat file: {}", line);
            continue;
        }
        let start = match fields[1].parse::<u32>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let end = match fields[2].parse::<u32>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let k = match fields[4].parse::<usize>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        map.insert(
            (fields[0].to_string(), start, end),
            RepeatEntry {
                name: fields[3].to_string(),
                k,
            },
        );
    }
    map
}

fn parse_k_spec(s: &str) -> KSpec {
    if s.eq_ignore_ascii_case("auto") {
        return KSpec::Auto;
    }
    match s.parse::<usize>() {
        Ok(k) if k > 0 => KSpec::Fixed(k),
        _ => {
            eprintln!(
                "Error: --kmer must be a positive integer or 'auto', got {:?}",
                s
            );
            std::process::exit(2);
        }
    }
}

/// Returns the sequence of the allele whose length is the median across the pooled cohort.
/// Used by `-k auto` to cross-check the period inferred from REF (short but accurate) against a
/// typical-length allele (longer, more positions to test, but may carry sequence variation).
fn median_length_allele(alleles: &[vcf::Allele]) -> &str {
    debug_assert!(!alleles.is_empty());
    let mut by_len: Vec<&vcf::Allele> = alleles.iter().collect();
    by_len.sort_by_key(|a| a.seq.len());
    &by_len[by_len.len() / 2].seq
}

/// Classify a sample's genotype at a locus from its allele sequences: "hom" when all alleles are
/// identical, "het" when they differ, "hemi" when only a single allele was called (haploid chrX/Y
/// or a single-allele genotype). Compares by sequence so a compound het of two equal-length but
/// different motifs is correctly "het".
fn zygosity(seqs: &[&str]) -> &'static str {
    match seqs.split_first() {
        None => "NA",
        Some((_, [])) => "hemi",
        Some((first, rest)) => {
            if rest.iter().all(|s| s == first) {
                "hom"
            } else {
                "het"
            }
        }
    }
}

/// Raise the soft `RLIMIT_NOFILE` to the hard limit (or a sane cap), so the all-files-open merge
/// can handle large cohorts without the user touching `ulimit`. A process may always raise its own
/// soft limit up to the hard limit without privileges. Best-effort; failures are ignored.
#[cfg(unix)]
fn raise_fd_limit() {
    // If the hard limit is "unlimited", target a concrete value: setting the cur limit above the
    // kernel's fs.nr_open (default 1048576 on Linux) is rejected, so don't ask for infinity.
    const TARGET_CAP: libc::rlim_t = 1 << 20;
    unsafe {
        let mut rl = std::mem::MaybeUninit::<libc::rlimit>::uninit();
        if libc::getrlimit(libc::RLIMIT_NOFILE, rl.as_mut_ptr()) != 0 {
            return;
        }
        let mut rl = rl.assume_init();
        let target = if rl.rlim_max == libc::RLIM_INFINITY {
            TARGET_CAP
        } else {
            rl.rlim_max
        };
        if rl.rlim_cur < target {
            rl.rlim_cur = target;
            // Ignore the result: if the kernel rejects it the original (lower) limit still holds
            // and VcfMerger::open will surface a clear EMFILE message if we then run short.
            let _ = libc::setrlimit(libc::RLIMIT_NOFILE, &rl);
        }
    }
}

#[cfg(not(unix))]
fn raise_fd_limit() {}
