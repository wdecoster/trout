use clap::Parser;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

mod features;
mod outlier;
mod plot;
mod vcf;

#[derive(Parser)]
#[command(
    name = "trout",
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
        help = "K-mer length for sequence composition features [default: 2]"
    )]
    kmer: usize,

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
        help = "Number of parallel threads for VCF loading and per-locus processing \
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

struct LocusResult {
    chrom: String,
    pos: u32,
    end: u32,
    label: Option<String>,
    outlier_samples: Vec<String>,
    top_axis_names: Vec<String>,
    n_outliers: usize,
    feat_rows: Vec<(String, usize, Vec<f64>, bool)>,
    plot_data: Vec<plot::LocusPlotData>,
}

fn main() {
    env_logger::init();
    let args = Args::parse();

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

    let repeat_map: Option<RepeatMap> = args.repeat.as_deref().map(parse_repeat_file);

    // Collect all unique k values: global default plus any per-locus k from the repeat file.
    let all_ks: HashSet<usize> = {
        let mut ks = HashSet::new();
        ks.insert(args.kmer);
        if let Some(ref rm) = repeat_map {
            for e in rm.values() {
                ks.insert(e.k);
            }
        }
        ks
    };
    let multi_k = all_ks.len() > 1;

    // Precompute KmerTable and feature names for every k that will be needed.
    let kmer_tables: HashMap<usize, features::KmerTable> =
        all_ks.iter().map(|&k| (k, features::KmerTable::new(k))).collect();
    let kmer_names: HashMap<usize, Vec<String>> = kmer_tables
        .iter()
        .map(|(&k, t)| (k, features::feature_names(t)))
        .collect();
    // Names for the global k — used for the features-out header in single-k mode.
    let global_names = kmer_names.get(&args.kmer).unwrap();

    eprintln!(
        "Loading {} VCFs (k={}, eps={}, min_samples={}, length_weight={})",
        n_vcfs, args.kmer, args.eps, min_samples, args.length_weight
    );

    let locus_map = vcf::read_vcfs(&args.vcfs);
    eprintln!("Loaded {} loci", locus_map.len());

    // Optional feature matrix writer — opened before the parallel section.
    let mut feat_writer: Option<BufWriter<Box<dyn Write>>> = args.features_out.as_ref().map(|p| {
        if multi_k {
            eprintln!(
                "Note: --features-out with multiple k values (from --repeat): \
                 k-mer frequency columns will be omitted (length only)."
            );
        }
        let f: Box<dyn Write> =
            Box::new(std::fs::File::create(p).expect("Cannot create features-out file"));
        let mut w = BufWriter::new(f);
        if repeat_map.is_some() {
            write!(w, "chrom\tstart\tend\tname\tsample\tis_outlier\tlength").unwrap();
        } else {
            write!(w, "chrom\tstart\tend\tsample\tis_outlier\tlength").unwrap();
        }
        if !multi_k {
            for name in global_names.iter().skip(1) {
                write!(w, "\t{}", name).unwrap();
            }
        }
        writeln!(w).unwrap();
        w
    });
    let need_feat_rows = feat_writer.is_some();
    let collect_plots = args.plot.is_some();

    // Copy scalars so the parallel closure doesn't need to borrow Args.
    let kmer = args.kmer;
    let eps = args.eps;
    let length_weight = args.length_weight;
    let min_axis_dev = args.min_axis_dev;
    let min_locus_samples = args.min_locus_samples;
    let min_length = args.min_length;
    let samples_of_interest: Option<HashSet<String>> =
        args.samples.as_deref().map(parse_samples);
    let has_repeat = repeat_map.is_some();

    if has_repeat {
        println!("chrom\tstart\tend\tname\tsamples\taxes");
    } else {
        println!("chrom\tstart\tend\tsamples\taxes");
    }

    let mut loci: Vec<_> = locus_map.into_iter().collect();
    loci.sort_unstable_by(|a, b| a.0.cmp(&b.0));

    // Each locus is fully independent: feature computation, DBSCAN, and plot data
    // generation are all embarrassingly parallel. Results are collected in sorted
    // order (rayon preserves input order for into_par_iter + collect).
    let results: Vec<LocusResult> = loci
        .into_par_iter()
        .filter_map(|((chrom, pos, end), alleles)| {
            let unique_samples: HashSet<&str> =
                alleles.iter().map(|a| a.sample.as_str()).collect();
            if unique_samples.len() < min_locus_samples {
                return None;
            }

            let max_len = alleles.iter().map(|a| a.seq.len()).max().unwrap_or(1);
            if max_len == 0 {
                return None;
            }

            // Per-locus k, table, and names from --repeat; fall back to global defaults.
            let repeat_entry =
                repeat_map.as_ref().and_then(|rm| rm.get(&(chrom.clone(), pos, end)));
            let locus_k = repeat_entry.map(|e| e.k).unwrap_or(kmer);
            let locus_table = kmer_tables.get(&locus_k).unwrap();
            let locus_names = kmer_names.get(&locus_k).unwrap();
            let label = repeat_entry.map(|e| e.name.clone());

            let points: Vec<Vec<f64>> = alleles
                .iter()
                .map(|a| {
                    features::build_feature_vector(
                        &a.seq,
                        locus_table,
                        locus_k,
                        max_len,
                        length_weight,
                    )
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
            let locus_top_axes =
                top_axes(&points, &is_noise, &cluster_mean, locus_names, 2, length_weight);

            if let Some(thresh) = min_axis_dev {
                if locus_top_axes
                    .first()
                    .map(|(_, dev)| *dev < thresh)
                    .unwrap_or(true)
                {
                    is_noise.iter_mut().for_each(|n| *n = false);
                }
            }

            let top_axis_names: Vec<String> =
                locus_top_axes.into_iter().map(|(n, _)| n).collect();

            // Collect outlier samples, restricted to --samples list when provided.
            let mut seen: HashSet<&str> = HashSet::new();
            let mut outlier_samples: Vec<String> = Vec::new();
            for (allele, &noise) in alleles.iter().zip(is_noise.iter()) {
                if noise && seen.insert(allele.sample.as_str()) {
                    let in_list = samples_of_interest
                        .as_ref()
                        .map(|s| s.contains(&allele.sample))
                        .unwrap_or(true);
                    if in_list {
                        outlier_samples.push(allele.sample.clone());
                    }
                }
            }
            let n_outliers = outlier_samples.len();

            let feat_rows = if need_feat_rows {
                let mut rows: Vec<(String, usize, Vec<f64>, bool)> = alleles
                    .iter()
                    .zip(points.iter())
                    .zip(is_noise.iter())
                    .map(|((a, pt), &noise)| (a.sample.clone(), a.seq.len(), pt.clone(), noise))
                    .collect();
                rows.sort_by(|(a, _, _, _), (b, _, _, _)| a.cmp(b));
                rows
            } else {
                Vec::new()
            };

            let title_prefix = label.as_deref()
                .map(|s| s.to_string())
                .unwrap_or_else(|| format!("{}:{}-{}", chrom, pos, end));

            let plot_data = if collect_plots && !outlier_samples.is_empty() {
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
                )
            } else {
                Vec::new()
            };

            Some(LocusResult {
                chrom,
                pos,
                end,
                label,
                outlier_samples,
                top_axis_names,
                n_outliers,
                feat_rows,
                plot_data,
            })
        })
        .collect();

    // Sequential output: order is preserved from the sorted loci input.
    let mut total_outliers = 0usize;
    let mut plot_loci: Vec<plot::LocusPlotData> = Vec::new();

    for result in results {
        if !result.outlier_samples.is_empty() {
            if has_repeat {
                let name = result.label.as_deref().unwrap_or(".");
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    result.chrom,
                    result.pos,
                    result.end,
                    name,
                    result.outlier_samples.join(","),
                    result.top_axis_names.join(",")
                );
            } else {
                println!(
                    "{}\t{}\t{}\t{}\t{}",
                    result.chrom,
                    result.pos,
                    result.end,
                    result.outlier_samples.join(","),
                    result.top_axis_names.join(",")
                );
            }
            total_outliers += result.n_outliers;
        }

        if let Some(ref mut w) = feat_writer {
            let name_col = result.label.as_deref().unwrap_or(".");
            for (sample, raw_len, point, noise) in &result.feat_rows {
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
                if !multi_k {
                    for v in point.iter().skip(1) {
                        write!(w, "\t{:.6}", v).unwrap();
                    }
                }
                writeln!(w).unwrap();
            }
        }

        plot_loci.extend(result.plot_data);
    }

    eprintln!("Done: {} outlier calls", total_outliers);

    if let Some(ref path) = args.plot {
        plot::render_scatter_plots(&plot_loci, path);
    }
}

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
) -> Vec<plot::LocusPlotData> {
    let mut pts: Vec<(&str, usize, &Vec<f64>, bool)> = alleles
        .iter()
        .zip(points.iter())
        .zip(is_noise.iter())
        .map(|((a, pt), &noise)| (a.sample.as_str(), a.seq.len(), pt, noise))
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

    let mut plot_data: Vec<plot::LocusPlotData> = Vec::new();

    for axis_name in &unique_axes {
        let y_feat_idx = if axis_name == "length" {
            let group_points: Vec<&Vec<f64>> = pts
                .iter()
                .zip(outlier_top_axes.iter())
                .filter(|((_, _, _, is_out), ax)| {
                    *is_out && ax.as_deref() == Some("length")
                })
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
        let title = format!("{} [{}]", title_prefix, axis_name);

        let mut normal: Vec<(f64, f64, String)> = Vec::new();
        let mut outlier_pts: Vec<(f64, f64, String)> = Vec::new();
        let mut other_outlier_pts: Vec<(f64, f64, String)> = Vec::new();

        for ((sample, raw_len, point, is_out), top_ax) in
            pts.iter().zip(outlier_top_axes.iter())
        {
            let y_val = point[y_feat_idx];
            if *is_out && top_ax.as_deref() == Some(axis_name.as_str()) {
                // Listed outlier for this axis → red
                outlier_pts.push((*raw_len as f64, y_val, sample.to_string()));
            } else if *is_out && top_ax.is_some() {
                // Listed outlier attributed to a different axis → orange (companion plot exists)
                other_outlier_pts.push((*raw_len as f64, y_val, sample.to_string()));
            } else {
                // Normal allele, control outlier, or unattributed outlier → blue
                normal.push((*raw_len as f64, y_val, sample.to_string()));
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
            eprintln!("Warning: could not read repeat file {}: {}", path.display(), e);
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
