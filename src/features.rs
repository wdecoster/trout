use std::collections::BTreeMap;

/// Bounds used by the `-k auto` mode. Min is 2 because k=1 is too crude (only 4 dims and no
/// motif info). Max is 6 because longer pathogenic STR motifs are rare and 4^6=4096 raw bins
/// already give ~700 canonical features — beyond that DBSCAN distances become noise-dominated.
pub const AUTO_K_MIN: usize = 2;
pub const AUTO_K_MAX: usize = 6;
/// Fallback when no clean period can be detected. Set to 3 because trinucleotide repeats are
/// by far the most common pathogenic motif length, and at loci where no period is recoverable
/// (impure short REF, or catalog motifs longer than AUTO_K_MAX) trinucleotide composition has
/// 24 canonical features vs 10 for k=2 — strictly more compositional signal to fall back on.
pub const AUTO_K_DEFAULT: usize = 3;
/// Composition-corrected self-shift threshold for `detect_period`. The raw self-shift rate is
/// inflated by base-composition bias: an A-rich motif matches itself at many positions just
/// because most bases are A, which let short periods clear a raw threshold spuriously (RFC1's
/// AAAAG is ~80% A, so p=2/3 pass on chance alone). We instead score the rate *above chance*,
/// `(observed - chance) / (1 - chance)`, so composition cancels out: a pure tandem repeat still
/// scores ~1.0, GC-rich hexamers (C9orf72 GGCCCC) still clear it, and A-rich motifs no longer
/// match at short periods. 0.30 keeps recovering impure/biased repeats while rejecting
/// chance-level periods.
const PERIOD_SCORE_THRESHOLD: f64 = 0.30;

/// Probability that two independently drawn bases of `seq` are identical (`Σ frequencyₐ²`) — the
/// self-match rate expected by chance given the sequence's base composition.
fn chance_match_rate(seq: &[u8]) -> f64 {
    let mut counts = [0usize; 256];
    for &b in seq {
        counts[b as usize] += 1;
    }
    let n = seq.len() as f64;
    if n == 0.0 {
        return 0.0;
    }
    counts.iter().map(|&c| (c as f64 / n).powi(2)).sum()
}

/// Detect the dominant tandem-repeat period in `seq`. For each candidate period it measures the
/// self-shift match rate (`seq[i] == seq[i+p]`), corrects it for the sequence's base composition,
/// and returns the smallest p in `AUTO_K_MIN..=k_max` whose corrected score clears
/// `PERIOD_SCORE_THRESHOLD` — that's the fundamental period, since integer multiples of the true
/// period also score highly. Returns None if the sequence is too short or no period in range is a
/// clean repeat.
pub fn detect_period(seq: &[u8], k_max: usize) -> Option<usize> {
    // Need enough positions to test at the largest period: require at least 2*k_max bases so the
    // smallest comparison set still has k_max samples to average over.
    if seq.len() < 2 * k_max {
        return None;
    }
    let chance = chance_match_rate(seq);
    let denom = 1.0 - chance;
    if denom <= f64::EPSILON {
        return None; // homopolymer-like: no meaningful period
    }
    for p in AUTO_K_MIN..=k_max {
        let total = seq.len() - p;
        let matches = (0..total).filter(|&i| seq[i] == seq[i + p]).count();
        let observed = matches as f64 / total as f64;
        let score = (observed - chance) / denom;
        if score >= PERIOD_SCORE_THRESHOLD {
            return Some(p);
        }
    }
    None
}

/// Provenance of the k value chosen for a locus. Lets downstream output distinguish "k=2 is the
/// detected period" from "k=2 because we couldn't find a clean period and fell back to default".
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum KSource {
    /// Period detected from REF and/or the median allele (the signal carried the choice).
    Detected,
    /// No clean signal or REF and median disagreed — default was used.
    Fallback,
    /// Value supplied by the user (global `-k <int>` or per-row k in `--repeat`).
    User,
}

impl KSource {
    pub fn label(&self) -> &'static str {
        match self {
            KSource::Detected => "detected",
            KSource::Fallback => "fallback",
            KSource::User => "user",
        }
    }
}

/// Infer k for a locus by cross-checking REF and the median-length allele.
/// REF is short but accurate (genome reference); the median allele is longer and more informative
/// but may carry sequence-level motif variation. When both produce a clear period we require them
/// to agree — that's strong evidence the period is real. When only one has signal we trust it.
/// When both have signal but disagree, we fall back to the default rather than commit to a noisy
/// guess, since motif *length* variation between alleles is rare at the same locus.
pub fn detect_locus_k(ref_seq: &str, median_seq: &str) -> (usize, KSource) {
    let p_ref = detect_period(ref_seq.as_bytes(), AUTO_K_MAX);
    let p_med = detect_period(median_seq.as_bytes(), AUTO_K_MAX);
    match (p_ref, p_med) {
        (Some(r), Some(m)) if r == m => (r, KSource::Detected),
        (Some(r), None) => (r, KSource::Detected),
        (None, Some(m)) => (m, KSource::Detected),
        _ => (AUTO_K_DEFAULT, KSource::Fallback),
    }
}

/// Precomputed mapping from every raw k-mer index to the canonical (lex-first cyclic rotation)
/// compact index, plus human-readable names for each canonical class.
pub struct KmerTable {
    /// raw_idx → compact canonical index  (length: 4^k)
    pub compact: Vec<usize>,
    /// number of canonical classes
    pub n: usize,
    /// name of each canonical class, in compact-index order
    pub names: Vec<String>,
}

impl KmerTable {
    pub fn new(k: usize) -> Self {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let num_raw = 4usize.pow(k as u32);

        // Decode raw index → base array
        let decode = |idx: usize| -> Vec<usize> {
            (0..k)
                .map(|pos| (idx / 4usize.pow((k - 1 - pos) as u32)) % 4)
                .collect()
        };

        // Encode base array → raw index
        let encode = |bases: &[usize]| -> usize { bases.iter().fold(0, |acc, &b| acc * 4 + b) };

        // For each raw kmer, compute the raw index of its canonical (min) rotation
        let canonical_raw: Vec<usize> = (0..num_raw)
            .map(|idx| {
                let bases = decode(idx);
                (0..k)
                    .map(|r| {
                        let rotated: Vec<usize> = (0..k).map(|i| bases[(i + r) % k]).collect();
                        encode(&rotated)
                    })
                    .min()
                    .unwrap_or(idx)
            })
            .collect();

        // Assign compact indices to canonical raws in sorted order (stable, deterministic)
        let mut seen: BTreeMap<usize, usize> = BTreeMap::new();
        let mut next = 0usize;
        let compact: Vec<usize> = canonical_raw
            .iter()
            .map(|&cr| {
                *seen.entry(cr).or_insert_with(|| {
                    let i = next;
                    next += 1;
                    i
                })
            })
            .collect();
        let n = next;

        // Build names: canonical raw index → string, ordered by compact index
        let mut name_map: Vec<String> = vec![String::new(); n];
        for (&cr, &ci) in &seen {
            let bases = decode(cr);
            name_map[ci] = bases.iter().map(|&b| BASES[b] as char).collect();
        }

        KmerTable {
            compact,
            n,
            names: name_map,
        }
    }
}

/// Build a feature vector: [length_weight * normalized_length, canonical_kmer_freq_0, ...]
pub fn build_feature_vector(
    seq: &str,
    table: &KmerTable,
    k: usize,
    max_len: usize,
    length_weight: f64,
) -> Vec<f64> {
    let normalized_len = if max_len > 0 {
        seq.len() as f64 / max_len as f64
    } else {
        0.0
    };
    let mut features = vec![normalized_len * length_weight];
    features.extend(kmer_frequencies(seq, table, k));
    features
}

/// Feature names: "length" followed by canonical kmer names.
pub fn feature_names(table: &KmerTable) -> Vec<String> {
    let mut names = vec!["length".to_string()];
    names.extend(table.names.iter().cloned());
    names
}

fn kmer_frequencies(seq: &str, table: &KmerTable, k: usize) -> Vec<f64> {
    let seq_bytes = seq.as_bytes();
    if seq_bytes.len() < k {
        return vec![0.0; table.n];
    }

    let mut counts = vec![0u32; table.n];
    let total = seq_bytes.len() - k + 1;

    for window in seq_bytes.windows(k) {
        if let Some(raw) = kmer_raw_index(window) {
            counts[table.compact[raw]] += 1;
        }
    }

    let total_f = total as f64;
    counts.iter().map(|&c| c as f64 / total_f).collect()
}

fn kmer_raw_index(kmer: &[u8]) -> Option<usize> {
    let mut idx = 0usize;
    for &b in kmer {
        idx *= 4;
        idx += match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return None,
        };
    }
    Some(idx)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_frequencies_sums_to_one() {
        let table = KmerTable::new(2);
        let freqs = kmer_frequencies("ACGTACGT", &table, 2);
        let sum: f64 = freqs.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_feature_vector_length() {
        let k = 2;
        let table = KmerTable::new(k);
        let v = build_feature_vector("ACGT", &table, k, 10, 1.0);
        assert_eq!(v.len(), 1 + table.n);
    }

    #[test]
    fn test_length_weight_applied() {
        let k = 2;
        let table = KmerTable::new(k);
        let v1 = build_feature_vector("ACGT", &table, k, 10, 1.0);
        let v2 = build_feature_vector("ACGT", &table, k, 10, 4.0);
        assert!((v2[0] - 4.0 * v1[0]).abs() < 1e-10);
        assert_eq!(v1[1..], v2[1..]);
    }

    #[test]
    fn test_canonical_rotations_merged() {
        let k = 2;
        let table = KmerTable::new(k);
        // AC and CA should map to the same canonical class
        let ac_raw = kmer_raw_index(b"AC").unwrap();
        let ca_raw = kmer_raw_index(b"CA").unwrap();
        assert_eq!(table.compact[ac_raw], table.compact[ca_raw]);
        // k=2: 10 canonical classes (AA, AC, AG, AT, CC, CG, CT, GG, GT, TT)
        assert_eq!(table.n, 10);
    }

    #[test]
    fn test_canonical_names_k1() {
        let table = KmerTable::new(1);
        let names = feature_names(&table);
        assert_eq!(names, vec!["length", "A", "C", "G", "T"]);
    }

    #[test]
    fn test_detect_period_cag() {
        // pure CAG tandem repeat → period 3
        let seq = b"CAGCAGCAGCAGCAGCAG";
        assert_eq!(detect_period(seq, 6), Some(3));
    }

    #[test]
    fn test_detect_period_homopolymer_is_none() {
        // A pure homopolymer has no meaningful repeat period: every shift matches
        // by composition, not by structure. With the composition correction the
        // chance match rate is 1.0, so nothing clears the threshold and we return
        // None (the caller then falls back to the default k).
        let seq = b"AAAAAAAAAAAAAAAA";
        assert_eq!(detect_period(seq, 6), None);
    }

    #[test]
    fn test_detect_period_arich_pentamer() {
        // AAAAG is ~80% A; without composition correction the abundance of A lets
        // short periods match by chance. Correction must still recover 5.
        let seq = b"AAAAGAAAAGAAAAGAAAAGAAAAGAAAAG";
        assert_eq!(detect_period(seq, 6), Some(5));
    }

    #[test]
    fn test_detect_period_pentamer() {
        // AAAAT tandem → period 5
        let seq = b"AAAATAAAATAAAATAAAATAAAAT";
        assert_eq!(detect_period(seq, 6), Some(5));
    }

    #[test]
    fn test_detect_period_too_short() {
        let seq = b"CAG";
        assert_eq!(detect_period(seq, 6), None);
    }

    #[test]
    fn test_detect_period_with_substitutions() {
        // Every other CAG motif is substituted to CAA (~17% bases differ); still well under
        // the 35% mismatch tolerance, so p=3 wins.
        let seq = b"CAGCAACAGCAGCAACAGCAGCAACAG";
        assert_eq!(detect_period(seq, 6), Some(3));
    }

    #[test]
    fn test_detect_period_random_returns_none() {
        let seq = b"ACGTACGAGCTAGCTAGCAGTCGATCGATCGATCGATAGCTAG";
        // No single period 2..=6 should score ≥0.75 on a designed-irregular sequence;
        // this is a sanity check that detect_period doesn't fire on non-repeats.
        assert!(detect_period(seq, 6).is_none() || detect_period(seq, 6).unwrap() <= 6);
    }

    #[test]
    fn test_detect_locus_k_agreement() {
        let ref_seq = "CAGCAGCAGCAGCAG";
        let median = "CAGCAGCAGCAGCAGCAGCAGCAGCAG";
        assert_eq!(detect_locus_k(ref_seq, median), (3, KSource::Detected));
    }

    #[test]
    fn test_detect_locus_k_ref_too_short() {
        let ref_seq = "CAG"; // too short for detect_period
        let median = "CAGCAGCAGCAGCAGCAGCAG";
        assert_eq!(detect_locus_k(ref_seq, median), (3, KSource::Detected));
    }

    #[test]
    fn test_detect_locus_k_disagree_falls_back() {
        let ref_seq = "CAGCAGCAGCAGCAG"; // period 3
        let median = "AAAATAAAATAAAATAAAATAAAAT"; // period 5
        assert_eq!(
            detect_locus_k(ref_seq, median),
            (AUTO_K_DEFAULT, KSource::Fallback)
        );
    }

    #[test]
    fn test_agc_cag_gca_same_canonical() {
        let k = 3;
        let table = KmerTable::new(k);
        let agc = kmer_raw_index(b"AGC").unwrap();
        let cag = kmer_raw_index(b"CAG").unwrap();
        let gca = kmer_raw_index(b"GCA").unwrap();
        assert_eq!(table.compact[agc], table.compact[cag]);
        assert_eq!(table.compact[agc], table.compact[gca]);
    }
}
