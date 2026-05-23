use std::collections::BTreeMap;

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
        let encode = |bases: &[usize]| -> usize {
            bases.iter().fold(0, |acc, &b| acc * 4 + b)
        };

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
            name_map[ci] = bases
                .iter()
                .map(|&b| BASES[b] as char)
                .collect();
        }

        KmerTable { compact, n, names: name_map }
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
