# trout internals

This document covers how trout represents alleles and detects outliers. It is
aimed at people who want to understand or tune the algorithm beyond what the
[README](../README.md) options describe. None of this is required to use trout.

## Feature representation

Every called allele of every sample becomes one point in a feature space. A
sample is diploid, so it contributes up to two points per locus. Each feature
vector is:

```
[ length_weight * normalized_length, kmer_freq_0, kmer_freq_1, ... ]
```

- **Normalized length** — the allele length in bp divided by the longest allele
  at that locus, so it lands in `[0, 1]` and is comparable across loci of very
  different sizes. It is then multiplied by `--length-weight` (see below).
- **k-mer frequencies** — the fraction of k-mer windows of each *canonical
  class*. These already sum to 1, so they live in `[0, 1]` too.

Putting length and composition on the same `[0, 1]` scale means a single
Euclidean distance in this space mixes "how long" and "what sequence" without
either dimension dominating purely because of its units. The `deviation` column
in the output is a distance in exactly this space.

### Canonical k-mer classes

k-mers that are cyclic rotations of each other describe the same tandem-repeat
motif viewed at a different phase — `CAG`, `AGC` and `GCA` are the same repeat.
trout collapses each k-mer to the lexicographically smallest rotation
(`KmerTable` in `src/features.rs`) so phase does not split one motif's signal
across several features. This is why the feature count is smaller than `4^k`:

| k | raw k-mers (`4^k`) | canonical classes |
|---|--------------------|-------------------|
| 2 | 16                 | 10                |
| 3 | 64                 | 24                |

Reverse complements are **not** collapsed: STRdust alleles are all oriented
relative to the reference strand, so a motif and its reverse complement are
genuinely different observations here, not phase variants of one motif.

## Outlier detection (DBSCAN)

trout runs [DBSCAN](https://en.wikipedia.org/wiki/DBSCAN) over the per-allele
feature vectors at each locus independently.

- A point is a **core point** if at least `--min-samples` points (including
  itself) lie within Euclidean distance `--eps` of it.
- Core points and the points reachable from them form **clusters** — the
  "normal" alleles.
- Points that belong to no cluster are **noise points**, and trout reports the
  samples carrying them as outliers.

`--eps` is therefore the primary sensitivity knob: a smaller radius makes
clusters tighter and flags more alleles; a larger radius is more permissive.
`--min-samples` defaults to `log2(n_samples) + 1` so that larger cohorts need
proportionally larger clusters to count as "normal", which keeps the false
positive rate roughly stable as cohort size grows.

A sample is reported once per locus even if both its alleles are noise; the
`zygosity` column carries the biallelic detail and `deviation`/`allele_length`
describe the more-deviant allele.

### Length weight

`--length-weight` multiplies only the normalized-length dimension before
distances are computed. With the default `1.0`, length carries the same weight
as any one k-mer class. Because there are `C` canonical k-mer classes, the
*combined* compositional dimensions can outweigh length when an allele is a pure
length expansion with unchanged composition. Setting `--length-weight` to the
k-mer-class count (e.g. ~4 for k=2) gives length roughly equal Euclidean mass to
the whole k-mer block — but this is aggressive with the default `--eps` and
tends to flag ordinary length variation, so it should be tuned together with
`--eps`.

## Per-locus k detection (`-k auto`)

With `-k auto`, trout chooses a k-mer length per locus from the sequence itself
(`detect_locus_k` / `detect_period` in `src/features.rs`).

### Self-shift period detection

For a candidate period `p`, the raw self-shift match rate is the fraction of
positions where `seq[i] == seq[i+p]`. A true tandem repeat of period `p` matches
itself almost perfectly at offset `p` (and at every integer multiple of `p`), so
trout scans `p` from 2 upward and takes the smallest `p` that clears the
threshold — that is the fundamental period.

### Base-composition correction

A raw match rate is biased by the sequence's base composition: an A-rich motif
matches itself at many offsets simply because most bases are A. RFC1's `AAAAG`
is ~80% A, so periods 2 and 3 would clear a naive raw threshold by chance alone.
trout corrects for this by scoring the match rate *above chance*:

```
score = (observed - chance) / (1 - chance)
chance = Σ frequency_base²        (probability two random bases of seq match)
```

Composition cancels out: a pure tandem repeat still scores ~1.0, GC-rich
hexamers like C9orf72's `GGGGCC` still clear the threshold, and A-rich motifs no
longer match at spuriously short periods. The acceptance threshold is `0.30`,
which keeps recovering impure/biased repeats while rejecting chance-level
periods. A homopolymer has `chance ≈ 1.0`, so the denominator collapses and no
period is detected — exactly the desired behaviour (there is no meaningful motif
length).

### REF / median agreement

`detect_locus_k` cross-checks two sequences:

- the **REF** allele — short but accurate (it is the genome reference);
- the **median-length allele** in the cohort — longer and more informative, but
  may carry sequence-level motif variation.

The chosen k is:

| REF period | median period | result |
|------------|---------------|--------|
| `r`        | `r` (agree)   | `r`, source `detected` |
| `r`        | none          | `r`, source `detected` |
| none       | `m`           | `m`, source `detected` |
| `r`        | `m` (disagree)| default `3`, source `fallback` |
| none       | none          | default `3`, source `fallback` |

Requiring agreement when both have signal guards against committing to a noisy
guess, since true motif-*length* variation between alleles at one locus is rare.
The fallback is k=3 because trinucleotide repeats are the most common pathogenic
motif length, and k=3 has 24 canonical features versus 10 at k=2 — strictly more
compositional signal to fall back on. Candidate periods are bounded to 2..=6
(`AUTO_K_MIN`/`AUTO_K_MAX`); loci whose catalog motif is longer than 6 bp
therefore land on the fallback. The `k_source` column in `--features-out`
records `detected`/`fallback`/`user` per locus, and a histogram of the chosen k
values is printed to stderr at the end of a run.
