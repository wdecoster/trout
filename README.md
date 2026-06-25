# trout

[![test](https://github.com/wdecoster/trout/actions/workflows/test.yml/badge.svg)](https://github.com/wdecoster/trout/actions/workflows/test.yml)

Tool to identify Tandem Repeat OUTliers based on sequence composition and length

## What it does

trout takes a cohort of per-sample VCF files produced by [STRdust](https://github.com/wdecoster/STRdust) and identifies samples that are outliers at any repeat locus — either because their alleles are unusually long, have an unusual sequence composition, or both.

Detection uses multidimensional DBSCAN. Each allele in each sample is represented as a feature vector combining its normalised length and k-mer frequency profile. Samples whose alleles fall outside all clusters (DBSCAN noise points) are reported as outliers. For the details of the feature representation and clustering, see [docs/internals.md](docs/internals.md).

## Installation

Download a prebuilt binary for your platform from the
[releases page](https://github.com/wdecoster/trout/releases) (Linux,
static Linux/musl, and macOS), make it executable, and put it on your `PATH`.

Or build from source with a [Rust toolchain](https://rustup.rs/) (1.85 or newer,
for the 2024 edition):

```bash
cargo install --path .
# or: cargo build --release  → target/release/trout
```

## Usage

```bash
trout [OPTIONS] sample1.vcf.gz sample2.vcf.gz ... > outliers.tsv
# or with a glob
trout [OPTIONS] cohort/*.vcf.gz > outliers.tsv
```

## Output

The main output is tab-separated to stdout, with **one row per outlier sample per locus**:

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome |
| `start` | Locus start (VCF POS, 1-based) |
| `end` | Locus end (from VCF INFO END) |
| `name` | Repeat name from `--repeat` file (only present when `--repeat` is used; `.` if locus has no entry) |
| `sample` | Sample carrying the outlier |
| `allele_length` | Length in bp of the sample's most-deviant flagged allele |
| `top_axis` | Feature axis with the largest deviation from the cluster mean for that allele (e.g. `length` for an expansion, or a k-mer name like `CGG` for a composition outlier) |
| `deviation` | Magnitude of that deviation in the normalized [0,1] feature space (length deviation is rescaled by `--length-weight` for fair comparison with k-mer deviations) |
| `zygosity` | Genotype zygosity at the locus: `hom` (alleles identical), `het` (alleles differ, including compound-het expansions), or `hemi` (a single allele was called, e.g. haploid chrX/Y) |
| `locus_count` | How many distinct samples are outliers at this locus, across any axis (the value is the same on every row of a locus). Lets you tell a one-off from a locus where many samples look unusual — a high `locus_count` often points to a difficult-to-genotype locus rather than per-sample biology. Sort/group by it to rank loci by recurrence. |

Rows are grouped by locus and within a locus sorted by `deviation` descending, so the most extreme calls appear first. A sample with both alleles flagged is reported as a single row — `zygosity` carries the biallelic signal (`hom` for a homozygous expansion, `het` for a compound het, where `allele_length`/`deviation` describe the more-deviant allele). Loci are emitted in genome/contig order (the `##contig` order of the input VCFs). Only loci with at least one outlier are printed.

## Options

### `-k / --kmer` (default: 2)

K-mer length for the sequence composition features. For k=2 there are 16 possible dimers (AA, AC, …, TT); for k=3 there are 64 trimers. Larger k captures more detailed composition information but increases the dimensionality of the DBSCAN feature space, which may require tuning `--eps`. k=2 is a good default for most repeat types.

Pass `-k auto` to let trout pick k per locus from the data, detecting the motif period (2–6 bp) from the repeat sequence itself. Loci with no clean signal (typically a motif longer than 6 bp) fall back to k=3, the most common pathogenic motif length. A histogram of the detected k values is printed to stderr at the end of the run, and an explicit `k` column in `--repeat` always overrides auto detection for that locus. See [docs/internals.md](docs/internals.md#per-locus-k-detection--k-auto) for how the period is detected (REF/median agreement and the base-composition correction).

### `--eps` (default: 0.2)

The DBSCAN neighbourhood radius. A point is considered a cluster member if at least `--min-samples` other points lie within this distance (Euclidean, in the normalised [0,1] feature space). **This is the primary sensitivity control:**

- Smaller eps → tighter clustering → more samples flagged as outliers
- Larger eps → looser clustering → fewer outliers

Typical range is 0.1–0.5. If you are getting too many calls (many loci with outliers), increase eps. If you are missing known expansions, decrease it.

### `--min-samples` (default: log₂(n\_samples) + 1)

Minimum number of nearby points needed to form a cluster. Points without enough neighbours are noise (outliers). The default scales with cohort size: larger cohorts require larger clusters to count as "normal", which keeps the false positive rate roughly stable. Override only if you have a specific reason — for a cohort of 430 samples the default is 9.

### `--length-weight` (default: 1.0)

A multiplier applied to the normalised length dimension before DBSCAN. With the default of 1.0, length and k-mer dimensions are treated equally in Euclidean distance. Increase this if you want pure length outliers (expansions with unchanged composition) to be more detectable relative to composition outliers. A value of 2^k (4 for k=2) gives length equal Euclidean mass to the combined k-mer space, but this can be too aggressive with the default eps and tends to flag natural length variation. Tune this together with `--eps`.

### `--min-locus-samples` (default: 5)

Loci represented in fewer than this many samples are skipped. Increase this if you want a minimum cohort size before calling outliers at a locus.

### `--features-out FILE`

Writes a TSV with the full feature matrix: one row per allele per sample per locus, with columns `chrom`, `start`, `end`, `sample`, `is_outlier`, `length` (raw bp), and one column per k-mer.

When k varies across loci (i.e. `-k auto`, or `--repeat` rows specifying different k values), the per-k-mer columns are replaced by three columns: `k` (the k value used at that locus), `k_source` (`detected` when picked by periodicity, `fallback` when no clean signal was found and the default was used, or `user` when supplied via `-k <int>` or `--repeat`), and `kmer_freqs` (a `;`-separated list of `KMER=FREQ` pairs in canonical-rotation order). Use this if you want to make custom scatter plots or inspect individual feature values.

### `--plot FILE`

Writes an SVG file with one scatter plot per outlier axis per locus. Each plot shows:
- **X axis**: allele length in bp
- **Y axis**: the top composition axis (or the best composition axis for length outliers)
- **Blue points**: normal alleles — shaded from light (one sample) to dark navy (many), see below
- **Red points**: outlier alleles for this axis
- **Orange points**: outlier alleles attributed to a different axis (context)

Outlier sample names are annotated with arrows and are hoverable. A shared legend panel is appended at the end of the grid.

By default the plot is built to **stay openable for large cohorts**: the normal cloud is the overwhelming majority of points and most of them overlap, so it is collapsed to one marker per position, each **shaded by how many samples it holds** (so a position with 300 samples reads as dark navy, not a lone point). Only outliers carry an individual hover label. This cuts the SVG to roughly a tenth the size and ~13× fewer DOM elements than drawing every sample as its own interactive marker — the difference between a file a browser can scroll and one it can't. For full per-point interactivity, see [`--interactive`](#--interactive).

### `--interactive`

Make the `--plot` SVG fully interactive: **hover any point** (normal or outlier) for its sample name, plus a **search box** to filter. This restores per-point detail at the cost of a much larger, heavier file — one DOM node per sample per locus — so it is only practical for **small cohorts or few loci**. Without it (the default) the normal cloud is collapsed and density-shaded as described under `--plot`. Off by default.

### `--jitter`

Add a small deterministic jitter to scatter-plot points so samples sharing the same length and composition — which otherwise stack into a single marker — fan out into a visible cloud. Purely cosmetic: it does not affect outlier calling. Because it relies on every point being drawn individually, it **turns off the default density-collapse** (like `--interactive`), so expect a larger file. Off by default.

### `--min-axis-dev THRESH`

Suppress outlier calls where the top axis deviation from the cluster mean is below THRESH. Values are in the normalized [0,1] feature space: 0.05 corresponds to a 5 percentage point difference in a k-mer frequency or 5% of the maximum allele length. Useful range: 0.05–0.15. Off by default.

This threshold is applied at two levels:
1. **Locus level**: if the average top-axis deviation across all outliers at a locus is below THRESH, the entire locus is suppressed.
2. **Per-outlier level** (plots only): individual outliers whose personal top-axis deviation is below THRESH are shown as orange context in other plots rather than generating their own dedicated axis group.

### `--min-length N`

Only report outliers where the flagged allele is at least N bp long. Useful to suppress noise from short alleles whose length variation is biological rather than pathological. Off by default.

### `--min-fold-length FOLD`

Only report length outliers whose flagged allele differs from the longest-allele cluster's mean length by at least FOLD-fold in either direction — e.g. `1.5` keeps alleles ≥1.5× longer (or ≤0.67× shorter) than that reference. On multi-modal loci the reference is the cluster with the longest alleles, not a pooled mean between modes. Suppresses length outliers with only a modest size change; composition (k-mer) outliers are unaffected. Combine with `--expansions-only` to keep large expansions only. Off by default.

### `--expansions-only`

Only report length outliers that are *longer* than the longest-allele cluster (expansions). A length outlier that is shorter than that reference cluster (a contraction) is suppressed. On multi-modal loci the reference is the cluster with the longest alleles. Composition (k-mer) outliers and length expansions are unaffected. Off by default.

### `--min-support N`

Drop alleles whose STRdust `SUP` (read support) field is below `N` before any clustering. Off by default. Useful for suppressing outlier calls driven by low-support variant calls. `SUP` is read per-allele from the VCF FORMAT (one comma-separated value per called GT allele). Alleles whose `SUP` is missing or unparseable are treated as zero support and dropped — if you set this and your VCFs lack `SUP`, all alleles will be dropped and trout will warn loudly via the "Dropped N alleles" line printed to stderr.

### `--summary FILE`

Write a per-sample QC TSV to FILE with columns `sample`, `in_samples`, `n_loci` (loci where the sample contributed at least one allele), `n_outlier` (loci where the sample was a DBSCAN noise point), `outlier_rate` (n_outlier / n_loci), and `mod_zscore`. Sorted by outlier count desc. Samples flagged at many loci are usually QC issues (low coverage, contamination, swap) rather than biologically interesting — useful as a first pass before investigating individual loci. The QC counts ignore `--samples`, so controls are included alongside the samples-of-interest.

**`mod_zscore`** answers "is this `n_outlier` a lot, or normal?". It is a *robust* z-score of `outlier_rate` across the whole cohort — `(rate − median) / (1.4826 × MAD)`, the Iglewicz–Hoaglin modified z-score, where MAD is the median absolute deviation. The median/MAD scale is used deliberately instead of mean/standard deviation: the high-rate samples we want to flag would otherwise inflate the mean and SD and mask themselves. By the usual convention, **`mod_zscore` > 3.5 marks a sample whose outlier rate is anomalously high** — almost always a technical problem rather than biology. (If more than half the cohort shares one rate so the MAD is zero, the scale falls back to the mean absolute deviation.)

Because `mod_zscore` scores the *rate*, a sample seen at only a handful of loci can still score high off a single flagged locus — read it together with `n_loci`. Samples flagged across many loci with a high `mod_zscore` are the strongest QC candidates.

### `--samples SAMPLES`

Restrict outlier reporting to a subset of samples of interest. Accepts a single sample name, a comma-separated list, or a path to a file with one sample name per line. Samples not in the list are treated as controls: they participate in DBSCAN and appear as blue dots in scatter plots, but are never reported as outliers.

### `--repeat FILE`

A tab-separated file with five columns: `chrom`, `start`, `end`, `name`, `k`. Each row overrides the default behaviour for the matching locus:

- **name** — replaces the coordinate-based label (`chrom:start-end`) in all output columns and scatter-plot titles.
- **k** — uses this k-mer length instead of the global `--kmer` for this locus.

Lines beginning with `#` and lines with fewer than five fields are ignored. Loci not present in the file use the global `--kmer` (or auto-detection, when `-k auto`) and keep their coordinate label. When multiple k values are in use and `--features-out` is also requested, per-k-mer columns are replaced by `k` and `kmer_freqs` columns (see `--features-out`).

### `--threads N`

Number of parallel threads for VCF loading and per-locus processing. Defaults to all available CPUs. Use `--threads 1` to disable parallelism (useful for reproducibility benchmarking or when running many trout processes simultaneously).

## Assumptions

- **Input VCFs are produced by STRdust** and follow its format: one sample per file, genotype in GT field, allele sequences in REF/ALT.
- **Samples are diploid**. Both alleles at each locus are included as separate data points in DBSCAN. A sample is flagged if either allele is an outlier. Both alleles appear as individual points in scatter plots and feature output.
- **Files are bgzip-compressed** (`.vcf.gz`). Plain `.vcf` files also work.
- **All input VCFs should come from the same repeat catalog** (i.e. the same STRdust `--pathogenic` or custom catalog run). Loci are matched by chromosome and position.

## Example

```bash
# Basic run — ~1 second for 430 samples
trout cohort/*.vcf.gz > outliers.tsv

# With scatter plots and feature matrix for follow-up
trout cohort/*.vcf.gz \
  --plot outliers.svg \
  --features-out features.tsv \
  > outliers.tsv

# More sensitive (smaller eps, higher length weight)
trout cohort/*.vcf.gz --eps 0.15 --length-weight 2 > outliers_sensitive.tsv
```

## Further reading

Implementation details — the feature representation, the DBSCAN clustering and
its parameters, per-locus k detection (`-k auto`), and the length-weight
rationale — are documented in [docs/internals.md](docs/internals.md).
