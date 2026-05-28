# trout
Tool to identify Tandem Repeat OUTliers based on sequence composition and length

## What it does

trout takes a cohort of per-sample VCF files produced by [STRdust](https://github.com/wdecoster/STRdust) and identifies samples that are outliers at any repeat locus — either because their alleles are unusually long, have an unusual sequence composition, or both.

Detection uses multidimensional DBSCAN. Each allele in each sample is represented as a feature vector combining its normalised length and k-mer frequency profile. Samples whose alleles fall outside all clusters (DBSCAN noise points) are reported as outliers.

## Usage

```bash
trout [OPTIONS] sample1.vcf.gz sample2.vcf.gz ... > outliers.tsv
# or with a glob
trout [OPTIONS] cohort/*.vcf.gz > outliers.tsv
```

## Output

The main output is tab-separated to stdout, with **one row per outlier allele**:

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome |
| `start` | Locus start (VCF POS, 1-based) |
| `end` | Locus end (from VCF INFO END) |
| `name` | Repeat name from `--repeat` file (only present when `--repeat` is used; `.` if locus has no entry) |
| `sample` | Sample carrying this outlier allele |
| `allele_length` | Length of the outlier allele in bp |
| `top_axis` | Feature axis with the largest deviation from the cluster mean for *this* allele (e.g. `length` for an expansion, or a k-mer name like `CGG` for a composition outlier) |
| `deviation` | Magnitude of that deviation in the normalized [0,1] feature space (length deviation is rescaled by `--length-weight` for fair comparison with k-mer deviations) |

Rows are grouped by locus and within a locus sorted by `deviation` descending, so the most extreme calls appear first. A sample with both alleles flagged produces two rows — itself a useful biallelic signal. Only loci with at least one outlier are printed.

## Options

### `-k / --kmer` (default: 2)

K-mer length for the sequence composition features. For k=2 there are 16 possible dimers (AA, AC, …, TT); for k=3 there are 64 trimers. Larger k captures more detailed composition information but increases the dimensionality of the DBSCAN feature space, which may require tuning `--eps`. k=2 is a good default for most repeat types.

Pass `-k auto` to let trout pick k per locus from the data. For each locus it runs a self-shift periodicity check on the REF allele and on the median-length allele in the cohort and accepts the period (in 2..=6) when they agree. REF is short but accurate; the median-length allele is longer but may carry motif-sequence variation — requiring agreement guards against picking a spurious period from either alone. When one signal is missing (REF too short, or median is too irregular) the available one is used; when both detect periods but disagree, or when neither has a clean signal (typically loci whose motif is longer than 6 bp), k falls back to 3 (the most common pathogenic motif length, with 24 canonical features vs only 10 at k=2). A summary histogram of the detected k values is printed to stderr at the end of the run. An explicit `k` column in `--repeat` always overrides auto detection for that locus.

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

Writes an interactive SVG file with one scatter plot per outlier axis per locus. Each plot shows:
- **X axis**: allele length in bp
- **Y axis**: the top composition axis (or the best composition axis for length outliers)
- **Blue points**: normal alleles
- **Red points**: outlier alleles for this axis
- **Orange points**: outlier alleles attributed to a different axis (context)

Outlier sample names are annotated with arrows. Hover over any point in a browser to see the sample name. A shared legend panel is appended at the end of the grid.

### `--min-axis-dev THRESH`

Suppress outlier calls where the top axis deviation from the cluster mean is below THRESH. Values are in the normalized [0,1] feature space: 0.05 corresponds to a 5 percentage point difference in a k-mer frequency or 5% of the maximum allele length. Useful range: 0.05–0.15. Off by default.

This threshold is applied at two levels:
1. **Locus level**: if the average top-axis deviation across all outliers at a locus is below THRESH, the entire locus is suppressed.
2. **Per-outlier level** (plots only): individual outliers whose personal top-axis deviation is below THRESH are shown as orange context in other plots rather than generating their own dedicated axis group.

### `--min-length N`

Only report outliers where the flagged allele is at least N bp long. Useful to suppress noise from short alleles whose length variation is biological rather than pathological. Off by default.

### `--min-support N`

Drop alleles whose STRdust `SUP` (read support) field is below `N` before any clustering. Off by default. Useful for suppressing outlier calls driven by low-coverage assemblies that produced an unusual sequence by accident. `SUP` is read per-allele from the VCF FORMAT (one comma-separated value per called GT allele). Alleles whose `SUP` is missing or unparseable are treated as zero support and dropped — if you set this and your VCFs lack `SUP`, all alleles will be dropped and trout will warn loudly via the "Dropped N alleles" line printed to stderr.

### `--summary FILE`

Write a per-sample QC TSV to FILE with columns `sample`, `n_loci` (loci where the sample contributed at least one allele), `n_outlier` (loci where the sample was a DBSCAN noise point), and `outlier_rate` (n_outlier / n_loci). Sorted by outlier count desc. Samples flagged at many loci are usually QC issues (low coverage, contamination, swap) rather than biologically interesting — useful as a first pass before investigating individual loci. The QC counts ignore `--samples`, so controls are included alongside the samples-of-interest.

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
