# OxoDam

`OxoDam` is a standalone BAM-based tool to profile transversion misincorporation patterns and summarize potential oxidative-style damage signals.

It computes:
- Global enrichment ratios (`G>T/G>C`, `C>A/C>G`, combined)
- End-specific summaries (`3p`, `5p`)
- Terminal-vs-interior fold change
- Per-position mismatch proportions along reads
- Optional PDF plots and TSV outputs

Unlike mapDamage-based workflows, `OxoDam` works directly from BAM + reference FASTA.

## Requirements

- Python 3.8+
- `pysam`
- R (`Rscript`) only if you want PDF plots

Install `pysam` if needed:

```bash
pip install pysam
```

## Files

- Main script: `OxoDam.py`

## Quick Start

Single sample:

```bash
python3 OxoDam.py \
  /path/to/sample.bam \
  --reference /path/to/reference.fa \
  --sample-label sample1 \
  --plot-log-y \
  --plot-pdf-out sample1_transversion_misincorporation.pdf \
  --pos-tsv-out sample1_pos.tsv
```

Batch mode:

```bash
python3 OxoDam.py \
  --sample-list sample_bams.tsv \
  --reference /path/to/reference.fa \
  --plot-log-y \
  --batch-summary-out batch_summary.tsv \
  --batch-plot-dir batch_plots \
  --batch-pos-dir batch_pos
```

`sample_bams.tsv` format (tab-separated):

```tsv
sample_name	/path/to/sample1.bam
sample2	/path/to/sample2.bam
```

Header row is optional (`sample_name` + `path/bam` is recognized and skipped).

## Output Meaning

### Console summary

For each sample, `OxoDam` prints:

- `Reads processed`: number of reads used after filter checks and read cap (`--max-reads`)
- `BAM decompression threads`: value of `--threads`
- `End binning mode`:
  - `strand-normalized` (default): reverse-mapped reads are flipped for 5p/3p assignment
  - `raw query-orientation`: end assignment follows query orientation directly
- `Overall G>T enrichment ratio [G>T / G>C]`:
  - primary GT-focused enrichment metric
  - `>1` suggests G>T occurs more than control G>C
- `Overall C>A enrichment ratio [C>A / C>G]`:
  - complementary metric for C-based transversions
- `Overall combined enrichment [(G>T + C>A) / (G>C + C>G)]`:
  - pooled oxidative-style signal over control transversions
- `3p/5p end ... enrichment ratio`:
  - same ratios computed per end
- `3p/5p end terminal-vs-interior fold (combined ratio)`:
  - `((positions 1..window combined ratio) / (positions >window combined ratio))`
  - `>1` means terminal enrichment
- `G>T/G mismatch proportion across plotted positions: mean, median, max`:
  - summary of per-position `G>T / G` for plotted range (`--plot-max-pos`)
- `C>A/C mismatch proportion across plotted positions: mean, median, max`:
  - summary of per-position `C>A / C` for plotted range

### Per-position TSV (`--pos-tsv-out`, `--batch-pos-dir`)

Columns:

- `End`: `3p` or `5p`
- `Pos`: distance from end (1-based)
- `G`, `C`: reference-base observations at that position
- `G>T`, `G>C`, `C>A`, `C>G`: mismatch counts
- `gt_over_gc`: `(G>T + pseudocount) / (G>C + pseudocount)`
- `ca_over_cg`: `(C>A + pseudocount) / (C>G + pseudocount)`
- `combined_over_other`: `(G>T + C>A + pseudocount) / (G>C + C>G + pseudocount)`
- `p_gt_over_g`: `G>T / G`
- `p_gc_over_g`: `G>C / G`
- `p_ca_over_c`: `C>A / C`
- `p_cg_over_c`: `C>G / C`

### Batch summary TSV (`--batch-summary-out`)

One row per sample with:

- sample/path
- overall ratios
- mean/median/max mismatch proportions
- 3p and 5p enrichment ratios
- 3p and 5p terminal/interior fold values

### PDF plot (`--plot-pdf-out`, `--batch-plot-dir`)

Two panels (`3p`, `5p`) showing per-position mismatch proportions:

- `G>T / G` (orange)
- `G>C / G` (blue)
- `C>A / C` (green)
- `C>G / C` (magenta)

Defaults to positions `1..30` (`--plot-max-pos 30`).

## Parameters

Run `python3 OxoDam.py --help` for full CLI help.

Core inputs:

- `bam` (positional): input BAM for single-sample mode
- `--sample-list`: 2-column TSV for batch mode
- `--reference` (required): reference FASTA

Output options:

- `--pos-tsv-out`: per-position TSV in single-sample mode
- `--plot-pdf-out`: PDF plot in single-sample mode
- `--batch-summary-out`: summary TSV for batch mode
- `--batch-plot-dir`: write one PDF per sample in batch mode
- `--batch-pos-dir`: write one per-position TSV per sample in batch mode

Read selection / filtering:

- `--max-reads` (default `1000000`):
  - max reads to process per sample
  - `<=0` means use all reads
- `--min-mapq` (default `30`)
- `--min-baseq` (default `30`)
- `--region`:
  - optional region string (`chr:start-end`) for testing/debug

Profiling settings:

- `--window` (default `10`):
  - terminal window size used for terminal vs interior fold
- `--max-pos` (default `70`):
  - maximum distance from read end counted internally
- `--plot-max-pos` (default `30`):
  - maximum distance shown in plot and used for mean/median/max mismatch summaries
- `--pseudocount` (default `0.5`):
  - stabilizes ratios when denominator counts are low/zero

Orientation/end handling:

- `--normalize-ends` (default):
  - recommended for robust 3p/5p interpretation
- `--no-normalize-ends`:
  - raw query-orientation end assignment

Performance:

- `--threads` (default `1`):
  - HTSlib BAM decompression threads

Plot axis options:

- `--plot-log-y`: log10 y-axis for mismatch proportion plots
- `--plot-y-max`: fixed y-axis maximum

## Recommended defaults for screening

For broad sample screening:

- Keep defaults for `--max-reads` and `--plot-max-pos`
- Use `--plot-log-y`
- Use `--normalize-ends` (default)
- Keep filters consistent across all samples (`--min-mapq`, `--min-baseq`)

Example:

```bash
python3 OxoDam.py \
  --sample-list sample_bams.tsv \
  --reference /path/to/reference.fa \
  --plot-log-y \
  --batch-summary-out batch_summary.tsv \
  --batch-plot-dir batch_plots
```

## Notes / Caveats

- `--max-reads` processes reads in BAM traversal order (deterministic, not random).
- Very low-count positions can produce unstable peaks; interpret per-position maxima with denominator counts in mind (`G`, `C` columns in TSV).
- Global enrichment ratios are usually more stable than single-position outliers.
