# ODam

`ODam` is a standalone BAM-based tool to profile transversion misincorporation patterns and summarize potential oxidative-style damage signals.

It computes:
- Global enrichment ratios (`G>T/G>C`, `C>A/C>G`, combined)
- End-specific summaries (`3p`, `5p`)
- Terminal-vs-interior fold change
- Per-position mismatch proportions along reads
- Optional PDF plots and TSV outputs

Unlike mapDamage-based workflows, `ODam` works directly from BAM + reference FASTA.

## Requirements

- Python 3.8+
- `pysam`
- R (`Rscript`) only if you want PDF plots

Install `pysam` if needed:

```bash
pip install pysam
```

## Files

- Main script: `ODam.py`

## Quick Start

Single sample:

```bash
python3 ODam.py \
  --bam /path/to/sample.bam \
  --reference /path/to/reference.fa \
  --sample-label sample1 \
  --output-dir odam_out \
  --plot-log-y \
```

Batch mode:

```bash
python3 ODam.py \
  --sample-list sample_bams.tsv \
  --output-dir odam_batch_out \
  --plot-log-y \
```

## Sample List Format

`--sample-list` expects a tab-separated text file with either two or three columns.

Recommended format, with one reference FASTA per sample:

```tsv
sample	bam	reference
sample1	/path/to/sample1.bam	/path/to/reference1.fa
sample2	/path/to/sample2.bam	/path/to/reference2.fa
```

Alternative format, when all BAMs use the same reference:

```tsv
sample	bam
sample1	/path/to/sample1.bam
sample2	/path/to/sample2.bam
```

For the two-column format, provide the shared reference on the command line:

```bash
python3 ODam.py \
  --sample-list sample_bams.tsv \
  --reference /path/to/reference.fa \
  --output-dir odam_batch_out
```

Header row is optional. Comment lines starting with `#` and blank lines are ignored.

Important details:

- Columns must be separated by tabs, not spaces.
- Column 1 is the sample label used in output filenames.
- Column 2 is the BAM path.
- Column 3 is optional and is the reference FASTA path.
- If column 3 is omitted, `--reference` must be supplied.
- If column 3 is present, different samples can use different references.

## Parameters

Run `python3 ODam.py --help` for full CLI help.

Core inputs:

- `--bam`: input BAM for single-sample mode
- `--sample-list`: 2- or 3-column TSV for batch mode
- `--reference`:
  - required in single-sample mode
  - optional in batch mode if sample list provides per-row reference (3rd column)
  - can be used as a global fallback for 2-column batch lists
- `--output-dir`:
  - optional base output directory
  - in single mode, if set and explicit outputs are omitted, `ODam` writes default summary/position/PDF files there
  - in batch mode, if set and explicit outputs are omitted, `ODam` writes `ODam_batch_summary.tsv`, `plots/`, and `positions/` under that directory

Output options:

- `--pos-tsv-out`: per-position TSV in single-sample mode
- `--plot-pdf-out`: PDF plot in single-sample mode
- `--summary-tsv-out`: one-row summary TSV in single-sample mode
- `--summary-tsv-append`: append row to `--summary-tsv-out` (header written if file is new)
- `--batch-summary-out`: summary TSV for batch mode
- `--batch-plot-dir`: write one PDF per sample in batch mode
- `--batch-pos-dir`: write one per-position TSV per sample in batch mode
- `--reuse-existing-pos`: if the expected per-position TSV already exists, rebuild the summary/plot from that file instead of rescanning the BAM

Read selection / filtering:

- `--max-reads` (default `1000000`):
  - randomly samples eligible reads without replacement, with every eligible read having an equal probability of selection
  - processes up to the specified number of reads per sample
  - sampling uses two BAM passes: one to count eligible reads and one to process the selected reads
  - `<=0` means use all reads
- `--random-seed` (default `1`):
  - RNG seed for reproducible random sampling
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

## Output Meaning

### Console summary

For each sample, `ODam` prints:

- `Reads processed`: number of reads used after filter checks and read cap (`--max-reads`)
- `Read sampling mode`:
  - random eligible-read sampling when `--max-reads` is positive
  - all eligible reads when `--max-reads <= 0`
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
- `G>C/G mismatch proportion across plotted positions: mean, median, max`:
  - summary of per-position `G>C / G` for plotted range
- `C>A/C mismatch proportion across plotted positions: mean, median, max`:
  - summary of per-position `C>A / C` for plotted range
- `C>G/C mismatch proportion across plotted positions: mean, median, max`:
  - summary of per-position `C>G / C` for plotted range

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
- `excess_gt_rate`: `max(0, (G>T - G>C) / G)`; control-adjusted per-position rate for adding `G>T` damage
- `excess_ca_rate`: `max(0, (C>A - C>G) / C)`; control-adjusted per-position rate for adding `C>A` damage

The excess-rate columns do not use pseudocounts. A negative signal-minus-control estimate is written as zero, because a negative mutation probability cannot be simulated. If no relevant `G` or `C` observations are available at a position, the corresponding value is `NA`. These are operational control-adjusted estimates rather than proof that the excess was caused by oxidation.

### Batch summary TSV (`--batch-summary-out`)

One row per sample with:

- sample/path/reference
- `status`: `ok`, `plot_error`, or `error`
- `error`: error message for failed samples or failed PDF generation
- `reads_processed`
- `missing_reference_contigs`: number of BAM contigs absent from the reference FASTA
- `reads_skipped_missing_reference`: reads skipped because their contig was absent from the reference FASTA
- `g_bases_observed`, `c_bases_observed`, `gc_bases_observed`: number of reference G/C bases contributing to the tracked read-end summaries
- `gt_count`, `gc_count`, `ca_count`, `cg_count`: raw tracked transversion counts
- `oxidative_transversion_count`: `G>T + C>A`
- `control_transversion_count`: `G>C + C>G`
- `total_tracked_transversion_count`: `G>T + C>A + G>C + C>G`
- overall ratios
- mean/median/max mismatch proportions
- 3p and 5p enrichment ratios
- 3p and 5p terminal/interior fold values
  - includes all four per-position proportions: `G>T/G`, `G>C/G`, `C>A/C`, `C>G/C`

In batch mode, the summary header is created at startup and each sample row is appended and flushed as soon as that sample finishes. This allows the table to be inspected while ODam is still running and preserves completed rows if the job is interrupted. One failed sample does not stop the whole run. Failed samples are retained in the summary table with `status=error`. If only PDF plotting fails, the computed metrics are still written with `status=plot_error`.

If a previous batch run created per-position TSVs but failed before writing the batch summary, rerun with `--reuse-existing-pos` and the same `--batch-pos-dir`. ODam will reuse existing position files where available and only scan BAMs for samples whose position files are missing. Reused rows are marked with `status=ok_reused_pos`. Use this only when the existing position files were generated with the same relevant settings (`--max-reads`, `--random-seed`, `--max-pos`, `--min-mapq`, `--min-baseq`, `--normalize-ends`, and related options).

### PDF plot (`--plot-pdf-out`, `--batch-plot-dir`)

Two panels (`3p`, `5p`) showing per-position mismatch proportions:

- `G>T / G` (orange)
- `G>C / G` (blue)
- `C>A / C` (green)
- `C>G / C` (magenta)

Defaults to positions `1..30` (`--plot-max-pos 30`).

## Recommended defaults for screening

For broad sample screening:

- Keep defaults for `--max-reads` and `--plot-max-pos`
- Use `--plot-log-y`
- Use `--normalize-ends` (default)
- Keep filters consistent across all samples (`--min-mapq`, `--min-baseq`)

Example:

```bash
python3 ODam.py \
  --sample-list sample_bams.tsv \
  --reference /path/to/reference.fa \
  --plot-log-y \
  --batch-summary-out batch_summary.tsv \
  --batch-plot-dir batch_plots
```

## Notes / Caveats

- Positive `--max-reads` values randomly sample eligible reads without replacement. Every eligible read has an equal probability of selection. The default `--random-seed 1` makes repeated runs with the same input and settings reproducible.
- Random sampling requires two passes through each BAM and may therefore take longer than processing the first reads encountered.
- Very low-count positions can produce unstable peaks; interpret per-position maxima with denominator counts in mind (`G`, `C` columns in TSV).
- Global enrichment ratios are usually more stable than single-position outliers.
- Reads aligned to BAM contigs absent from the reference FASTA are skipped and counted in the summary diagnostics.
