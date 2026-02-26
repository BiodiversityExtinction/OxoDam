#!/usr/bin/env python3
"""Compute transversion misincorporation profiles directly from BAM files.

Outputs:
- concise per-sample summary (same core metrics as oxidative_report.R)
- optional batch summary TSV
- optional per-position TSV
- optional PDF plots (via base R, no extra R packages required)
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
import tempfile
from collections import defaultdict
from statistics import mean, median
from typing import Dict, Iterable, List, Optional, Tuple

import pysam


TransvCounts = Dict[str, int]
EndPosCounts = Dict[str, Dict[int, TransvCounts]]


def _new_bucket() -> TransvCounts:
    return {"G": 0, "C": 0, "G>T": 0, "G>C": 0, "C>A": 0, "C>G": 0}


def parse_sample_list(path: str) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    with open(path, "r", encoding="utf-8") as fh:
        for i, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"{path}:{i}: expected at least 2 tab-separated columns")
            s, p = parts[0].strip(), parts[1].strip()
            if i == 1 and s.lower() in {"sample", "sample_name"} and re.search(r"path|bam", p.lower()):
                continue
            out.append((s, p))
    if not out:
        raise ValueError(f"No samples found in {path}")
    return out


def compute_counts(
    bam_path: str,
    ref_fasta_path: str,
    max_pos: int,
    min_mapq: int,
    min_baseq: int,
    max_reads: int,
    threads: int,
    normalize_ends: bool,
    region: Optional[str] = None,
) -> Tuple[EndPosCounts, int]:
    counts: EndPosCounts = {"5p": defaultdict(_new_bucket), "3p": defaultdict(_new_bucket)}
    reads_used = 0
    with pysam.AlignmentFile(bam_path, "rb", threads=max(1, threads)) as bam, pysam.FastaFile(ref_fasta_path) as ref:
        read_iter = bam.fetch(region=region) if region else bam.fetch(until_eof=True)
        for aln in read_iter:
            if max_reads > 0 and reads_used >= max_reads:
                break
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.is_duplicate or aln.is_qcfail:
                continue
            if aln.mapping_quality < min_mapq:
                continue
            if aln.query_sequence is None or aln.query_qualities is None:
                continue
            read_len = aln.query_length
            if not read_len or read_len <= 0:
                continue
            ref_start = aln.reference_start
            ref_end = aln.reference_end
            if ref_start is None or ref_end is None or ref_end <= ref_start:
                continue
            ref_seq = ref.fetch(aln.reference_name, ref_start, ref_end).upper()
            if not ref_seq:
                continue
            reads_used += 1
            for qpos, rpos in aln.get_aligned_pairs(matches_only=True):
                if qpos is None or rpos is None:
                    continue
                if qpos >= len(aln.query_sequence) or qpos >= len(aln.query_qualities):
                    continue
                if aln.query_qualities[qpos] < min_baseq:
                    continue
                offset = rpos - ref_start
                if offset < 0 or offset >= len(ref_seq):
                    continue
                ref_base = ref_seq[offset]
                if ref_base not in {"G", "C"}:
                    continue
                read_base = aln.query_sequence[qpos].upper()
                if read_base not in {"A", "C", "G", "T"}:
                    continue
                pos5 = qpos + 1
                pos3 = read_len - qpos
                if normalize_ends and aln.is_reverse:
                    # Flip end assignment for reverse-mapped reads so all reads are
                    # represented in a strand-normalized molecule-end frame.
                    pos_5p = pos3
                    pos_3p = pos5
                else:
                    pos_5p = pos5
                    pos_3p = pos3
                if pos_5p <= max_pos:
                    update_bucket(counts["5p"][pos_5p], ref_base, read_base)
                if pos_3p <= max_pos:
                    update_bucket(counts["3p"][pos_3p], ref_base, read_base)
    return counts, reads_used


def update_bucket(bucket: TransvCounts, ref_base: str, read_base: str) -> None:
    if ref_base == "G":
        bucket["G"] += 1
        if read_base == "T":
            bucket["G>T"] += 1
        elif read_base == "C":
            bucket["G>C"] += 1
    elif ref_base == "C":
        bucket["C"] += 1
        if read_base == "A":
            bucket["C>A"] += 1
        elif read_base == "G":
            bucket["C>G"] += 1


def aggregate_counts(buckets: Iterable[TransvCounts]) -> TransvCounts:
    out = _new_bucket()
    for b in buckets:
        for k in out:
            out[k] += b[k]
    return out


def summarize_bucket(bucket: TransvCounts, pseudocount: float) -> Dict[str, float]:
    gt = float(bucket["G>T"])
    gc = float(bucket["G>C"])
    ca = float(bucket["C>A"])
    cg = float(bucket["C>G"])
    oxidative = gt + ca
    other = gc + cg
    return {
        "gt_over_gc": (gt + pseudocount) / (gc + pseudocount),
        "ca_over_cg": (ca + pseudocount) / (cg + pseudocount),
        "combined": (oxidative + pseudocount) / (other + pseudocount),
    }


def per_position_rows(counts: EndPosCounts, pseudocount: float) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for end in ("3p", "5p"):
        for pos in sorted(counts[end].keys()):
            b = counts[end][pos]
            g = b["G"]
            c = b["C"]
            gt = b["G>T"]
            gc = b["G>C"]
            ca = b["C>A"]
            cg = b["C>G"]
            oxidative = gt + ca
            other = gc + cg
            rows.append(
                {
                    "End": end,
                    "Pos": pos,
                    "G": g,
                    "C": c,
                    "G>T": gt,
                    "G>C": gc,
                    "C>A": ca,
                    "C>G": cg,
                    "gt_over_gc": (gt + pseudocount) / (gc + pseudocount),
                    "ca_over_cg": (ca + pseudocount) / (cg + pseudocount),
                    "combined_over_other": (oxidative + pseudocount) / (other + pseudocount),
                    "p_gt_over_g": (gt / g) if g > 0 else float("nan"),
                    "p_gc_over_g": (gc / g) if g > 0 else float("nan"),
                    "p_ca_over_c": (ca / c) if c > 0 else float("nan"),
                    "p_cg_over_c": (cg / c) if c > 0 else float("nan"),
                }
            )
    return rows


def fold_terminal_vs_interior(counts_by_pos: Dict[int, TransvCounts], window: int, pseudocount: float) -> float:
    term = aggregate_counts(v for p, v in counts_by_pos.items() if p <= window)
    interior = aggregate_counts(v for p, v in counts_by_pos.items() if p > window)
    s_term = summarize_bucket(term, pseudocount)
    s_int = summarize_bucket(interior, pseudocount)
    if s_int["combined"] <= 0:
        return float("nan")
    return s_term["combined"] / s_int["combined"]


def safe_mean(vals: List[float]) -> float:
    vals = [v for v in vals if not math.isnan(v)]
    return mean(vals) if vals else float("nan")


def safe_median(vals: List[float]) -> float:
    vals = [v for v in vals if not math.isnan(v)]
    return median(vals) if vals else float("nan")


def safe_max(vals: List[float]) -> float:
    vals = [v for v in vals if not math.isnan(v)]
    return max(vals) if vals else float("nan")


def write_pos_tsv(path: str, rows: List[Dict[str, object]]) -> None:
    fields = [
        "End",
        "Pos",
        "G",
        "C",
        "G>T",
        "G>C",
        "C>A",
        "C>G",
        "gt_over_gc",
        "ca_over_cg",
        "combined_over_other",
        "p_gt_over_g",
        "p_gc_over_g",
        "p_ca_over_c",
        "p_cg_over_c",
    ]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            out = dict(r)
            for k in ("p_gt_over_g", "p_gc_over_g", "p_ca_over_c", "p_cg_over_c"):
                if isinstance(out[k], float) and math.isnan(out[k]):
                    out[k] = "NA"
            w.writerow(out)


def plot_pdf_from_rows(
    rows: List[Dict[str, object]],
    out_pdf: str,
    sample_label: str,
    plot_max_pos: int,
    plot_log_y: bool,
    plot_y_max: Optional[float],
) -> None:
    # Base-R plotting keeps dependencies light.
    with tempfile.TemporaryDirectory(prefix="bam_tv_plot_") as td:
        tsv = os.path.join(td, "pos.tsv")
        rscript = os.path.join(td, "plot.R")
        write_pos_tsv(tsv, rows)
        y_max_expr = "NA" if plot_y_max is None else f"{plot_y_max:.10g}"
        r_code = f"""\
df <- read.delim("{tsv}", sep="\\t", header=TRUE, check.names=FALSE)
df <- df[df$Pos <= {plot_max_pos}, ]
if (nrow(df) == 0) quit(save='no', status=1)
to_num <- function(x) suppressWarnings(as.numeric(x))
for (k in c("p_gt_over_g","p_gc_over_g","p_ca_over_c","p_cg_over_c","Pos")) df[[k]] <- to_num(df[[k]])
ends <- sort(unique(df$End))
vals <- c(df$p_gt_over_g, df$p_gc_over_g, df$p_ca_over_c, df$p_cg_over_c)
vals <- vals[is.finite(vals)]
y_top <- if (is.na({y_max_expr})) max(vals) * 1.08 else {y_max_expr}
if (!is.finite(y_top) || y_top <= 0) y_top <- 1
y_floor <- min(1e-6, y_top/1e6)
pdf("{out_pdf}", width=10, height=5)
par(mfrow=c(1,length(ends)), mar=c(4.5,4.5,3,1))
for (e in ends) {{
  d <- df[df$End==e, ]
  d <- d[order(d$Pos), ]
  y_gt <- pmax(d$p_gt_over_g, y_floor)
  y_gc <- pmax(d$p_gc_over_g, y_floor)
  y_ca <- pmax(d$p_ca_over_c, y_floor)
  y_cg <- pmax(d$p_cg_over_c, y_floor)
  plot(
    d$Pos, y_gt, type="l", lwd=2, col="#D55E00",
    ylim=if ({str(plot_log_y).upper()}) c(y_floor, y_top) else c(0, y_top),
    log=if ({str(plot_log_y).upper()}) "y" else "",
    yaxt="n",
    xlab=paste0(e, " position from end"),
    ylab=if ({str(plot_log_y).upper()}) "Mismatch proportion at reference base (log10)" else "Mismatch proportion at reference base",
    main=paste0("{sample_label}", " - ", e)
  )
  axis(2, las=2)
  lines(d$Pos, y_gc, lwd=1.5, col="#0072B2")
  lines(d$Pos, y_ca, lwd=1.5, col="#009E73")
  lines(d$Pos, y_cg, lwd=1.5, col="#CC79A7")
  legend(
    "topright",
    legend=c("G>T / G", "G>C / G", "C>A / C", "C>G / C"),
    col=c("#D55E00", "#0072B2", "#009E73", "#CC79A7"),
    lwd=c(2,1.5,1.5,1.5), bty="n", cex=0.8
  )
}}
invisible(dev.off())
"""
        with open(rscript, "w", encoding="utf-8") as fh:
            fh.write(r_code)
        subprocess.run(["Rscript", rscript], check=True)


def analyze_sample(
    sample: str,
    bam_path: str,
    ref_fasta_path: str,
    args: argparse.Namespace,
    pos_tsv_out: Optional[str] = None,
    plot_pdf_out: Optional[str] = None,
) -> Dict[str, object]:
    counts = compute_counts(
        bam_path=bam_path,
        ref_fasta_path=ref_fasta_path,
        max_pos=args.max_pos,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        max_reads=args.max_reads,
        threads=args.threads,
        normalize_ends=args.normalize_ends,
        region=args.region,
    )
    counts, reads_used = counts
    total = aggregate_counts(list(counts["3p"].values()) + list(counts["5p"].values()))
    s_overall = summarize_bucket(total, args.pseudocount)
    s3 = summarize_bucket(aggregate_counts(counts["3p"].values()), args.pseudocount)
    s5 = summarize_bucket(aggregate_counts(counts["5p"].values()), args.pseudocount)
    fold3 = fold_terminal_vs_interior(counts["3p"], args.window, args.pseudocount)
    fold5 = fold_terminal_vs_interior(counts["5p"], args.window, args.pseudocount)
    rows = per_position_rows(counts, args.pseudocount)
    rows_plot = [r for r in rows if int(r["Pos"]) <= args.plot_max_pos]
    gt_vals = [float(r["p_gt_over_g"]) for r in rows_plot]
    ca_vals = [float(r["p_ca_over_c"]) for r in rows_plot]
    gt_mean, gt_median, gt_max = safe_mean(gt_vals), safe_median(gt_vals), safe_max(gt_vals)
    ca_mean, ca_median, ca_max = safe_mean(ca_vals), safe_median(ca_vals), safe_max(ca_vals)

    print(f"[{sample}] Sample summary (terminal window: positions 1-{args.window}; pseudocount: {args.pseudocount:g})")
    if args.max_reads > 0:
        print(f"  Reads processed: {reads_used} (max_reads={args.max_reads})")
    else:
        print(f"  Reads processed: {reads_used} (max_reads=all)")
    print(f"  BAM decompression threads: {args.threads}")
    print(f"  End binning mode: {'strand-normalized' if args.normalize_ends else 'raw query-orientation'}")
    print(f"  Plot/report range: positions 1-{args.plot_max_pos} from each read end")
    print("  Ratios > 1 indicate enrichment of putative oxidative transversions over control transversions.")
    print(f"  Overall G>T enrichment ratio [G>T / G>C]: {s_overall['gt_over_gc']:.3f}")
    print(f"  Overall C>A enrichment ratio [C>A / C>G]: {s_overall['ca_over_cg']:.3f}")
    print(f"  Overall combined enrichment [(G>T + C>A) / (G>C + C>G)]: {s_overall['combined']:.3f}")
    print(f"  3p end G>T enrichment ratio [G>T / G>C]: {s3['gt_over_gc']:.3f}")
    print(f"  3p end C>A enrichment ratio [C>A / C>G]: {s3['ca_over_cg']:.3f}")
    print(f"  3p end terminal-vs-interior fold (combined ratio): {fold3:.3f}")
    print(f"  5p end G>T enrichment ratio [G>T / G>C]: {s5['gt_over_gc']:.3f}")
    print(f"  5p end C>A enrichment ratio [C>A / C>G]: {s5['ca_over_cg']:.3f}")
    print(f"  5p end terminal-vs-interior fold (combined ratio): {fold5:.3f}")
    print(
        f"\n  G>T/G mismatch proportion across plotted positions: "
        f"mean={gt_mean:.5f}, median={gt_median:.5f}, max={gt_max:.5f}"
    )
    print(
        f"  C>A/C mismatch proportion across plotted positions: "
        f"mean={ca_mean:.5f}, median={ca_median:.5f}, max={ca_max:.5f}\n"
    )

    if pos_tsv_out:
        write_pos_tsv(pos_tsv_out, rows)
        print(f"Wrote per-position metrics TSV: {pos_tsv_out}")
    if plot_pdf_out:
        plot_pdf_from_rows(
            rows=rows,
            out_pdf=plot_pdf_out,
            sample_label=sample,
            plot_max_pos=args.plot_max_pos,
            plot_log_y=args.plot_log_y,
            plot_y_max=args.plot_y_max,
        )
        print(f"Wrote transversion-misincorporation profile PDF: {plot_pdf_out}")

    return {
        "sample": sample,
        "path": bam_path,
        "overall_gt_over_gc": s_overall["gt_over_gc"],
        "overall_ca_over_cg": s_overall["ca_over_cg"],
        "overall_combined": s_overall["combined"],
        "mean_p_gt_over_g": gt_mean,
        "median_p_gt_over_g": gt_median,
        "max_p_gt_over_g": gt_max,
        "mean_p_ca_over_c": ca_mean,
        "median_p_ca_over_c": ca_median,
        "max_p_ca_over_c": ca_max,
        "3p_gt_over_gc": s3["gt_over_gc"],
        "3p_ca_over_cg": s3["ca_over_cg"],
        "3p_terminal_interior_fold": fold3,
        "5p_gt_over_gc": s5["gt_over_gc"],
        "5p_ca_over_cg": s5["ca_over_cg"],
        "5p_terminal_interior_fold": fold5,
    }


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("bam", nargs="?", help="Input BAM (single-sample mode)")
    ap.add_argument("--reference", required=True, help="Reference FASTA used for mapping")
    ap.add_argument("--sample-label", default="sample", help="Label for single-sample mode")
    ap.add_argument("--sample-list", help="2-column TSV: sample_name<TAB>bam_path")
    ap.add_argument("--batch-summary-out", help="Write batch summary TSV")
    ap.add_argument("--batch-plot-dir", help="Output directory for batch PDFs")
    ap.add_argument("--batch-pos-dir", help="Output directory for batch per-position TSV files")
    ap.add_argument("--pos-tsv-out", help="Write per-position TSV (single-sample mode)")
    ap.add_argument("--plot-pdf-out", help="Write PDF plot (single-sample mode)")
    ap.add_argument("--window", type=int, default=10, help="Terminal window for fold metric (default: 10)")
    ap.add_argument("--max-pos", type=int, default=70, help="Max distance-from-end to count (default: 70)")
    ap.add_argument("--plot-max-pos", type=int, default=30, help="Max distance-from-end to plot/stats (default: 30)")
    ap.add_argument("--pseudocount", type=float, default=0.5, help="Pseudocount for ratio stability")
    ap.add_argument("--min-mapq", type=int, default=30, help="Minimum mapping quality (default: 30)")
    ap.add_argument("--min-baseq", type=int, default=30, help="Minimum base quality (default: 30)")
    ap.add_argument("--max-reads", type=int, default=1000000, help="Maximum reads to process per sample; <=0 means all (default: 1000000)")
    ap.add_argument("--threads", type=int, default=1, help="HTSlib BAM decompression threads (default: 1)")
    ap.add_argument(
        "--normalize-ends",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Strand-normalize 5p/3p end assignment (default: true)",
    )
    ap.add_argument("--plot-log-y", action="store_true", help="Use log10 y-axis in PDF plot")
    ap.add_argument("--plot-y-max", type=float, help="Force y-axis max in PDF plot")
    ap.add_argument("--region", help="Optional region (chr:start-end) for testing")
    args = ap.parse_args()

    if bool(args.bam) == bool(args.sample_list):
        ap.error("Use exactly one mode: positional BAM or --sample-list.")
    return args


def main() -> None:
    args = parse_args()
    if args.sample_list:
        samples = parse_sample_list(args.sample_list)
        if args.batch_plot_dir:
            os.makedirs(args.batch_plot_dir, exist_ok=True)
        if args.batch_pos_dir:
            os.makedirs(args.batch_pos_dir, exist_ok=True)
        seen: Dict[str, int] = defaultdict(int)
        out_rows: List[Dict[str, object]] = []
        for sample, bam_path in samples:
            seen[sample] += 1
            suffix = f"_{seen[sample]}" if seen[sample] > 1 else ""
            safe = re.sub(r"[^A-Za-z0-9._-]", "_", sample) or "sample"
            pdf_out = (
                os.path.join(args.batch_plot_dir, f"{safe}{suffix}_transversion_misincorporation.pdf")
                if args.batch_plot_dir
                else None
            )
            pos_out = os.path.join(args.batch_pos_dir, f"{safe}{suffix}_pos.tsv") if args.batch_pos_dir else None
            out_rows.append(analyze_sample(sample, bam_path, args.reference, args, pos_out, pdf_out))
        if args.batch_summary_out:
            fields = list(out_rows[0].keys())
            with open(args.batch_summary_out, "w", encoding="utf-8", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
                w.writeheader()
                w.writerows(out_rows)
            print(f"Wrote batch summary TSV: {args.batch_summary_out}")
    else:
        analyze_sample(
            sample=args.sample_label,
            bam_path=args.bam,
            ref_fasta_path=args.reference,
            args=args,
            pos_tsv_out=args.pos_tsv_out,
            plot_pdf_out=args.plot_pdf_out,
        )


if __name__ == "__main__":
    main()
