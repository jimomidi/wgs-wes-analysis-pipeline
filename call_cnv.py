"""Binned read-depth CNV calling from a BAM file.

Counts reads in fixed-size genomic bins, converts counts to a ratio
against either a matched control BAM or the sample's own genome-wide
median, smooths per chromosome, and merges consecutive bins crossing the
gain/loss threshold into segments. Output matches the ratio/segment file
format Control-FREEC produces, so it drops into tooling built around that
format.

This is a simple reference implementation (no mappability/GC correction,
no blacklist handling) -- fine for exploring a dataset, not a replacement
for CNVkit/Control-FREEC on anything going into a paper.
"""

import argparse

import numpy as np
import pandas as pd
import pysam


def binned_counts(bam_path, bin_size):
    rows = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for contig, length in zip(bam.references, bam.lengths):
            for start in range(0, length, bin_size):
                end = min(start + bin_size, length)
                rows.append((contig, start, end, bam.count(contig, start, end)))
    return pd.DataFrame(rows, columns=["Chromosome", "Start", "End", "ReadCount"])


def add_ratio_columns(bins, control_bins, ploidy, smoothing_window=5):
    if control_bins is not None:
        # Normalize both tracks by their own median depth before dividing.
        # Without this, any systematic difference in overall sequencing
        # depth between the sample and control BAM (near-universal in
        # practice) shows up as a spurious genome-wide gain or loss instead
        # of cancelling out of the ratio.
        sample_baseline = bins.loc[bins["ReadCount"] > 0, "ReadCount"].median()
        control_baseline = control_bins.loc[control_bins["ReadCount"] > 0, "ReadCount"].median()
        if not sample_baseline or not control_baseline or pd.isna(sample_baseline) or pd.isna(control_baseline):
            raise ValueError("Cannot normalize: sample or control BAM has no covered bins.")

        ctrl = control_bins.set_index(["Chromosome", "Start"])["ReadCount"].replace(0, np.nan)
        bins = bins.set_index(["Chromosome", "Start"])
        bins["Ratio"] = (bins["ReadCount"] / sample_baseline) / (ctrl / control_baseline)
        bins = bins.reset_index()
    else:
        baseline = bins.loc[bins["ReadCount"] > 0, "ReadCount"].median()
        bins["Ratio"] = bins["ReadCount"] / baseline if baseline else 0.0

    bins["MedianRatio"] = bins.groupby("Chromosome")["Ratio"].transform(
        lambda s: s.rolling(smoothing_window, center=True, min_periods=1).median()
    )
    bins["CopyNumber"] = (bins["MedianRatio"] * ploidy).round()
    return bins.dropna(subset=["Ratio"])


def merge_segments(bins, gain_threshold, loss_threshold, min_bins):
    bins = bins.copy()
    # Classify on the rolling-median-smoothed ratio, not the raw per-bin
    # ratio -- smoothing exists specifically to reduce single-bin noise
    # before thresholding, so segment calling has to use it too.
    bins["log2r"] = np.log2(bins["MedianRatio"].clip(lower=1e-6))
    bins["direction"] = "neutral"
    bins.loc[bins["log2r"] >= gain_threshold, "direction"] = "gain"
    bins.loc[bins["log2r"] <= loss_threshold, "direction"] = "loss"

    segments = []
    for chrom, chrom_bins in bins.groupby("Chromosome", sort=False):
        chrom_bins = chrom_bins.sort_values("Start").reset_index(drop=True)
        run_start = 0
        for i in range(1, len(chrom_bins) + 1):
            end_of_run = i == len(chrom_bins) or chrom_bins.loc[i, "direction"] != chrom_bins.loc[run_start, "direction"]
            if not end_of_run:
                continue
            run = chrom_bins.iloc[run_start:i]
            direction = run.iloc[0]["direction"]
            if direction != "neutral" and len(run) >= min_bins:
                segments.append((chrom, int(run["Start"].min()), int(run["End"].max()),
                                  int(run["CopyNumber"].median()), direction))
            run_start = i
    return pd.DataFrame(segments, columns=["Chromosome", "Start", "End", "CopyNumber", "CNV_Type"])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True)
    p.add_argument("--control-bam", default=None)
    p.add_argument("--bin-size", type=int, default=50_000)
    p.add_argument("--ploidy", type=int, default=2)
    p.add_argument("--gain-threshold", type=float, default=0.321928, help="log2 ratio, default = log2(1.25)")
    p.add_argument("--loss-threshold", type=float, default=-0.415037, help="log2 ratio, default = log2(0.75)")
    p.add_argument("--min-segment-bins", type=int, default=3)
    p.add_argument("--output-ratio", required=True)
    p.add_argument("--output-segments", required=True)
    args = p.parse_args()

    bins = binned_counts(args.bam, args.bin_size)
    control_bins = binned_counts(args.control_bam, args.bin_size) if args.control_bam else None
    bins = add_ratio_columns(bins, control_bins, args.ploidy)
    segments = merge_segments(bins, args.gain_threshold, args.loss_threshold, args.min_segment_bins)

    bins[["Chromosome", "Start", "Ratio", "MedianRatio", "CopyNumber"]].to_csv(args.output_ratio, sep="\t", index=False)
    segments.to_csv(args.output_segments, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
