"""On-target depth summary for a BAM restricted to a capture region BED.

Pipes `samtools depth -a -b <bed> <bam>` and aggregates mean depth plus
percent of target bases at or above each requested threshold. Kept
separate from the Snakemake rule because the threshold aggregation is
easier to get right in Python than in a shell one-liner.
"""

import argparse
import subprocess

import pandas as pd


def summarize(bam_path, target_bed, thresholds):
    proc = subprocess.Popen(
        ["samtools", "depth", "-a", "-b", target_bed, bam_path],
        stdout=subprocess.PIPE, text=True,
    )
    total_bases, depth_sum = 0, 0
    above = {t: 0 for t in thresholds}

    for line in proc.stdout:
        depth = int(line.rstrip("\n").split("\t")[2])
        total_bases += 1
        depth_sum += depth
        for t in thresholds:
            if depth >= t:
                above[t] += 1

    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError("samtools depth failed")

    row = {
        "target_bases": total_bases,
        "mean_target_depth": depth_sum / total_bases if total_bases else 0.0,
    }
    row.update({f"pct_bases_ge_{t}x": (above[t] / total_bases * 100 if total_bases else 0.0) for t in thresholds})
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True)
    p.add_argument("--target-bed", required=True)
    p.add_argument("--thresholds", type=int, nargs="+", default=[10, 20, 30])
    p.add_argument("--sample-id", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    row = {"sample_id": args.sample_id, **summarize(args.bam, args.target_bed, args.thresholds)}
    pd.DataFrame([row]).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
