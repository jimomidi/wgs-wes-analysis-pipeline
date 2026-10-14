"""One-time reference genome preparation.

Builds the BWA index, the samtools FASTA index, and the GATK sequence
dictionary that both `wgs.smk` and `wes.smk` expect to already exist next
to the reference FASTA. Run once per reference before either workflow.

    python prepare_reference.py --reference reference/genome.fasta
"""

import argparse
import subprocess
from pathlib import Path


def run(cmd):
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", type=Path, required=True)
    args = p.parse_args()

    run(["bwa", "index", str(args.reference)])
    run(["samtools", "faidx", str(args.reference)])
    run(["gatk", "CreateSequenceDictionary", "-R", str(args.reference),
         "-O", str(args.reference.with_suffix(".dict"))])


if __name__ == "__main__":
    main()
