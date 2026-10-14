# WGS/WES Analysis Pipeline

Snakemake workflow that takes paired-end FASTQ reads from a whole-genome
or whole-exome sequencing run through alignment, variant calling, and
QC. Two entry points, `wgs.smk` and `wes.smk`, share one rule file
(`common.smk`) for everything up through variant filtering and diverge
only where the biology actually diverges: WGS gets genome-wide coverage
and binned-read-depth CNV calling, WES gets on-target coverage against a
capture BED instead.

```
FASTQ --fastqc--> QC report
FASTQ --fastp--> trimmed FASTQ --bwa mem | samtools sort--> BAM
BAM --gatk MarkDuplicates--> dedup BAM
dedup BAM --gatk HaplotypeCaller--> raw VCF
raw VCF --GATK hard filters (SNP/INDEL separately)--> filtered VCF
filtered VCF --ANNOVAR (optional)--> exonic_variant_function

WGS only: dedup BAM --samtools coverage--> genome coverage
          dedup BAM --call_cnv.py--> ratio + segment files
WES only: dedup BAM --target_coverage_summary.py--> on-target coverage
```

## Getting it running

```bash
conda env create -f environment.yml   # bwa, samtools, fastp, fastqc, gatk4, snakemake
conda activate wgs-wes-analysis-pipeline
```

(`requirements.txt` covers the pure-Python side if you'd rather manage
the bioinformatics tools yourself and just `pip install` the rest.)

Index the reference once:

```bash
python prepare_reference.py --reference reference/genome.fasta
```

Edit `config.yaml` in place -- it's committed with placeholder paths, not
copied from an `.example` file. At minimum set `reference_fasta` and
`samples` (a CSV path). WES runs also need `target_bed`, either as a
top-level path in the config or a per-sample `target_bed` column in the
sample sheet if different samples used different capture kits.

Sample sheet columns:

- always: `sample_id`, `fastq_r1`, `fastq_r2`
- WGS, optional: `control_bam` -- an already-aligned, indexed BAM to use
  as the CNV ratio denominator instead of the sample's own genome-wide
  median
- WES, optional: `target_bed` -- overrides the config's `target_bed` for
  that row

Then run whichever workflow applies:

```bash
snakemake -s wgs.smk --configfile config.yaml --cores 8
snakemake -s wes.smk --configfile config.yaml --cores 8
```

`snakemake -s wgs.smk --configfile config.yaml -n` gives a dry run --
useful for checking the DAG and sample sheet resolve correctly before
committing compute to it. Outputs land under `results/<sample_id>/...`
(fastqc/, trimmed/, aligned/, dedup/, variants/, and coverage/ or cnv/
depending on which workflow you ran); nothing here writes outside
`results/`.

## On the parts that aren't full production pipelines

ANNOVAR isn't bundled -- it's registration-licensed and can't be pulled
in through conda. `rule annotate_variants` in `common.smk` only enters
the DAG if `config["annovar"]["humandb_dir"]` is set, so leaving it
blank just skips annotation rather than failing the run.

Variant filtering uses the standard GATK hard-filter thresholds
(`SNP_FILTER`/`INDEL_FILTER` in `common.smk`) rather than VQSR, since
VQSR needs population-scale truth sets that are organism/build-specific
and don't generalize into a template repo. Same reasoning for BQSR --
not included, add it between `mark_duplicates` and `call_variants` if
you have a known-sites resource for your reference.

`call_cnv.py` is a plain binned read-depth caller: count reads per bin,
ratio against a control BAM or the sample's own median, rolling-median
smoothing per chromosome, merge adjacent bins past a log2 threshold into
segments. No mappability or GC correction, no blacklist. It's enough to
see gains and losses in exploratory work, and the output columns match
Control-FREEC's ratio/segment format so it slots into tooling built
around that, but for anything going into a publication use CNVkit or
Control-FREEC directly. It's also why CNV calling is WGS-only -- exome
capture coverage isn't uniform enough across the genome for a bin-count
approach to mean anything.

## License

MIT -- see `LICENSE`.
