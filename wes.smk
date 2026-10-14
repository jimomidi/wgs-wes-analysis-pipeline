def wes_targets(wildcards):
    targets = []
    for sample in samples["sample_id"]:
        has_row_bed = "target_bed" in samples.columns and pd.notna(samples.loc[sample, "target_bed"])
        if not has_row_bed and not config.get("target_bed"):
            raise ValueError(
                f"No target_bed for sample '{sample}': set a per-sample target_bed column "
                "or a top-level 'target_bed' in config.yaml."
            )
        targets.append(f"results/{sample}/variants/{sample}.filtered.vcf.gz")
        targets.append(f"results/{sample}/coverage/{sample}.target_coverage.csv")
        if config.get("annovar", {}).get("humandb_dir"):
            targets.append(f"results/{sample}/annotation/{sample}.exonic_variant_function")
    return targets


# rule all must stay the first rule in the file (Snakemake's default target
# when none is given on the command line) -- hence defining it before the
# `include`, even though wes_targets isn't evaluated until the DAG is built,
# by which point common.smk's `samples` table already exists.
rule all:
    input:
        wes_targets,


include: "common.smk"


rule target_coverage:
    input:
        bam="results/{sample}/dedup/{sample}.dedup.bam",
        bai="results/{sample}/dedup/{sample}.dedup.bam.bai",
        target_bed=target_bed_input,
    output:
        "results/{sample}/coverage/{sample}.target_coverage.csv",
    params:
        target=target_bed_for,
        thresholds=config.get("coverage_depth_thresholds", [10, 20, 30]),
    shell:
        "python {workflow.basedir}/target_coverage_summary.py --bam {input.bam} "
        "--target-bed {params.target} --thresholds {params.thresholds} "
        "--sample-id {wildcards.sample} --output {output}"
