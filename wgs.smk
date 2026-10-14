def wgs_targets(wildcards):
    targets = []
    for sample in samples["sample_id"]:
        targets.append(f"results/{sample}/variants/{sample}.filtered.vcf.gz")
        targets.append(f"results/{sample}/coverage/{sample}.genome_coverage.txt")
        targets.append(f"results/{sample}/cnv/{sample}.ratio.txt")
        targets.append(f"results/{sample}/cnv/{sample}.segments.txt")
        if config.get("annovar", {}).get("humandb_dir"):
            targets.append(f"results/{sample}/annotation/{sample}.exonic_variant_function")
    return targets


# rule all must stay the first rule in the file (Snakemake's default target
# when none is given on the command line) -- hence defining it before the
# `include`, even though wgs_targets isn't evaluated until the DAG is built,
# by which point common.smk's `samples` table already exists.
rule all:
    input:
        wgs_targets,


include: "common.smk"


rule genome_coverage:
    input:
        bam="results/{sample}/dedup/{sample}.dedup.bam",
        bai="results/{sample}/dedup/{sample}.dedup.bam.bai",
    output:
        "results/{sample}/coverage/{sample}.genome_coverage.txt",
    shell:
        "samtools coverage {input.bam} > {output}"


def control_bam_input(wc):
    # Declared as a rule `input`, not just referenced from `params`, so that
    # when control_bam points at another sample's own dedup BAM (the
    # documented use case) Snakemake actually schedules that BAM -- and its
    # index -- to be built first, instead of racing it under parallel jobs.
    if "control_bam" in samples.columns and pd.notna(samples.loc[wc.sample, "control_bam"]):
        path = samples.loc[wc.sample, "control_bam"]
        return [path, f"{path}.bai"]
    return []


rule call_cnv:
    input:
        bam="results/{sample}/dedup/{sample}.dedup.bam",
        bai="results/{sample}/dedup/{sample}.dedup.bam.bai",
        control=control_bam_input,
    output:
        ratio="results/{sample}/cnv/{sample}.ratio.txt",
        segments="results/{sample}/cnv/{sample}.segments.txt",
    params:
        control=lambda wc: (
            f"--control-bam {samples.loc[wc.sample, 'control_bam']}"
            if "control_bam" in samples.columns and pd.notna(samples.loc[wc.sample, "control_bam"])
            else ""
        ),
        bin_size=config.get("bin_size", 50_000),
        ploidy=config.get("ploidy", 2),
        gain=config.get("gain_threshold_log2", 0.321928),
        loss=config.get("loss_threshold_log2", -0.415037),
        min_bins=config.get("min_segment_bins", 3),
    shell:
        "python {workflow.basedir}/call_cnv.py --bam {input.bam} {params.control} "
        "--bin-size {params.bin_size} --ploidy {params.ploidy} "
        "--gain-threshold {params.gain} --loss-threshold {params.loss} "
        "--min-segment-bins {params.min_bins} "
        "--output-ratio {output.ratio} --output-segments {output.segments}"
