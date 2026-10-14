import pandas as pd

samples = pd.read_csv(config["samples"], dtype=str).set_index("sample_id", drop=False)


def fastq_r1(wc):
    return samples.loc[wc.sample, "fastq_r1"]


def fastq_r2(wc):
    return samples.loc[wc.sample, "fastq_r2"]


def target_bed_for(wc):
    if "target_bed" in samples.columns and pd.notna(samples.loc[wc.sample, "target_bed"]):
        return samples.loc[wc.sample, "target_bed"]
    return config.get("target_bed", "")


def target_bed_input(wc):
    # Used as a rule `input`, not just a `params` string, so that Snakemake
    # actually tracks the target BED as a dependency: it fails at DAG-build
    # time if the file is missing, and re-runs the rule if the BED changes.
    bed = target_bed_for(wc)
    return [bed] if bed else []


SNP_FILTER = "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0"
INDEL_FILTER = "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0"


rule fastqc:
    input:
        r1=fastq_r1,
        r2=fastq_r2,
    output:
        directory("results/{sample}/fastqc"),
    shell:
        "mkdir -p {output} && fastqc -o {output} {input.r1} {input.r2}"


rule trim:
    input:
        r1=fastq_r1,
        r2=fastq_r2,
    output:
        r1="results/{sample}/trimmed/{sample}_R1.fastq.gz",
        r2="results/{sample}/trimmed/{sample}_R2.fastq.gz",
        json="results/{sample}/trimmed/{sample}.fastp.json",
        html="results/{sample}/trimmed/{sample}.fastp.html",
    threads: config.get("threads", 4)
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "--json {output.json} --html {output.html} --thread {threads}"


rule align:
    input:
        r1="results/{sample}/trimmed/{sample}_R1.fastq.gz",
        r2="results/{sample}/trimmed/{sample}_R2.fastq.gz",
        ref=config["reference_fasta"],
    output:
        "results/{sample}/aligned/{sample}.sorted.bam",
    threads: config.get("threads", 4)
    params:
        rg=lambda wc: f"@RG\\tID:{wc.sample}\\tSM:{wc.sample}\\tPL:ILLUMINA",
    shell:
        "bwa mem -t {threads} -R '{params.rg}' {input.ref} {input.r1} {input.r2} "
        "| samtools sort -@ {threads} -o {output} -"


rule index_bam:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    shell:
        "samtools index {input}"


rule mark_duplicates:
    input:
        bam="results/{sample}/aligned/{sample}.sorted.bam",
        bai="results/{sample}/aligned/{sample}.sorted.bam.bai",
    output:
        bam="results/{sample}/dedup/{sample}.dedup.bam",
        metrics="results/{sample}/dedup/{sample}.dedup_metrics.txt",
    shell:
        "gatk MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics}"


rule call_variants:
    input:
        bam="results/{sample}/dedup/{sample}.dedup.bam",
        bai="results/{sample}/dedup/{sample}.dedup.bam.bai",
        ref=config["reference_fasta"],
        target_bed=target_bed_input,
    output:
        "results/{sample}/variants/{sample}.raw.vcf.gz",
    params:
        target=lambda wc: f"-L {target_bed_for(wc)}" if target_bed_for(wc) else "",
    shell:
        "gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output} {params.target}"


rule select_snps:
    input:
        vcf="results/{sample}/variants/{sample}.raw.vcf.gz",
        ref=config["reference_fasta"],
    output:
        "results/{sample}/variants/{sample}.snp.vcf.gz",
    shell:
        "gatk SelectVariants -R {input.ref} -V {input.vcf} --select-type-to-include SNP -O {output}"


rule filter_snps:
    input:
        vcf="results/{sample}/variants/{sample}.snp.vcf.gz",
        ref=config["reference_fasta"],
    output:
        "results/{sample}/variants/{sample}.snp.filtered.vcf.gz",
    shell:
        'gatk VariantFiltration -R {input.ref} -V {input.vcf} '
        '--filter-expression "{SNP_FILTER}" --filter-name SNP_HARD_FILTER -O {output}'


rule select_indels:
    input:
        vcf="results/{sample}/variants/{sample}.raw.vcf.gz",
        ref=config["reference_fasta"],
    output:
        "results/{sample}/variants/{sample}.indel.vcf.gz",
    shell:
        "gatk SelectVariants -R {input.ref} -V {input.vcf} --select-type-to-include INDEL -O {output}"


rule filter_indels:
    input:
        vcf="results/{sample}/variants/{sample}.indel.vcf.gz",
        ref=config["reference_fasta"],
    output:
        "results/{sample}/variants/{sample}.indel.filtered.vcf.gz",
    shell:
        'gatk VariantFiltration -R {input.ref} -V {input.vcf} '
        '--filter-expression "{INDEL_FILTER}" --filter-name INDEL_HARD_FILTER -O {output}'


rule merge_filtered_variants:
    input:
        snp="results/{sample}/variants/{sample}.snp.filtered.vcf.gz",
        indel="results/{sample}/variants/{sample}.indel.filtered.vcf.gz",
    output:
        "results/{sample}/variants/{sample}.filtered.vcf.gz",
    shell:
        "gatk MergeVcfs -I {input.snp} -I {input.indel} -O {output}"


# Only reachable if config["annovar"]["humandb_dir"] is set -- see wgs_targets()/wes_targets()
# in wgs.smk / wes.smk, which don't request this output otherwise.
rule annotate_variants:
    input:
        vcf="results/{sample}/variants/{sample}.filtered.vcf.gz",
    output:
        "results/{sample}/annotation/{sample}.exonic_variant_function",
    params:
        avinput="results/{sample}/annotation/{sample}.avinput",
        outprefix="results/{sample}/annotation/{sample}",
        convert=config.get("annovar", {}).get("convert2annovar_script", ""),
        annotate=config.get("annovar", {}).get("annotate_variation_script", ""),
        buildver=config.get("annovar", {}).get("buildver", ""),
        humandb=config.get("annovar", {}).get("humandb_dir", ""),
    shell:
        "perl {params.convert} -format vcf4 {input.vcf} --outfile {params.avinput} --includeinfo && "
        "perl {params.annotate} -geneanno -buildver {params.buildver} {params.avinput} {params.humandb} "
        "-outfile {params.outprefix}"
