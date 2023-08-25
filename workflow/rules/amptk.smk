import os
from os.path import basename
import lib


localrules:
    amptk_stats_paired,
    amptk_multiqc,


rule amptk_merge_trim:
    params:
        err_rate=lambda w: cfg[w.workflow]["settings"]["primers"]["trim_settings"]["max_error_rate"],
        min_len=lambda w: cfg[w.workflow]["settings"]["filter"]["min_length"],
        program=lambda w: cfg[w.workflow]["settings"]["usearch"]["merge"]["program"],
        usearch_bin=config["software"]["usearch"]["binary"],
    input:
        primers_yaml="processing/primers/primers.yaml",
        fq=lambda wildcards: expand_samples(
            "processing/{workflow}/input/{technology}/{layout}/{run}/{sample}_R{read}.fastq.gz",
            technology="illumina",
            layout="paired",
            pool=cfg[wildcards.workflow]["settings"]["pool_raw"],
            **wildcards
        ),
    output:
        demux="processing/{workflow}/amptk/{run}/paired/{marker}__{f_primer}...{r_primer}/illumina.demux.fq.gz",
        mapping="processing/{workflow}/amptk/{run}/paired/{marker}__{f_primer}...{r_primer}/illumina.mapping_file.txt",
        log="processing/{workflow}/amptk/{run}/paired/{marker}__{f_primer}...{r_primer}/illumina.amptk-demux.log",
    log:
        "logs/{workflow}/amptk/{run}/paired/{marker}__{f_primer}...{r_primer}/trim_merge.log",
    group:
        "prepare"
    conda:
        config["software"]["amptk"]["conda_env"]
    threads: workflow.cores
    resources:
        mem_mb=10000,
        runtime=24 * 60,
    script:
        "../scripts/amptk_trim_paired.py"


rule amptk_denoise:
    params:
        method=lambda wildcards: wildcards.method,
        usearch_par=lambda wildcards: cfg[wildcards.workflow]["settings"]["usearch"],
        dada2_par=lambda wildcards: cfg[wildcards.workflow]["settings"]["dada2"],
        usearch_bin=config["software"]["usearch"]["binary"],
    input:
        demux="processing/{workflow}/amptk/{run}/paired/{primers}/illumina.demux.fq.gz",
    output:
        denoised="results/{workflow}/workflow_amptk_{method}/{run}_paired/{primers}/denoised.fasta",
        tab="results/{workflow}/workflow_amptk_{method}/{run}_paired/{primers}/denoised_otutab.txt.gz",
    log:
        "logs/{workflow}/amptk/{run}/paired/{primers}/{method}.log",
    conda:
        config["software"]["amptk"]["conda_env"]
    group:
        "denoise"
    threads: max(10, workflow.cores)  # dereplication/clustering only use one core, only mapping uses all -> don't claim too much (will be slower, though)
    resources:
        mem_mb=30000,
        runtime=36 * 60,
    script:
        "../scripts/amptk_denoise_paired.py"


##########################
#### QC
##########################


rule amptk_multiqc:
    input:
        rules.multiqc_fastqc.output,
    output:
        "results/{workflow}/validation/multiqc_amptk/multiqc_report.html",
    log:
        "logs/{workflow}/multiqc_amptk.log",
    shell:
        "ln -srf {input} {output} 2> {log}"


# TODO: not implemented
rule amptk_stats_paired:
    output:
        touch("results/{workflow}/workflow_amptk_{cluster}/{run}_paired/sample_report.tsv"),
    log:
        "logs/{workflow}/workflow_amptk_{cluster}/{run}/paired/sample_report.log",
