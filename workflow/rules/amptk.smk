import os
from os.path import basename
import lib


localrules:
    amptk_stats_paired,


rule amptk_merge_trim:
    params:
        err_rate=lambda w: cfg[w.workflow]["settings"]["primers"]["trim_settings"]["max_error_rate"],
        min_len=lambda w: cfg[w.workflow]["settings"]["usearch"]["filter"]["min_length"] \
            if cfg[w.workflow]["cluster"] == "amptk_unoise3" \
            else cfg[w.workflow]["settings"]["dada2"]["min_length"],
        program=lambda w: with_default(cfg[w.workflow]["settings"]["usearch"], "merge", "program"),
        usearch_bin=config["software"]["usearch"]["binary"],
    input:
        primers_yaml="workdir/primers/primers.yaml",
        sample_tab=lambda wildcards: "workdir/{{workflow}}/input/sample_config/{technology}/paired/{run}/samples.tsv".format(
            **cfg.get_run_data(layout="paired", **wildcards)
        ),
        fq=lambda wildcards: expand_input_files(layout="paired", **wildcards),
    output:
        demux="workdir/{workflow}/{run}_paired/{marker}__{f_primer}...{r_primer}/illumina.demux.fq.gz",
        mapping="workdir/{workflow}/{run}_paired/{marker}__{f_primer}...{r_primer}/illumina.mapping_file.txt",
        log="workdir/{workflow}/{run}_paired/{marker}__{f_primer}...{r_primer}/illumina.amptk-demux.log",
    log:
        "logs/{workflow}/{run}_paired/{marker}__{f_primer}...{r_primer}/amptk_trim_merge.log",
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

# ruleorder:
#     amptk_denoise > otutab_to_biom

rule amptk_denoise:
    params:
        method=lambda wildcards: wildcards.cluster,
        usearch_par=lambda wildcards: cfg[wildcards.workflow]["settings"]["usearch"],
        dada2_par=lambda wildcards: cfg[wildcards.workflow]["settings"]["dada2"],
        unoise_program=lambda w: with_default(cfg[w.workflow]["settings"]["usearch"], "unoise3", "program"),
        usearch_bin=config["software"]["usearch"]["binary"],
    input:
        demux="workdir/{workflow}/{run}_paired/{primers}/illumina.demux.fq.gz",
    output:
        denoised="results/{workflow}/workflow_amptk_{cluster}/{run}_paired/{primers}/denoised.fasta",
        tab="results/{workflow}/workflow_amptk_{cluster}/{run}_paired/{primers}/denoised_otutab.txt.gz",
    log:
        "logs/{workflow}/{run}_paired/{primers}/{cluster}_amptk_denoise.log",
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


rule amptk_combine_logs:
    input:
        merge_trim=expand(
            "logs/{{workflow}}/{{run}}_paired/{primers}/amptk_trim_merge.log",
            primers=cfg.primer_combinations_flat
        ),
        denoise=expand(
            "logs/{{workflow}}/{{run}}_paired/{primers}/{{cluster}}_amptk_denoise.log",
            primers=cfg.primer_combinations_flat
        ),
    output:
        "logs/{workflow}/{run}_paired_amptk_{cluster}_all.log",
    shell:
        """
        exec 1> "{output}"
        echo "Read merging and primer trimming"
        echo "================================"
        cat {input.merge_trim:q}
        printf "\n\n\n"
        echo "Denoising"
        echo "=========="
        cat {input.denoise:q}
        """


##########################
#### QC
##########################



# TODO: not implemented
rule amptk_stats_paired:
    output:
        touch("results/{workflow}/workflow_amptk_{cluster}/{run}_paired/sample_report.tsv"),
    log:
        "logs/{workflow}/{run}_paired/amptk_{cluster}_sample_report.log",
