import os
from os.path import basename
import lib


localrules:
    amptk_collect,
    amptk_stats_paired,


rule amptk_collect:
    input:
        expand(
            "input/grouped/paired/{sample}/{sample}_R{read}.fastq.gz",
            sample=cfg.sample_names["paired"],
            read=[1, 2],
        ),
    output:
        expand(
            "processing/{{name}}/amptk/input/grouped/paired/{sample}_R{read}.fastq.gz",
            sample=cfg.sample_names["paired"],
            read=[1, 2],
        ),
    log:
        "logs/{name}/amptk/collect.log",
    script:
        "../scripts/amptk_collect.py"


rule amptk_merge_trim:
    params:
        err_rate=lambda w: cfg[w.name]["settings"]["primers"]["trim_settings"]["max_error_rate"],
        min_len=lambda w: cfg[w.name]["settings"]["filter"]["min_length"],
        program=lambda w: cfg[w.name]["settings"]["usearch"]["merge"]["program"],
        usearch_bin=config['software']['usearch']['binary']
    input:
        primers_yaml='processing/primers/primers.yaml',
        fq=expand(
            "processing/{{name}}/amptk/input/grouped/paired/{sample}_R{read}.fastq.gz",
            sample=cfg.sample_names["paired"],
            read=[1, 2],
        ),
    output:
        demux="processing/{name}/amptk/analysis/paired/{marker}__{f_primer}...{r_primer}/illumina.demux.fq.gz",
        mapping="processing/{name}/amptk/analysis/paired/{marker}__{f_primer}...{r_primer}/illumina.mapping_file.txt",
        log="processing/{name}/amptk/analysis/paired/{marker}__{f_primer}...{r_primer}/illumina.amptk-demux.log",
    log:
        "logs/{name}/amptk/paired/{marker}__{f_primer}...{r_primer}/trim_merge.log",
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


def amptk_denoise_params(method, settings, usearch_bin=None):
    upar = settings["usearch"]
    maxee = upar["merge"]["expected_length"] * upar["filter"]["max_error_rate"]
    if method == "unoise3":
        args = [
                "--usearch",
                usearch_bin,
                "--maxee",
                str(maxee),
                "--method",
                upar["unoise"]["program"],
                "--minsize",
                str(upar["unoise"]["min_size"]),
            ]
        if upar["unoise"]["program"] == "usearch":
            assert usearch_bin is not None
            args += ["--usearch", usearch_bin]
        return " ".join(args)
    if method == "dada2":
        par = settings["dada2"]
        out = [
            "--maxee",
            str(maxee),
            "--chimera_method",
            par["chimera_method"],
        ]
        p = par.get("pooling_method", "independent")
        if p == "pooled":
            out.append("--pool")
        elif p == "pseudo":
            out.append("--pseudopool")
        return " ".join(out)
    raise Exception("Unknown / unimplemented Amptk command " + method)


rule amptk_denoise:
    params:
        args=lambda w: amptk_denoise_params(w.method, cfg[w.name]["settings"], usearch_bin=config['software']['usearch']['binary']),
    input:
        demux="processing/{name}/amptk/analysis/paired/{primers}/illumina.demux.fq.gz",
    output:
        denoised="results/{name}/pipeline_amptk_{method}/{primers}/paired/denoised.fasta",
        tab="results/{name}/pipeline_amptk_{method}/{primers}/paired/denoised_otutab.txt.gz",
    log:
        "logs/{name}/amptk/paired/{primers}/{method}.log",
    conda:
        config["software"]["amptk"]["conda_env"]
    group:
        "denoise"
    threads: max(10, workflow.cores)  # dereplication/clustering only use one core, only mapping uses all -> don't claim too much (will be slower, though)
    resources:
        mem_mb=30000,
        runtime=36 * 60,
    shell:
        """
        outdir=$(dirname {input.demux})
        (  cd $outdir &&
            amptk {wildcards.method} -i $(basename {input.demux}) \
                -o {wildcards.method} \
                {params.args} \
                --cpus {threads}  ) 2> {log} >/dev/null
        # copy files
        mkdir -p $(dirname {output.denoised})
        cp $outdir/{wildcards.method}.ASVs.fa {output.denoised}
        gzip -nc $outdir/{wildcards.method}.otu_table.txt > {output.tab}
        """


##########################
#### QC
##########################


rule amptk_multiqc_paired:
    input:
        rules.multiqc_fastqc.output,
    output:
        "results/{name}/pipeline_amptk_{method}/_validation/multiqc/multiqc_report.html",
    shell:
        "ln -srf {input} {output}"


rule amptk_stats_paired:
    input:
        merge=expand(
            "processing/{{name}}/amptk/analysis/paired/{primers}/illumina.amptk-demux.log",
            primers=cfg.primer_combinations_flat,
        ),
    output:
        "results/{name}/pipeline_amptk_{method}/_validation/sample_report.tsv",
    log:
        "logs/{name}/amptk/sample_report_{method}.log",
    run:
        # TODO: not implemented
        with open(output[0], "w") as out:
            pass
