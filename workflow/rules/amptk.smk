import os
from os.path import basename
import lib


localrules:
    amptk_collect,
    amptk_stats_paired,


# TODO: usearch not checked in conda environment


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
    run:
        with lib.file_logging(log) as out:
            for source, target in zip(input, output):
                assert basename(source) == basename(target)
                if os.path.exists(target):
                    os.remove(target)
                outdir = os.path.dirname(target)
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                os.symlink(os.path.abspath(source), os.path.abspath(target))
                source = os.path.relpath(source, ".")
                target = os.path.relpath(target, ".")
                print("{} > {}".format(source, target), file=out)


rule amptk_merge_trim:
    params:
        f_primer_seq=lambda w: cfg.primers_consensus[w.marker]["forward"][w.f_primer],
        r_primer_seq=lambda w: cfg.primers_consensus[w.marker]["reverse"][w.r_primer],
        # note: we limit the max. number of primer mismatches to 2.
        # The reason is that Amptk apparently does **another** primer trimming after merging 
        # the already trimmed reads.
        # Too liberal mismatch thresholds lead to many unspecific primer matches
        # and consequently to unwanted trimming of reads.
        pmismatch=lambda w: min(2, round((
                len(cfg.primers_consensus[w.marker]['forward'][w.f_primer])
                + len(cfg.primers_consensus[w.marker]['reverse'][w.r_primer])
            ) / 2
            * cfg[w.name]["settings"]["primers"]["trim_settings"]["max_error_rate"]
        )),
        min_len=lambda w: cfg[w.name]["settings"]["filter"]["min_length"],
        program=lambda w: cfg[w.name]["settings"]["usearch"]["merge"]["program"],
    input:
        expand(
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
    shell:
        """
        prefix={output.demux}
        prefix="${{prefix%.demux.fq.gz}}"
        outdir=$(dirname $prefix)
        mkdir -p "$outdir"
        (   indir=$(realpath --relative-to=$outdir processing/{wildcards.name}/amptk/input/grouped/paired)
            cd $outdir
            amptk illumina -i $indir -o $(basename $prefix) \
                -f {params.f_primer_seq} -r {params.r_primer_seq} \
                --min_len {params.min_len} \
                --trim_len 10000000 `# high enough to never be longer`  \
                --cpus {threads} \
                --cleanup \
                --require_primer=on \
                --rescue_forward=off \
                --primer_mismatch {params.pmismatch} \
                --merge_method {params.program} \
                --usearch $(which usearch)  ) 2> {log} >/dev/null
        """


def amptk_denoise_params(method, settings):
    upar = settings["usearch"]
    maxee = upar["merge"]["expected_length"] * upar["filter"]["max_error_rate"]
    if method == "unoise3":
        return " ".join(
            [
                "--usearch",
                "$(which usearch)",
                "--maxee",
                str(maxee),
                "--method",
                upar["unoise"]["program"],
                "--minsize",
                str(upar["unoise"]["min_size"]),
            ]
        )
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
        args=lambda w: amptk_denoise_params(w.method, cfg[w.name]["settings"]),
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
