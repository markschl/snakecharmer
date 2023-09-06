
rule fastqc:
    input:
        sample_dir=rules.collect_sample_files.output.sample_dir,
    output:
        qc_dir=directory("input/fastqc/{technology}/{layout}/{demux}/{run}")
    log:
        "logs/qc/{technology}/{layout}/{demux}/{run}.log",
    threads: workflow.cores,
    group:
        "run"
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        mkdir -p {output.qc_dir}
        fastqc -q -f fastq -t {threads} -o {output.qc_dir} {input.sample_dir}/*.fastq.gz 2> {log}
        """


rule multiqc_fastqc:
    input:
        fastqc=lambda wildcards: rules.fastqc.output.qc_dir.format(
            demux="demux",
            **cfg.get_run_data(pooled=False, **wildcards)
        ),
    output:
        "results/_qc/multiqc_{run}_{layout}/multiqc_report.html",
    log:
        "logs/qc/multiqc_{run}_{layout}.log",
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        outdir="$(dirname "{output}")"
        multiqc -fm fastqc -o "$outdir" "{input.fastqc}" 2> {log}
        (cd "$outdir" && zip -FSqr -rm multiqc_data.zip multiqc_data) 2> {log}
        """


ruleorder:
    multiqc_fastqc_pooled_workdir > multiqc_link_workdir


rule multiqc_link_workdir:
    input:
        fastqc=rules.multiqc_fastqc.input.fastqc,
    output:
        "results/{workflow}/workflow_{pipeline}/{run}_{layout}/_qc/multiqc_report.html",
    log:
        "logs/{workflow}/qiime/workflow_{pipeline}/{run}_{layout}/qc/multiqc_report.log",
    priority:
        -1
    shell:
        """
        ln -srf {input.existing} {output} 2> {log}
        """


rule multiqc_fastqc_pooled_workdir:
    input:
        fastqc=lambda wildcards: rules.fastqc.output.qc_dir.format(
            demux="demux",
            **cfg.get_run_data(**wildcards)
        ),
    output:
        "results/{workflow}/workflow_{pipeline}/{run}_{layout}/_qc/multiqc_report.html",
    log:
        "logs/{workflow}/qiime/workflow_{pipeline}/{run}_{layout}/qc/multiqc_report.log",
    priority:
        -1
    shell:
        """
        outdir="$(dirname "{output}")"
        multiqc -fm fastqc -o "$outdir" "{input.fastqc}" 2> {log}
        (cd "$outdir" && zip -FSqr -rm multiqc_data.zip multiqc_data) 2> {log}
        """
