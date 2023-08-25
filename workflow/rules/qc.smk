rule fastqc:
    input:
        lambda w: expand("input/{{technology}}/{{layout}}/{{demux_method}}/{{run}}/{{sample}}_R{read}.fastq.gz",
               read=[1, 2] if w.layout == "paired" else [1]),
    output:
        expand("input/{{technology}}/{{layout}}/{{demux_method}}/{{run}}/fastqc/{{sample}}_R{read}_fastqc.{ext}",
               read=[1, 2],
               ext=("html", "zip"),
               allow_missing=True),
    log:
        "logs/qc/{technology}/{layout}/{demux_method}/{run}/{sample}.log",
    group:
        "sample"
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        outdir=$(dirname "{output[0]}")
        mkdir -p "$outdir"
        fastqc -q -f fastq -t 1 -o "$outdir" {input} 2> {log}
        """


rule multiqc_fastqc:
    input:
        lambda wildcards: expand_samples(
            path="input/{technology}/{layout}/{demux_method}/{run}/fastqc/{sample}_R{read}_fastqc.html",
            demux_method="demux",
            **wildcards
        ),
            # demux_method=cfg[wildcards.workflow]["demux_method"],
    output:
        "results/{workflow}/validation/multiqc/multiqc_report.html",
    log:
        "logs/{workflow}/multiqc.log",
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        indir=input
        outdir="$(dirname {output})"
        multiqc -fm fastqc -o "$outdir" "$indir" 2> {log}
        (cd "$outdir" && zip -FSqr -rm multiqc_data.zip multiqc_data) 2> {log}
        """
