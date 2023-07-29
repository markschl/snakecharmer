
localrules:
    clean_fastqc,


rule fastqc:
    input:
        "input/grouped/{strategy}/{sample}/{prefix}.fastq.gz",
    output:
        html="results/_validation/fastqc/{strategy,[^/]+}/{sample,[^/]+}/{prefix}_fastqc.html",
        zip="results/_validation/fastqc/{strategy,[^/]+}/{sample,[^/]+}/{prefix}_fastqc.zip",
    log:
        "logs/fastqc/{strategy}/{sample}/{prefix}.log",
    group:
        "qc"
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        outdir=$(dirname {output.html})
        mkdir -p $outdir
        fastqc -q -f fastq -t 1 -o "$outdir" {input} 2> {log}
        """


rule multiqc_fastqc:
    input:
        [
            join(
                "results",
                "_validation",
                "fastqc",
                splitext(splitext(join(path, name))[0])[0] + "_fastqc.html",
            )
            for _, path, name in link_paths_flat
        ],
    output:
        "results/_validation/multiqc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        outdir="$(dirname {output})"
        multiqc -fm fastqc -o "$outdir" results/_validation/fastqc 2> {log}
        (cd "$outdir" && zip -FSqr -rm multiqc_data.zip multiqc_data) 2> {log}
        """


rule clean_fastqc:
    shell:
        "rm -Rf results/_validation/fastqc"
