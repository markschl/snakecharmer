
localrules:
    combine_sample_reports,


rule combine_sample_reports:
    params:
        path_pattern=r"(?P<run>[^/]+?)_(?P<layout>single(\.rev)?|paired)/sample_report\.tsv",
    input:
        reports=lambda wildcards: run_results("/sample_report.tsv", workflows=[wildcards.workflow]),
    output:
        report="results/{workflow}/sample_report.tsv"
    log:
        "logs/{workflow}/combine_sample_reports.log"
    script:
        "../scripts/combine_sample_reports.py"


rule orient_otus:
    input:
        otus="results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/clusters.fasta",
    output:
        otus="results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/clusters_fwd.fasta",
    log:
        "logs/{workflow}/{run}_{layout}/{primers}/{cluster}_orient_otus.log"
    shell:
        """
        exec &> {log}
        set -xeuo pipefail
        if [ "{wildcards.layout}" = "single.rev" ]; then
            vsearch --fastx_revcomp "{input.otus}" --fastaout "{output.otus}"
        else
            ln -srf "{input.otus}" "{output.otus}"
        fi        
        """
