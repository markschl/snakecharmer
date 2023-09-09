

rule combine_sample_reports:
    params:
        path_pattern="(?P<run>[^/]+?)_(?P<layout>single|paired)/sample_report\.tsv",
    input:
        reports=lambda wildcards: run_results("/sample_report.tsv", workflows=[wildcards.workflow]),
    output:
        report="results/{workflow}/sample_report.tsv"
    log:
        "logs/{workflow}/combine_sample_reports.log"
    script:
        "../scripts/combine_sample_reports.py"
