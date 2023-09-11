

cfg.taxonomy_formats["idtaxa"] = "idtaxa"


rule idtaxa_train:
    input:
        params="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/idtaxa/conversion_config_{cnv_id}.yaml",
        seq="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/db_full_ranks.fasta",
    output:
        db="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/idtaxa/cnv_{cnv_id}",
    wildcard_constraints:
        source_id = "\w+",
        filter_id = "\w+",
        cnv_id = "\w+",
    cache: True
    conda:
        "envs/decipher.yaml"
    log:
        "logs/taxonomy/db_regular_{source_id}/flt_{filter_id}/idtaxa_train-{cnv_id}.log",
    resources:
        mem_mb=40000,
        runtime=24 * 60,
    script:
        "../scripts/taxonomy/idtaxa_train.R"


rule idtaxa_classify:
    params:
        par=lambda wildcards: cfg.tax_config(**wildcards)["assign"],
    input:
        seq="results/{workflow}/workflow_{cluster}/{run}/{marker}__{primers}/clusters.fasta",
        db=lambda wildcards: get_refdb_path(format="idtaxa", **wildcards),
    output:
        tax="results/{workflow}/workflow_{cluster}/{run}/{marker}__{primers}/taxonomy/{db_name}-idtaxa-{tax_method}.txt.gz",
    log:
        "logs/{workflow}/{run}/{marker}__{primers}/taxonomy/{cluster}_{db_name}-{tax_method}.log",
    conda:
        "envs/decipher.yaml"
    threads: workflow.cores,
    resources:
        mem_mb=50000,
        runtime=36 * 60,
    script:
        "../scripts/taxonomy/idtaxa_assign.R"
