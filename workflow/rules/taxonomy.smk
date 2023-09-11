"""
Rules for obtaining and workdir taxonomy databases, and some taxonomy-related
sequence workdir code.

Directory structure of 'refdb':
refdb/
  taxonomy/
    db_(regular|preformatted)_<source config hash>/
      db.fasta.zst
      db_config.yaml
      flt_<filter config hash>/
        filter_config.yaml
        db.fasta.zst
        db_full_ranks.fasta.zst 
        <database format>/
          cnv_<conversion config hash>.fasta.zst  (preformatted databases land here directly, without .fasta)

Ultimately, filtered and converted (to <format>) databases are used as input
for taxonomic assignments.
If no filtering is done, the 'unfiltered' is used as ID instead of the filter config hash.
If there are no conversion options, 'standard' is used as ID instead of the conversion option hash.
"""

localrules:
    clean_tax,
    write_taxdb_source_config,
    write_taxdb_filter_config,
    write_taxdb_conversion_config,


# these formats can be imported by 'obtain_taxdb'
cfg.imported_tax_formats = ("unite_otus", "midori", "gtdb", "utax", "qiime", "qiime_qza")


# This rule writes the taxonomy database configuration to a file used as input
# for obtain_taxdb, which is a cached rule, and thus cannot use wildcards.
# The config file names contain a SHA-256 hash of the database configuration,
# so after every change a new config file will be generated, leading to a new
# download of the database.
rule write_taxdb_source_config:
    params:
        dbconfig=lambda wildcards: cfg.taxdb_sources_by_hash[wildcards.source_id],
        exclude=["source_id", "name", "preformatted"]
    output:
        yml="refdb/taxonomy/db_{preformatted}_{source_id}/db_config.yaml",
    log:
        "logs/taxonomy/db_{preformatted}_{source_id}/write_db_config.log"
    conda:
        "envs/basic.yaml"
    group:
        "taxonomy"
    script:
        "../scripts/taxonomy/write_db_config.py"


rule obtain_taxdb:
    input:
        yml="refdb/taxonomy/db_regular_{source_id}/db_config.yaml",
    output:
        db="refdb/taxonomy/db_regular_{source_id}/db.fasta.zst",
    log:
        "logs/taxonomy/db_regular_{source_id}/obtain_taxdb.log",
    wildcard_constraints:
        source_id = "[a-z0-9]+",
    cache:
        True
    conda:
        "envs/taxonomy.yaml"
    group:
        "taxonomy"
    script:
        "../scripts/taxonomy/obtain_taxdb.py"


# This rule writes a database filter config file (using SHA-256 of keys/values).
# Similarly to write_taxdb_config, this ensures caching of databases
# as well as re-filtering if any setting changs.
rule write_taxdb_filter_config:
    params:
        dbconfig=lambda wildcards: cfg.taxdb_filter_by_hash[wildcards.filter_id],
        exclude=["filter_id", "db", "name"],
    output:
        yml="refdb/taxonomy/db_{preformatted}_{source_id}/flt_{filter_id}/filter_config.yaml",
    log:
        "logs/taxonomy/db_{preformatted}_{source_id}/write_filter_config_{filter_id}.log"
    conda:
        "envs/basic.yaml"
    group:
        "taxonomy"
    script:
        "../scripts/taxonomy/write_db_config.py"


rule filter_taxdb:
    input:
        all="refdb/taxonomy/db_regular_{source_id}/db.fasta.zst",
        cfg="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/filter_config.yaml",
    output:
        filtered="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/db.fasta.zst",
    log:
        "logs/taxonomy/db_regular_{source_id}/filter_db_{filter_id}.log",
    wildcard_constraints:
        source_id = "[a-z0-9]+",
        filter_id = "[a-z0-9]+",
    conda:
        "envs/taxonomy.yaml"
    group:
        "taxonomy"
    script:
        "../scripts/taxonomy/filter_taxdb.py"


# This rule writes a database conversion/training config file
# (using SHA-256 of keys/values).
rule write_taxdb_conversion_config:
    params:
        dbconfig=lambda wildcards: cfg.taxdb_training_cfg_by_hash[wildcards.cnv_id],
        exclude=["conversion_id"],
    output:
        yml="refdb/taxonomy/db_{preformatted}_{source_id}/flt_{filter_id}/{format}/conversion_config_{cnv_id}.yaml",
    log:
        "logs/taxonomy/db_{preformatted}_{source_id}/flt_{filter_id}/write_conversion_config_{format}_{cnv_id}.log"
    conda:
        "envs/basic.yaml"
    group:
        "taxonomy"
    script:
        "../scripts/taxonomy/write_db_config.py"


# TODO: this returns a compressed database, even though trained datasets are usually already compressed
rule obtain_taxdb_preformatted:
    input:
        params="refdb/taxonomy/db_preformatted_{source_id}/db_config.yaml",
    output:
        db="refdb/taxonomy/db_preformatted_{source_id}/flt_unfiltered/{format}/cnv_standard.zst",
    log:
        "logs/taxonomy/db_preformatted_{source_id}/flt_unfiltered/{format}/obtain_{format}_standard.log",
    wildcard_constraints:
        source_id = "[a-z0-9]+",
        filter_id = "[a-z0-9]+",
    cache:
        True
    conda:
        "envs/taxonomy.yaml"
    group:
        "taxonomy"
    script:
        "../scripts/taxonomy/obtain_taxdb_formatted.py"


rule uncompress_taxdb:
    input:
        db="refdb/taxonomy/{sub_path}.zst",
    output:
        db=temp("refdb/taxonomy/{sub_path}"),
    log:
        "logs/taxonomy/{sub_path}_uncompress.log",
    wildcard_constraints:
        sub_path=r".+?(?<!\.zst)$"
    group:
        "taxonomy"
    conda:
        "envs/basic.yaml"
    shell:
        """
        zstd -dcqf "{input.db}" > "{output.db}"
        """


rule taxdb_normalize_internal_ranks:
    input:
        db="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/db.fasta.zst",
    output:
        db=temp("refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/db_full_ranks.fasta.zst"),
    log:
        "logs/taxonomy/db_regular_{source_id}/flt_{filter_id}/normalize_internal_ranks.log",
    wildcard_constraints:
        source_id = "[a-z0-9]+",
        filter_id = "[a-z0-9]+",
    group:
        "taxonomy"
    conda:
        "envs/taxonomy.yaml"
    script:
        "../scripts/taxonomy/rank_propagate.py"


rule make_tax_fasta:
    input:
        fa="results/{workflow}/workflow_{cluster}/{run}/{primers}/clusters.fasta",
        tax="results/{workflow}/workflow_{cluster}/{run}/{primers}/taxonomy/{tax_name}.txt.gz",
    output:
        "results/{workflow}/workflow_{cluster}/{run}/{primers}/taxonomy/fasta/{tax_name}.fasta.gz",
    log:
        "logs/{workflow}/{run}/{primers}/taxonomy/{cluster}_make_tax_fasta_{tax_name}.log",
    conda:
        "envs/basic.yaml"
    group:
        "taxonomy"
    shell:
        """
        tax={input.tax}
        st set -ul <(gzip -dc "$tax") -d {{l:2}} "{input.fa}" |
          st replace -d '__' ':' |
          st replace -dr ' *; *' ' ' |
          gzip -nc > {output}
        """


rule make_tax_biom:
    input:
        biom="results/{workflow}/workflow_{cluster}/{run}/{primers}/otutab.biom",
        biom_hdf5="results/{workflow}/workflow_{cluster}/{run}/{primers}/otutab.hdf5.biom",
        tax="results/{workflow}/workflow_{cluster}/{run}/{primers}/taxonomy/{tax_name}.txt.gz",
    output:
        tax_tmp=temp(
            "workdir/{workflow}/workflow_{cluster}/{run}/{primers}/_tax_tmp/{tax_name}.txt"
        ),
        biom_tmp=temp(
            "workdir/{workflow}/{run}/{primers}/{cluster}_tax_tmp/{tax_name}.biom"
        ),
        biom="results/{workflow}/workflow_{cluster}/{run}/{primers}/taxonomy/{tax_name}.biom.gz",
        biom_hdf5="results/{workflow}/workflow_{cluster}/{run}/{primers}/taxonomy/{tax_name}.hdf5.biom.gz",
    log:
        "logs/{workflow}/{run}/{primers}/taxonomy/{cluster}_make_tax_biom_{tax_name}.log",
    conda:
        "envs/biom.yaml"
    group:
        "taxonomy"
    shell:
        """
        exec &> {log}
        set -xeuo pipefail
        mkdir -p $(dirname {output.tax_tmp})
        gzip -dc {input.tax} | 
        sed 's/Taxon/taxonomy/g' |
        sed 's/Feature ID/# Feature ID/g' > {output.tax_tmp}
        if [[ $(wc -l < "{output.tax_tmp}") -ge 2 ]]; then
                biom add-metadata -i {input.biom}  \
                    -o /dev/stdout \
                    --observation-metadata-fp {output.tax_tmp} \
                    --sc-separated taxonomy --float-fields Confidence --output-as-json |
                gzip -nc > {output.biom}
                biom add-metadata -i {input.biom_hdf5}  \
                    -o {output.biom_tmp} \
                    --observation-metadata-fp {output.tax_tmp} \
                    --sc-separated taxonomy --float-fields Confidence
                gzip -nc {output.biom_tmp} > {output.biom_hdf5}
        else
            # no taxa
            echo -n | gzip -nc > {output.biom}
            echo -n | gzip -nc > {output.biom_hdf5}
        fi
        """


rule clean_tax:
    shell:
        "rm -Rf results/*/workflow_*/*/*/taxonomy"
