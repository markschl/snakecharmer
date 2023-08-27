"""
Rules for obtaining and processing taxonomy databases, and some taxonomy-related
sequence processing code.

Directory structure of 'refdb':
refdb:
  taxonomy:
    db_(regular|preformatted)_<source config hash>
      <format>.fasta.zst
      <format>_db_config.yaml
      flt_<filter config hash>:
        filter_config.yaml
        <format.fasta.zst>
        <other format>.fasta.zst

Ultimately, filtered and converted (to <format>) databases are used as input
for taxonomic assignments. If not filter configuration was supplied, 
"""

localrules:
    clean_tax,
    write_taxdb_config,
    write_taxdb_filter_config,


# This rule writes the taxonomy database configuration to a file used as input
# for obtain_taxdb, which is a cached rule, and thus cannot use wildcards.
# The config file names contain a SHA-256 hash of the database configuration,
# so after every change a new config file will be generated, leading to a new
# download of the database.
rule write_taxdb_config:
    params:
        dbconfig=lambda wildcards: cfg.taxdb_sources_by_hash[wildcards.source_id],
        exclude=["source_id", "name", "preformatted"]
    output:
        yml="refdb/taxonomy/db_{preformatted}_{source_id}/{type}_db_config.yaml",
    log:
        "logs/taxonomy/db_{preformatted}_{source_id}/write_db_config_{type}.log"
    wildcard_constraints:
        preformatted = "\w+",
        source_id = "\w+",
        type = "\w+",
    conda:
        "envs/basic.yaml"
    group:
        "taxonomy"
    script:
        "../scripts/taxonomy/write_db_config.py"


rule obtain_taxdb:
    input:
        yml="refdb/taxonomy/db_{preformatted}_{source_id}/{type}_db_config.yaml",
    output:
        db="refdb/taxonomy/db_{preformatted}_{source_id}/{type}.fasta.zst",
    log:
        "logs/taxonomy/db_{preformatted}_{source_id}/obtain_{type}.log",
    wildcard_constraints:
        preformatted = "\w+",
        source_id = "\w+",
        type = "\w+",
    cache:
        True
    conda:
        "envs/taxonomy.yaml"
    group:
        "taxonomy"
    script:
        "../scripts/taxonomy/obtain_taxdb.py"


# This rule writes a database filter config file (using SHA-256 of keys/values).
# Similarly to write_taxdb_config, this ensures caching of trained databases.
rule write_taxdb_filter_config:
    params:
        dbconfig=lambda wildcards: cfg.taxdb_filter_by_hash[wildcards.filter_id],
        exclude=["filter_id", "source", "db", "name"],
    output:
        yml="refdb/taxonomy/db_{preformatted}_{source_id}/flt_{filter_id}/filter_config.yaml",
    log:
        "logs/taxonomy/db_{preformatted}_{source_id}/write_filter_config_{filter_id}.log"
    wildcard_constraints:
        preformatted = "\w+",
        source_id = "\w+",
        filter_id = "\w+",
    conda:
        "envs/basic.yaml"
    group:
        "taxonomy"
    script:
        "../scripts/taxonomy/write_db_config.py"


rule filter_taxdb:
    input:
        all="refdb/taxonomy/db_{source_id}/qiime.fasta.zst",
        cfg="refdb/taxonomy/db_{source_id}/flt_{filter_id}/filter_config.yaml",
    output:
        filtered="refdb/taxonomy/db_{source_id}/flt_{filter_id}/qiime.fasta.zst",
    log:
        "logs/taxonomy/db_{source_id}/filter_db_{filter_id}.log",
    wildcard_constraints:
        preformatted = "\w+",
        source_id = "\w+",
        filter_id = "\w+",
        dbfile = "[^/]+",
    conda:
        "envs/taxonomy.yaml"
    group:
        "taxonomy"
    script:
        "../scripts/taxonomy/filter_taxdb.py"


rule make_tax_fasta:
    input:
        fa="results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/denoised.fasta",
        tax="results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/taxonomy/{tax_name}.txt.gz",
    output:
        "results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/taxonomy/fasta/{tax_name}.fasta.gz",
    log:
        "logs/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/taxonomy/fasta/make_tax_fasta_{tax_name}.log",
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
        biom="results/{name}/{workflow}/{primers}/{layout}/denoised.biom",
        biom_hdf5="results/{name}/{workflow}/{primers}/{layout}/denoised.hdf5.biom",
        tax="results/{name}/{workflow}/{primers}/{layout}/taxonomy/{tax_name}.txt.gz",
    output:
        tax_tmp=temp(
            "processing/{name}/{workflow}/{primers}/{layout}/_tax_tmp/{tax_name}.txt"
        ),
        biom_tmp=temp(
            "processing/{name}/{workflow}/{primers}/{layout}/_tax_tmp/{tax_name}.biom"
        ),
        biom="results/{name}/{workflow}/{primers}/{layout}/taxonomy/{tax_name}.biom.gz",
        biom_hdf5="results/{name}/{workflow}/{primers}/{layout}/taxonomy/{tax_name}.hdf5.biom.gz",
    log:
        "logs/{name}/other/{layout}/{workflow}/{primers}/make_tax_biom/{tax_name}.log",
    conda:
        "envs/biom.yaml"
    group:
        "taxonomy"
    shell:
        """
        {{
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
        }} 2> {log}
        """


rule clean_tax:
    shell:
        "rm -Rf results/*/workflow_*/*/*/taxonomy"
