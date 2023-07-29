
localrules:
    clean_taxdb,
    clean_tax,


rule obtain_taxdb:
    params:
        par=lambda w: config["taxonomy_dbs"][w.marker][w.db],
    output:
        "refdb/taxonomy/{marker}/{db}/all/seqs.fasta.zst",
    log:
        "logs/taxdb/obtain/{marker}/{db}.log",
    conda:
        "envs/basic.yaml"
    # cache:
    #     True
    group:
        "taxonomy"
    shell:
        """
        $PIPELINE_DIR/workflow/scripts/taxonomy/obtain_taxdb.sh {params.par[type]} {params.par[url]} {output} 2> {log}
        """


rule filter_taxdb:
    input:
        "refdb/taxonomy/{marker}/{db}/all/seqs.fasta.zst",
    output:
        "refdb/taxonomy/{marker}/{db}/filtered/{defined}.fasta.zst",
    log:
        "logs/taxdb/filter/{marker}/{db}/{defined}.log",
    conda:
        "envs/basic.yaml"
    group:
        "taxonomy"
    shell:
        """
        if [[ "{wildcards.defined}" == "_all" ]]; then
          # nothing to filter
          mkdir -p $(dirname {output})
          ln -srf {input} {output}
        else
          $PIPELINE_DIR/workflow/scripts/taxonomy/filter_taxdb.sh {input} {wildcards.defined} {output} 2> {log}
        fi
        """


rule make_tax_fasta:
    input:
        fa="results/{name}/{pipeline}/{primers}/{strategy}/denoised.fasta",
        tax="results/{name}/{pipeline}/{primers}/{strategy}/taxonomy/{tax_name}.txt.gz",
    output:
        "results/{name}/{pipeline}/{primers}/{strategy}/taxonomy/fasta/{tax_name}.fasta.gz",
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
        biom="results/{name}/{pipeline}/{primers}/{strategy}/denoised.biom",
        tax="results/{name}/{pipeline}/{primers}/{strategy}/taxonomy/{tax_name}.txt.gz",
    output:
        tax_tmp=temp(
            "processing/{name}/{pipeline}/{primers}/{strategy}/_tax_tmp/{tax_name}.txt"
        ),
        biom="results/{name}/{pipeline}/{primers}/{strategy}/taxonomy/{tax_name}.biom.gz",
    log:
        "logs/{name}/other/{strategy}/{pipeline}/{primers}/make_tax_biom/{tax_name}.log",
    conda:
        "envs/biom.yaml"
    group:
        "taxonomy"
    shell:
        """
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
        else
            echo -n | gzip -nc > {output.biom}
        fi
        """


rule clean_taxdb:
    shell:
        "rm -Rf refdb"


rule clean_tax:
    shell:
        "rm -Rf results/*/pipeline_*/*/*/taxonomy"
