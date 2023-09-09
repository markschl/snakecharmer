
localrules:
    qiime_make_manifest,
    qiime_combine_sample_reports,


rule qiime_make_manifest:
    params:
        qiime_style=True,
        subdir="nested"
    input:
        tab=lambda wildcards: "workdir/{{workflow}}/input/sample_config/{technology}/{layout}/{run}/samples.tsv".format(**cfg.get_run_data(**wildcards)),
        sample_dir=lambda wildcards: "workdir/{{workflow}}/input/{technology}/{layout}/{run}".format(**cfg.get_run_data(**wildcards)),
    output:
        tab="workdir/{workflow}/{run}_{layout}/qiime_manifest.txt",
    log:
        "logs/{workflow}/{run}_{layout}/qiime_make_manifest.log",
    script:
        "../scripts/make_new_sample_tab.py"


rule qiime_import:
    params:
        type=lambda wildcards: "SequencesWithQuality"
          if wildcards.layout == "single"
          else "PairedEndSequencesWithQuality",
        format=lambda wildcards: "SingleEndFastqManifestPhred33V2"
          if wildcards.layout == "single"
          else "PairedEndFastqManifestPhred33V2",
    input:
        manifest=rules.qiime_make_manifest.output.tab,
    output:
        "workdir/{workflow}/{run}_{layout}/demux.qza",
    log:
        "logs/{workflow}/{run}_{layout}/qiime_import.log",
    group:
        "prepare"
    resources:
        runtime=12 * 60,
    conda:
        config["software"]["qiime"]["conda_env"]
    shell:
        """
        qiime tools import \
            --type 'SampleData[{params.type}]' \
            --input-path {input.manifest} \
            --output-path {output} \
            --input-format {params.format} &> {log}
        """


rule qiime_trim_paired:
    params:
        err_rate=lambda w: cfg[w.workflow]["settings"]["primers"]["trim_settings"]["max_error_rate"],
        min_length=lambda w: cfg[w.workflow]["settings"]["dada2"]["min_length"],
    input:
        yaml="workdir/primers/primers.yaml",
        demux="workdir/{workflow}/{run}_paired/demux.qza",
    output:
        qza="workdir/{workflow}/{run}_paired/{marker}__{f_primer}...{r_primer}/trim.qza",
    log:
        "logs/{workflow}/{run}_paired/{marker}__{f_primer}...{r_primer}/qiime_trim.log",
    group:
        "prepare"
    conda:
        config["software"]["qiime"]["conda_env"]
    threads: workflow.cores
    resources:
        runtime=12 * 60,
    script:
        "../scripts/qiime_trim_paired.py"


rule qiime_dada2_paired:
    params:
        par=lambda w: cfg[w.workflow]["settings"]["dada2"],
    input:
        trim="workdir/{workflow}/{run}_paired/{primers}/trim.qza",
    output:
        denoised0="workdir/{workflow}/{run}_paired/{primers}/dada2.qza",
        tab0="workdir/{workflow}/{run}_paired/{primers}/dada2_tab.qza",
        stats="workdir/{workflow}/{run}_paired/{primers}/dada2_stats.qza",
    log:
        "logs/{workflow}/{run}_paired/{primers}/dada2_qiime.log",
    conda:
        config["software"]["qiime"]["conda_env"]
    group:
        "denoise"
    threads: workflow.cores
    resources:
        mem_mb=30000,
        runtime=36 * 60,
    script:
        "../scripts/qiime_denoise_paired.py"


ruleorder:
    qiime_export > otutab_to_biom > biom_to_hdf5


rule qiime_export:
    input:
        denoised="workdir/{workflow}/{run}/{primers}/{cluster}.qza",
        tab="workdir/{workflow}/{run}/{primers}/{cluster}_tab.qza",
        stats="workdir/{workflow}/{run}/{primers}/dada2_stats.qza",
    output:
        # directory("results/{workflow}/workflow_qiime_{cluster}/{run}/{primers}"),
        denoised="results/{workflow}/workflow_qiime_{cluster}/{run}/{primers}/denoised.fasta",
        tab="results/{workflow}/workflow_qiime_{cluster}/{run}/{primers}/denoised_otutab.txt.gz",
        biom_json="results/{workflow}/workflow_qiime_{cluster}/{run}/{primers}/denoised.biom",
        biom_hdf5="results/{workflow}/workflow_qiime_{cluster}/{run}/{primers}/denoised.hdf5.biom",
        stats=temp("results/{workflow}/workflow_qiime_{cluster}/{run}/{primers}/sample_report.tsv"),
        tmp=temp(
            directory("workdir/{workflow}/{run}/{primers}/{cluster}_export_tmp")
        ),
    log:
        "logs/{workflow}/{run}/{primers}/{cluster}_qiime_export.log",
    conda:
        config["software"]["qiime"]["conda_env"]
    group:
        "denoise"
    shell:
        """
        exec &> {log}
        # export table
        mkdir -p {output.tmp}
        tab={output.tab}
        tab=${{tab%.gz}}
        qiime tools export \
            --input-path {input.tab} \
            --output-path {output.tmp}
        mv {output.tmp}/feature-table.biom {output.biom_hdf5}
        biom convert -i {output.biom_hdf5}  \
            -o {output.biom_json} \
            --to-json --table-type "OTU table"
        biom convert -i {output.biom_hdf5}  \
            -o $tab \
            --to-tsv --table-type "OTU table"
        sed -i '1,1d' $tab
        gzip -nf $tab

        # export seqs
        qiime tools export \
            --input-path {input.denoised} \
            --output-path {output.tmp}
        cat {output.tmp}/dna-sequences.fasta > {output.denoised}

        # export stats
        qiime metadata tabulate \
            --m-input-file {input.stats} \
            --o-visualization {output.tmp}/stats.qzv
        qiime tools export \
            --input-path {output.tmp}/stats.qzv \
            --output-path {output.tmp}/stats
        statfile={output.tmp}/stats/metadata.tsv
        head -n1 $statfile > {output.stats}
        tail -n+3 $statfile >> {output.stats}
        """


rule qiime_combine_sample_reports:
    params:
        path_pattern="/(?P<primers>[^/]+?)/sample_report\.tsv",
    input:
        reports=expand(
            "results/{{workflow}}/workflow_qiime_{{cluster}}/{{run}}/{primers}/sample_report.tsv", 
            primers=cfg.primer_combinations_flat,
        ),
    output:
        report="results/{workflow}/workflow_qiime_{cluster}/{run}/sample_report.tsv"
    log:
        "logs/{workflow}/{run}/{cluster}_qiime_combine_sample_reports.log"
    script:
        "../scripts/combine_sample_reports.py"


rule qiime_combine_logs:
    input:
        trim=expand(
            "logs/{{workflow}}/{{run}}/{primers}/qiime_trim.log",
            primers=cfg.primer_combinations_flat
        ),
        denoise=expand(
            "logs/{{workflow}}/{{run}}/{primers}/{{cluster}}_qiime.log",
            primers=cfg.primer_combinations_flat
        ),
    output:
        "logs/{workflow}/{run}_qiime_{cluster}_all.log",
    shell:
        """
        exec 1> "{output}"
        echo "Primer trimming"
        echo "==============="
        cat {input.trim:q}
        printf "\n\n\n"
        echo "Denoising"
        echo "=========="
        cat {input.denoise:q}
        """


##########################
#### Taxonomy
##########################


rule qiime_taxdb_import:
    input:
        seq="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/qiime.fasta.zst",
    output:
        tmp=temp(directory("workdir/qiime_taxdb/db_regular_{source_id}/flt_{filter_id}")),
        seq="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/qiime-seqdb.qza",
        tax="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/qiime-taxonomy.qza",
    wildcard_constraints:
        source_id = "\w+",
        filter_id = "\w+",
    conda:
        config["software"]["qiime"]["conda_env"]
    log:
        "logs/taxdb/convert/db_regular_{source_id}/flt_{filter_id}-import-qiime.log",
    group:
        "taxonomy"
    shell:
        """
        exec &> {log}
        mkdir -p {output.tmp}
        zstd -dcqf {input.seq} > {output.tmp}/seq.fasta 2> {log}
        # extract taxonomic lineages from FASTA for import in QIIME
        zstd -dcqf {input} | st . --to-tsv id,desc > {output.tmp}/tax.txt
        qiime tools import \
            --type 'FeatureData[Sequence]' \
            --input-path {output.tmp}/seq.fasta \
            --output-path {output.seq}
        qiime tools import \
            --type 'FeatureData[Taxonomy]' \
            --input-format HeaderlessTSVTaxonomyFormat \
            --input-path {output.tmp}/tax.txt \
            --output-path {output.tax}
        """


rule qiime_taxdb_train_nb:
    input:
        seq="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/qiime-seqdb.qza",
        tax="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/qiime-taxonomy.qza",
    output:
        trained="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/qiime_nb.qza",
    wildcard_constraints:
        source_id = "\w+",
        filter_id = "\w+",
    cache: True
    conda:
        config["software"]["qiime"]["conda_env"]
    log:
        "logs/taxdb/convert/db_regular_{source_id}/flt_{filter_id}-train_qiime_nb.log",
    resources:
        mem_mb=40000,
        runtime=24 * 60,
    shell:
        """
        qiime feature-classifier fit-classifier-naive-bayes \
            --i-reference-reads {input.seq} \
            --i-reference-taxonomy {input.tax} \
            --o-classifier {output.trained} 2> {log}
        """


rule assign_taxonomy_qiime_sklearn:
    params:
        par=lambda wildcards: cfg.tax_config(**wildcards)["assign"],
    input:
        seq="results/{workflow}/workflow_{cluster}/{run}/{marker}__{primers}/denoised.fasta",
        db=lambda wildcards: "refdb/taxonomy/db_{preformatted}_{source_id}/flt_{filter_id}/qiime_nb.qza".format(
            **cfg.tax_config(**wildcards)
        ),
    output:
        tmp=temp(
            directory(
                "workdir/{workflow}/workflow_{cluster}/{run}/{marker}__{primers}/qiime_sklearn_{db_name}-{tax_method}"
            )
        ),
        tax="results/{workflow}/workflow_{cluster}/{run}/{marker}__{primers}/taxonomy/{db_name}-qiime_sklearn-{tax_method}.txt.gz",
    log:
        "logs/{workflow}/{run}/{marker}__{primers}/taxonomy/{cluster}_{db_name}-{tax_method}.log",
    conda:
        config["software"]["qiime"]["conda_env"]
    threads: 1  # needs a LOT of memory depending on the database
    resources:
        mem_mb=50000,
        runtime=36 * 60,
    shell:
        """
        exec &> {log}
        mkdir -p {output.tmp}
        # lowercase letters cause problems -> convert to uppercase
        # (cannot use 'st upper' because seqtool is not in conda environment,
        # therefore using awk)
        awk '/^>/ {{print($0)}}; /^[^>]/ {{print(toupper($0))}}' {input.seq} > {output.tmp}/input.fasta
        qiime tools import \
            --input-path {output.tmp}/input.fasta \
            --type 'FeatureData[Sequence]' \
            --input-format DNAFASTAFormat \
            --output-path {output.tmp}/seqs.qza

        qiime feature-classifier classify-sklearn \
            --i-classifier {input.db} \
            --i-reads {output.tmp}/seqs.qza \
            --o-classification {output.tmp}/classified.qza \
            --p-reads-per-batch 1000 \
            --p-confidence {params.par[confidence]} \
            --p-n-jobs {threads} \
            --verbose

        qiime tools export \
        --input-path {output.tmp}/classified.qza \
        --output-path {output.tmp}

        gzip -nc {output.tmp}/taxonomy.tsv > {output.tax}
        """
