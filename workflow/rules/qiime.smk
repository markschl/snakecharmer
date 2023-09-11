
cfg.cluster_capabilities["qiime"] = [
    ("illumina", "single"),
    ("illumina", "single.rev"),
    ("illumina", "paired"),
]
cfg.taxonomy_formats["qiime_sklearn"] = "qiime_nb"


localrules:
    qiime_make_manifest,
    qiime_combine_sample_reports,


##########################
#### Prepare
##########################


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
          if wildcards.layout.startswith("single")
          else "PairedEndSequencesWithQuality",
        format=lambda wildcards: "SingleEndFastqManifestPhred33V2"
          if wildcards.layout.startswith("single")
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


rule qiime_trim_single:
    params:
        err_rate=lambda w: cfg[w.workflow]["settings"]["primers"]["trim_settings"]["max_error_rate"],
        min_length=lambda w: cfg[w.workflow]["settings"]["primers"]["trim_settings"]["min_length"],
        rev_read=lambda wildcards: wildcards.suffix == ".rev",
    input:
        yaml="workdir/primers/primers.yaml",
        demux="workdir/{workflow}/{run}_single{suffix}/demux.qza",
    output:
        qza="workdir/{workflow}/{run}_single{suffix}/{marker}__{f_primer}...{r_primer}/trim.qza",
    log:
        "logs/{workflow}/{run}_single{suffix}/{marker}__{f_primer}...{r_primer}/qiime_trim.log",
    wildcard_constraints:
        suffix = r"(|\.rev)",
    group:
        "prepare"
    conda:
        config["software"]["qiime"]["conda_env"]
    threads: workflow.cores
    resources:
        runtime=12 * 60,
    script:
        "../scripts/qiime_trim_single.py"


rule qiime_trim_paired:
    params:
        err_rate=lambda w: cfg[w.workflow]["settings"]["primers"]["trim_settings"]["max_error_rate"],
        min_length=lambda w: cfg[w.workflow]["settings"]["primers"]["trim_settings"]["min_length"],
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


rule qiime_dada2:
    params:
        par=lambda w: cfg[w.workflow]["settings"]["dada2"],
    input:
        trim="workdir/{workflow}/{run}_{layout}/{primers}/trim.qza",
    output:
        asvs="results/{workflow}/workflow_qiime_dada2/{run}_{layout}/{primers}/clusters.qza",
        tab="results/{workflow}/workflow_qiime_dada2/{run}_{layout}/{primers}/otutab.qza",
        stats="workdir/{workflow}/{run}_{layout}/{primers}/dada2_stats.qza",
    log:
        "logs/{workflow}/{run}_{layout}/{primers}/dada2_qiime.log",
    conda:
        config["software"]["qiime"]["conda_env"]
    group:
        "cluster"
    threads: workflow.cores
    resources:
        mem_mb=30000,
        runtime=36 * 60,
    script:
        "../scripts/qiime_dada2.py"



##########################
#### After clustering
##########################


ruleorder:
    qiime_export > otutab_to_biom > biom_to_hdf5


rule qiime_export:
    input:
        clustered="results/{workflow}/workflow_qiime_{cluster}/{run}/{primers}/clusters.qza",
        tab="results/{workflow}/workflow_qiime_{cluster}/{run}/{primers}/otutab.qza",
        stats="workdir/{workflow}/{run}/{primers}/dada2_stats.qza",
    output:
        # directory("results/{workflow}/workflow_qiime_{cluster}/{run}/{primers}"),
        clustered="results/{workflow}/workflow_qiime_{cluster}/{run}/{primers}/clusters.fasta",
        tab="results/{workflow}/workflow_qiime_{cluster}/{run}/{primers}/otutab.txt.gz",
        biom_json="results/{workflow}/workflow_qiime_{cluster}/{run}/{primers}/otutab.biom",
        biom_hdf5="results/{workflow}/workflow_qiime_{cluster}/{run}/{primers}/otutab.hdf5.biom",
        stats=temp("results/{workflow}/workflow_qiime_{cluster}/{run}/{primers}/sample_report.tsv"),
        tmp=temp(
            directory("workdir/{workflow}/{run}/{primers}/{cluster}_export_tmp")
        ),
    log:
        "logs/{workflow}/{run}/{primers}/{cluster}_qiime_export.log",
    wildcard_constraints:
        layout = r"(paired|single(_rev)?)",
    conda:
        config["software"]["qiime"]["conda_env"]
    group:
        "cluster"
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
            --input-path {input.clustered} \
            --output-path {output.tmp}
        cat {output.tmp}/dna-sequences.fasta > {output.clustered}

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
        cluster=expand(
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
        cat {input.cluster:q}
        """


##########################
#### Taxonomy
##########################


rule qiime_taxdb_train_nb:
    input:
        params="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/qiime_nb/conversion_config_{cnv_id}.yaml",
        seq="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/db.fasta.zst",
    output:
        "refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/qiime_nb/cnv_{cnv_id}.zst",
    wildcard_constraints:
        source_id = "\w+",
        filter_id = "\w+",
        cnv_id = "\w+",
    cache: True
    conda:
        config["software"]["qiime"]["conda_env"]
    log:
        "logs/taxonomy/db_regular_{source_id}/flt_{filter_id}/qiime_train_nb-{cnv_id}.log",
    resources:
        mem_mb=40000,
        runtime=24 * 60,
    shell:
        """
        exec &> {log}
        par=$(cat {input.params)
        if [ "$par" != "{}" ]; then
            echo "Unknown QIIME naive-bayes classifier config: $par"
            exit 1
        fi
        mkdir -p {output.tmp}
        zstd -dcqf {input.seq} > {output.tmp}/seq.fasta 2> {log}
        # extract taxonomic lineages from FASTA for import in QIIME
        zstd -dcqf {input} | st . --to-tsv id,desc > {output.tmp}/tax.txt
        qiime tools import \
            --type 'FeatureData[Sequence]' \
            --input-path {output.tmp}/seq.fasta \
            --output-path {output.tmp}/seq.qza
        qiime tools import \
            --type 'FeatureData[Taxonomy]' \
            --input-format HeaderlessTSVTaxonomyFormat \
            --input-path {output.tmp}/tax.txt \
            --output-path {output.tmp}/tax.qza
        qiime feature-classifier fit-classifier-naive-bayes \
            --i-reference-reads {output.tmp}/seq.qza \
            --i-reference-taxonomy {output.tmp}/tax.qza \
            --o-classifier {output.tmp}/trained.qza
        zstd -dcq {output.tmp}/trained.qza > {output.trained}
        """


rule qiime_classify_sklearn:
    params:
        par=lambda wildcards: cfg.tax_config(**wildcards)["assign"],
    input:
        seq="results/{workflow}/workflow_{cluster}/{run}/{marker}__{primers}/clusters.fasta",
        db=lambda wildcards: get_refdb_path(format="qiime_nb", **wildcards)
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
