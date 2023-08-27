
localrules:
    qiime_make_manifest,
    qiime_multiqc,
    qiime_stats_paired,


rule qiime_make_manifest:
    # This is a bit complicated, but serves to obtain the config and output
    # directory paths, inferring the technology wildcard (which is otherwise
    # not known)
    input:
        run_config=lambda wildcards: expand_runs(
            "input/sample_config/{technology}/{layout}/{run}",
            pool=cfg[wildcards.workflow]["settings"]["pool_raw"],
            **wildcards
        )[0],
        fq=lambda wildcards: expand_runs(
            "processing/{workflow}/input/{technology}/{layout}/{run}",
            pool=cfg[wildcards.workflow]["settings"]["pool_raw"],
            **wildcards
        )[0],
    output:
        manifest="processing/{workflow}/qiime/manifest/{run}/{layout}/manifest.txt",
    log:
        "logs/{workflow}/qiime/{run}/{layout}/make_manifest.log",
    script:
        "../scripts/qiime_make_manifest.py"


rule qiime_import:
    params:
        type=lambda wildcards: "SequencesWithQuality"
          if wildcards.layout == "single"
          else "PairedEndSequencesWithQuality",
        format=lambda wildcards: "SingleEndFastqManifestPhred33V2"
          if wildcards.layout == "single"
          else "PairedEndFastqManifestPhred33V2",
    input:
        manifest="processing/{workflow}/qiime/manifest/{run}/{layout}/manifest.txt",
    output:
        "processing/{workflow}/qiime/{run}/{layout}/demux.qza",
    log:
        "logs/{workflow}/qiime/manifest/{run}/{layout}/import.log",
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
        min_len=lambda w: cfg[w.workflow]["settings"]["filter"]["min_length"],
    input:
        yaml="processing/primers/primers.yaml",
        demux="processing/{workflow}/qiime/{run}/paired/demux.qza",
    output:
        qza="processing/{workflow}/qiime/{run}/paired/{marker}__{f_primer}...{r_primer}/trim.qza",
    log:
        "logs/{workflow}/qiime/{run}/paired/{marker}__{f_primer}...{r_primer}/trim.log",
    group:
        "prepare"
    conda:
        config["software"]["qiime"]["conda_env"]
    threads: workflow.cores
    resources:
        runtime=12 * 60,
    script:
        "../scripts/qiime_trim_paired.py"


rule qiime_denoise_paired:
    params:
        par=lambda w: cfg[w.workflow]["settings"]["dada2"],
    input:
        trim="processing/{workflow}/qiime/{run}/paired/{primers}/trim.qza",
    output:
        denoised0="processing/{workflow}/qiime/dada2/{run}/paired/{primers}/dada2.qza",
        tab0="processing/{workflow}/qiime/dada2/{run}/paired/{primers}/dada2_tab.qza",
        stats="processing/{workflow}/qiime/dada2/{run}/paired/{primers}/dada2_stats.qza",
    log:
        "logs/{workflow}/qiime/dada2/{run}/paired/{primers}/dada2.log",
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
    qiime_denoised_export > tsv_to_biom


rule qiime_denoised_export:
    input:
        denoised0="processing/{workflow}/qiime/dada2/{run}/{layout}/{primers}/dada2.qza",
        tab0="processing/{workflow}/qiime/dada2/{run}/{layout}/{primers}/dada2_tab.qza",
    output:
        denoised="results/{workflow}/workflow_qiime_dada2/{run}_{layout}/{primers}/denoised.fasta",
        tab="results/{workflow}/workflow_qiime_dada2/{run}_{layout}/{primers}/denoised_otutab.txt.gz",
        biom_json="results/{workflow}/workflow_qiime_dada2/{run}_{layout}/{primers}/denoised.biom",
        biom_hdf5="results/{workflow}/workflow_qiime_dada2/{run}_{layout}/{primers}/denoised.hdf5.biom",
        tmp=temp(
            directory("processing/{workflow}/qiime/dada2/{run}/{layout}/{primers}/denoised_convert_tmp")
        ),
    log:
        "logs/{workflow}/qiime/dada2/{run}/{layout}/{primers}/denoised_convert.log",
    conda:
        config["software"]["qiime"]["conda_env"]
    group:
        "denoise"
    shell:
        """
        # export table
        mkdir -p {output.tmp}
        tab={output.tab}
        tab=${{tab%.gz}}
        qiime tools export \
            --input-path {input.tab0} \
            --output-path {output.tmp} &> {log}
        mv {output.tmp}/feature-table.biom {output.biom_hdf5}
        biom convert -i {output.biom_hdf5}  \
            -o {output.biom_json} \
            --to-json --table-type "OTU table" &> {log}
        biom convert -i {output.biom_hdf5}  \
            -o $tab \
            --to-tsv --table-type "OTU table" &> {log}
        sed -i '1,1d' $tab
        gzip -nf $tab

        # export seqs
        qiime tools export \
            --input-path {input.denoised0} \
            --output-path {output.tmp} &> {log}
        cat {output.tmp}/dna-sequences.fasta > {output.denoised}
        """


##########################
#### Taxonomy
##########################


rule qiime_taxdb_import:
    input:
        seq="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/qiime.fasta.zst",
    output:
        tmp=temp(directory("processing/qiime_taxdb/db_regular_{source_id}/flt_{filter_id}")),
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
        mkdir -p {output.tmp}
        zstd -dcqf {input.seq} > {output.tmp}/seq.fasta 2> {log}
        # extract taxonomic lineages from FASTA for import in QIIME
        zstd -dcqf {input} | st . --to-tsv id,desc > {output.tmp}/tax.txt
        qiime tools import \
            --type 'FeatureData[Sequence]' \
            --input-path {output.tmp}/seq.fasta \
            --output-path {output.seq} &> {log}
        qiime tools import \
            --type 'FeatureData[Taxonomy]' \
            --input-format HeaderlessTSVTaxonomyFormat \
            --input-path {output.tmp}/tax.txt \
            --output-path {output.tax} &> {log}
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
            --o-classifier {output.trained}
        """


rule assign_taxonomy_qiime_sklearn:
    params:
        par=lambda w: cfg[w.workflow]["taxonomy"][w.marker][(w.db_name, w.tax_method)]["assign"],
    input:
        seq="results/{workflow}/workflow_{cluster}/{run}_{layout}/{marker}__{primers}/denoised.fasta",
        db=lambda w: "refdb/taxonomy/db_{preformatted}_{source_id}/flt_{filter_id}/qiime_nb.qza".format(
            **cfg[w.workflow]["taxonomy"][w.marker][(w.db_name, w.tax_method)]
        ),
    output:
        tmp=temp(
            directory(
                "processing/{workflow}/workflow_{cluster}/qiime_sklearn/{run}/{layout}/{marker}__{primers}/{db_name}-{tax_method}"
            )
        ),
        tax="results/{workflow}/workflow_{cluster}/{run}_{layout}/{marker}__{primers}/taxonomy/{db_name}-qiime_sklearn-{tax_method}.txt.gz",
    log:
        "logs/processing/{workflow}/workflow_{cluster}/qiime_sklearn/{run}/{layout}/{marker}__{primers}/{db_name}-{tax_method}.log",
    conda:
        config["software"]["qiime"]["conda_env"]
    threads: 1  # needs a LOT of memory depending on the database
    resources:
        mem_mb=50000,
        runtime=36 * 60,
    shell:
        """
        mkdir -p {output.tmp}
        # lowercase letters cause problems -> convert to uppercase
        # (cannot use 'st upper' because seqtool is not in conda environment,
        # therefore using awk)
        awk '/^>/ {{print($0)}}; /^[^>]/ {{print(toupper($0))}}' {input.seq} > {output.tmp}/input.fasta
        qiime tools import \
            --input-path {output.tmp}/input.fasta \
            --type 'FeatureData[Sequence]' \
            --input-format DNAFASTAFormat \
            --output-path {output.tmp}/seqs.qza &> {log}

        qiime feature-classifier classify-sklearn \
            --i-classifier {input.db} \
            --i-reads {output.tmp}/seqs.qza \
            --o-classification {output.tmp}/classified.qza \
            --p-reads-per-batch 1000 \
            --p-confidence {params.par[confidence]} \
            --p-n-jobs {threads} \
            --verbose &> {log}

        qiime tools export \
        --input-path {output.tmp}/classified.qza \
        --output-path {output.tmp} &> {log}

        gzip -nc {output.tmp}/taxonomy.tsv > {output.tax}
        """


##########################
#### QC
##########################


rule qiime_multiqc:
    input:
        rules.multiqc_fastqc.output,
    output:
        "results/{workflow}/validation/multiqc_qiime/multiqc_report.html",
    log:
        "logs/{workflow}/multiqc_qiime.log",
    shell:
        "ln -srf {input} {output} 2> {log}"


# TODO: not implemented
rule qiime_stats_paired:
    output:
        touch("results/{workflow}/workflow_qiime_{cluster}/{run}_paired/sample_report.tsv"),
    log:
        "logs/{workflow}/workflow_qiime_{cluster}/{run}_paired/sample_report.log",
