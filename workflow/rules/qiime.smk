from os.path import basename
from copy import deepcopy
from sys import stderr


localrules:
    qiime_make_manifest_paired,
    # qiime_multiqc_paired,
    qiime_stats_paired,


rule qiime_make_manifest_paired:
    params:
        samples=cfg.sample_names["paired"],
    input:
        expand(
            "input/grouped/paired/{sample}/{sample}_R{read}.fastq.gz",
            sample=cfg.sample_names["paired"],
            read=[1, 2],
        ),
    output:
        "processing/qiime/paired/manifest.txt",
    log:
        "logs/qiime/make_manifest_paired.log",
    script:
        "../scripts/qiime_manifest_paired.py"


rule qiime_import:
    params:
        type=lambda w: "SequencesWithQuality"
        if w.strategy == "single"
        else "PairedEndSequencesWithQuality",
        format=lambda w: "SingleEndFastqManifestPhred33V2"
        if w.strategy == "single"
        else "PairedEndFastqManifestPhred33V2",
    input:
        manifest="processing/qiime/{strategy}/manifest.txt",
        seq=lambda w: expand(
            "input/grouped/paired/{sample}/{sample}_R{read}.fastq.gz",
            sample=cfg.sample_names[w.strategy],
            read=[1, 2] if w.strategy == "paired" else [1],
        ),
    output:
        "processing/{name}/qiime/{strategy}/demux.qza",
    log:
        "logs/{name}/qiime/import_{strategy}.log",
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
        err_rate=lambda w: cfg[w.name]["settings"]["primers"]["trim_settings"]["max_error_rate"],
        min_len=lambda w: cfg[w.name]["settings"]["filter"]["min_length"],
    input:
        yaml="processing/primers/primers.yaml",
        demux="processing/{name}/qiime/paired/demux.qza",
    output:
        qza="processing/{name}/qiime/paired/{marker}__{f_primer}...{r_primer}/trim.qza",
    log:
        "logs/{name}/qiime/paired/{marker}__{f_primer}...{r_primer}/trim.log",
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
        par=lambda w: cfg[w.name]["settings"]["dada2"],
    input:
        trim="processing/{name}/qiime/paired/{primers}/trim.qza",
    output:
        denoised0="processing/{name}/qiime/paired/{primers}/dada2.qza",
        tab0="processing/{name}/qiime/paired/{primers}/dada2_tab.qza",
        stats="processing/{name}/qiime/paired/{primers}/dada2_stats.qza",
    log:
        "logs/{name}/qiime/paired/{primers}/dada2.log",
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


rule qiime_denoised_convert:
    input:
        denoised0="processing/{name}/qiime/paired/{primers}/dada2.qza",
        tab0="processing/{name}/qiime/paired/{primers}/dada2_tab.qza",
    output:
        denoised="results/{name}/pipeline_qiime_dada2/{primers}/paired/denoised.fasta",
        tab="results/{name}/pipeline_qiime_dada2/{primers}/paired/denoised_otutab.txt.gz",
        tmp=temp(
            directory("processing/{name}/qiime/paired/{primers}/denoised_convert_tmp")
        ),
    log:
        "logs/{name}/qiime/paired/{primers}/denoised_convert.log",
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
        biom convert -i {output.tmp}/feature-table.biom  \
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


rule taxdb_extract_taxonomy:
    input:
        "refdb/taxonomy/{db}/filtered/{defined}.fasta.zst",
    output:
        "refdb/taxonomy/{db}/filtered/{defined}.txt.zst",
    conda:
        "envs/basic.yaml"
    group:
        "taxonomy"
    shell:
        """
        zstd -dcqf {input} | st . --to-tsv id,desc | zstd -cq > {output}
        """


rule qiime_taxdb_import:
    input:
        seq="refdb/taxonomy/{db}/filtered/{defined}.fasta.zst",
        tax="refdb/taxonomy/{db}/filtered/{defined}.txt.zst",
    output:
        tmp=temp(directory("processing/qiime_taxdb/{db}/{defined}")),
        seq="refdb/taxonomy/{db}/formatted/qiime_qza/{defined}.qza",
        tax="refdb/taxonomy/{db}/formatted/qiime_qza/{defined}-taxonomy.qza",
    conda:
        config["software"]["qiime"]["conda_env"]
    log:
        "logs/taxdb/convert/{db}/qiime_import/{defined}.log",
    group:
        "taxonomy"
    shell:
        """
        mkdir -p {output.tmp}
        zstd -dcqf {input.seq} > {output.tmp}/seq.fasta 2> {log}
        zstd -dcqf {input.tax} > {output.tmp}/tax.txt 2> {log}
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
        seq="refdb/taxonomy/{db}/formatted/qiime_qza/{defined}.qza",
        tax="refdb/taxonomy/{db}/formatted/qiime_qza/{defined}-taxonomy.qza",
    output:
        "refdb/taxonomy/{db}/formatted/qiime_nb/{defined}.qza",
    conda:
        config["software"]["qiime"]["conda_env"]
    log:
        "logs/taxdb/convert/{db}/qiime_nb/{defined}.log",
    resources:
        mem_mb=40000,
        runtime=24 * 60,
    shell:
        """
        qiime feature-classifier fit-classifier-naive-bayes \
            --i-reference-reads {input.seq} \
            --i-reference-taxonomy {input.tax} \
            --o-classifier {output}
        """


rule assign_taxonomy_qiime_sklearn:
    params:
        par=lambda w: cfg[w.name]["taxonomy"][w.marker][(w.db_name, w.tax_method)],
    input:
        seq="results/{name}/{pipeline}/{marker}__{primers}/{strategy}/denoised.fasta",
        db=lambda w: "refdb/taxonomy/{{marker}}/{db_name}/formatted/qiime_nb/{defined}.qza".format(
            **cfg[w.name]["taxonomy"][w.marker][(w.db_name, w.tax_method)]
        ),
    output:
        tmp=temp(
            directory(
                "processing/{name}/{pipeline}/{marker}__{primers}/{strategy}/taxonomy/{db_name}-{tax_method}/denoised_sklearn"
            )
        ),
        tax="results/{name}/{pipeline}/{marker}__{primers}/{strategy}/taxonomy/{db_name}-qiime_sklearn-{tax_method}.txt.gz",
    log:
        "logs/{name}/{pipeline}/{strategy}/taxonomy_sklearn/{marker}__{primers}/{db_name}-{tax_method}.log",
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


rule qiime_multiqc_paired:
    input:
        rules.multiqc_fastqc.output,
    output:
        "results/{name}/pipeline_qiime_{method}/_validation/multiqc/multiqc_report.html",
    shell:
        "ln -srf {input} {output}"


rule qiime_stats_paired:
    output:
        "results/{name}/pipeline_qiime_{method}/_validation/sample_report.tsv",
    log:
        "logs/{name}/qiime/sample_report_{method}.log",
    run:
        # TODO: not implemented
        with open(output[0], "w") as out:
            pass
