from os.path import join, abspath, exists
from lib import get_repo_location, download_repo

cfg.pipeline_capabilities["usearch"] = [
    ("illumina", "paired"),
]


# def get_uvsnake_source():
#     uvcfg = usearch_cfg["uvsnake"]
#     path = uvcfg.get("path", None)
#     if path is not None:
#         return join(path, "workflow", "Snakefile")
#     return github(
#         uvcfg["github"],
#         path="workflow/Snakefile",
#         tag=uvcfg.get("tag", None),
#         commit=uvcfg.get("commit", None)
#     )


# paths
uvsnake_workdir = "workdir/{workflow}/{run}_paired"
uvsnake_workdir_q = "workdir/{{workflow}}/{{run}}_paired"
uvsnake_path, uvsnake_url, uvsnake_id = get_repo_location(**config["software"]["uvsnake"]["repo"])
if uvsnake_path is None:
    uvsnake_path = f"workdir/_uvsnake/{uvsnake_id}"
    if not exists(uvsnake_path):
        download_repo(uvsnake_url, uvsnake_path)
uvsnakefile_path = abspath(f"{uvsnake_path}/workflow/Snakefile")


rule uvsnake_gen_config:
    params:
        primer_config=lambda _: config["primers"],
        usearch_config=lambda wildcards: cfg[wildcards.workflow]["settings"]["usearch"],
        snakefile=uvsnakefile_path,
    input:
        sample_tab=lambda wildcards: "workdir/{{workflow}}/input/sample_config/{technology}/paired/{run}/samples.tsv".format(
            **cfg.get_run_data(layout="paired", **wildcards)
        ),
    output:
        config=uvsnake_workdir + "/config.yaml", 
        snakefile=uvsnake_workdir + "/Snakefile", 
    log:
        "logs/{workflow}/{run}_paired/uvsnake_gen_config.log",
    conda:
        "envs/basic.yaml"
    group:
        "cluster"
    script:
        "../scripts/uvsnake_gen_config.py"


rule uvsnake_prepare:
    params:
        command="prepare"
    input:
        snakefile=rules.uvsnake_gen_config.output.snakefile,
    output:
        trim_dir=directory(uvsnake_workdir + "/workdir/prepare_paired/2_trim"),
        stats=uvsnake_workdir + "/results/sample_report.tsv",
    log:
        "logs/{workflow}/{run}_paired/uvsnake_trim.log",
    conda:
        "snakemake"
    group:
        "cluster"
    threads:
        workflow.cores
    script:
        "../scripts/uvsnake_run.py"


rule uvsnake_cluster:
    params:
        command=lambda wildcards: wildcards.cluster_method,
    input:
        snakefile=rules.uvsnake_gen_config.output.snakefile,
        _trim=rules.uvsnake_prepare.output.trim_dir,
    output:
        results=expand(
            uvsnake_workdir_q + "/results/{primers}/{{cluster_method}}{rest}", 
            primers=cfg.primer_combinations_nomarker,
            rest=[".fasta", "_otutab.txt.gz", ".biom"],
        ),
        combined_log=uvsnake_workdir + "/logs/cluster_{cluster_method}_all.log",
    log:
        "logs/{workflow}/{run}_paired/{cluster_method}_uvsnake.log",
    conda:
        "snakemake"
    group:
        "cluster"
    threads:
        workflow.cores
    script:
        "../scripts/uvsnake_run.py"


rule uvsnake_copy_results:
    input:
        results=rules.uvsnake_cluster.output.results,
        stats=uvsnake_workdir + "/results/sample_report.tsv",
        log=rules.uvsnake_cluster.output.combined_log,
    output:
        results=expand(
            "results/{{workflow}}/workflow_usearch_{{cluster_method}}/{{run}}_paired/{primers}/{what}", 
            primers=cfg.primer_combinations_flat,
            what=["clusters.fasta", "otutab.txt.gz", "otutab.biom"],
        ),
        stats="results/{workflow}/workflow_usearch_{cluster_method}/{run}_paired/sample_report.tsv",
        log="logs/{workflow}/{run}_paired_usearch_{cluster_method}_all.log",
    log:
        "logs/{workflow}/{run}_paired/uvsearch_copy_{cluster_method}.log",
    group:
        "cluster"
    shell:
        """
        indir="$(dirname "{input.results[0]}")"
        outdir="$(dirname "{output.results[0]}")"
        cp -f "$indir/{wildcards.cluster_method}.fasta" "$outdir/clusters.fasta"
        cp -f "$indir/{wildcards.cluster_method}_otutab.txt.gz" "$outdir/otutab.txt.gz"
        cp -f "$indir/{wildcards.cluster_method}.biom" "$outdir/otutab.biom"
        cp -f "{input.stats}" "{output.stats}"
        cp -f "{input.log}" "{output.log}"
        """

ruleorder:
    uvsnake_copy_results > otutab_to_biom

# TODO: priority?
ruleorder:
    usearch_multiqc > multiqc_fastqc_pooled_workdir > multiqc_link_workdir



##########################
#### QC
##########################


rule usearch_multiqc:
    input:
        fastqc=lambda wildcards: rules.fastqc.output.qc_dir.format(
            demux="demux",
            **cfg.get_run_data(**wildcards)
        ),
        cutadapt_dir=rules.uvsnake_prepare.output.trim_dir,
    output:
        "results/{workflow}/workflow_usearch_{cluster}/{run}_paired/_qc/multiqc_report.html",
    log:
        "logs/{workflow}/{run}_paired/usearch_{cluster}_paired_multiqc.log",
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        outdir="$(dirname "{output}")"
        multiqc -f -m fastqc -m cutadapt -o "$outdir" "{input.fastqc}" "{input.cutadapt_dir}" &> {log}
        """


##########################
#### Taxonomy
##########################


rule convert_taxdb_utax:
    input:
        db="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/qiime.fasta.zst"
    output:
        db="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/utax.fasta.zst"
    log:
        "logs/taxdb/filter/db_regular_{source_id}/flt_{filter_id}/convert_utax.log",
    wildcard_constraints:
        source_id = "\w+",
        filter_id = "\w+",
    conda:
        "envs/taxonomy.yaml"
    group:
        "taxonomy"
    script:
        "../scripts/taxonomy/convert_utax.py"



# We directly reuse the sintax assignment rule from UVSnake by specifying
# the pipeline as module
module uvsnake:
    snakefile: uvsnakefile_path
    # provide only a minimal "dummy" configuration, since the 
    # rule params/input/output is overridden anyway
    config: {"sintax": None}
    skip_validation: True 


use rule assign_taxonomy_sintax from uvsnake with:
    params:
        confidence=lambda wildcards: cfg.tax_config(**wildcards)["assign"]["confidence"],
        program=lambda wildcards: cfg.tax_config(**wildcards)["assign"]["program"],
        usearch_bin=config["software"]["usearch"]["binary"],
        maxaccepts=1,
        maxrejects=1,
    input:
        fa="results/{workflow}/workflow_{cluster}/{run}/{marker}__{primers}/clusters.fasta",
        db=lambda wildcards: "refdb/taxonomy/db_{source[preformatted]}_{source[source_id]}/flt_{filter_id}/utax.fasta".format(
            **cfg.tax_config(**wildcards)
        ),
    output:
        taxtab="results/{workflow}/workflow_{cluster}/{run}/{marker}__{primers}/taxonomy/{db_name}-sintax_usearch-{tax_method}.txt.gz",
        sintax="results/{workflow}/workflow_{cluster}/{run}/{marker}__{primers}/taxonomy/sintax/{db_name}-sintax_usearch-{tax_method}.txt.gz",
    log:
        "logs/{workflow}/{run}/{marker}__{primers}/taxonomy/{cluster}_{db_name}-{tax_method}.log",
    group:
        "taxonomy"
    conda:
        "envs/vsearch.yaml"
    threads:
        lambda wildcards: workflow.cores \
        if cfg.tax_config(**wildcards)["assign"]["program"] == "vsearch" \
        else 1
 