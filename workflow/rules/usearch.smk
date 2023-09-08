from os.path import join, abspath, exists

usearch_cfg = config["software"]["usearch"]


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

def get_uvsnake_location():
    uvcfg = usearch_cfg["uvsnake"]
    path = uvcfg.get("path", None)
    tag = uvcfg.get("tag", None)
    commit = uvcfg.get("commit", None)
    # remote
    base_url = "https://github.com/{github}/archive".format(**uvcfg)
    if tag is not None:
        url = f"{base_url}/refs/tags/{tag}.zip"
        id_ = tag
    elif commit is not None:
        url = f"{base_url}/{commit}.zip"
        id_ = commit
    assert url is not None or path is not None, \
        "Either tag or commit or path must be defined with uvsnake source"
    return path, url, id_


def download_uvsnake(url, target_dir):
    """
    Downloads the 'uvsnake' pipeline to the working directory.
    This solution was chosen over specifying a remote
    Snakefile with github(...) due to Python modules not being
    included (see also https://github.com/snakemake/snakemake/issues/1632)
    """
    from urllib.request import urlopen
    from io import BytesIO
    import zipfile
    import sys
    import os
    import shutil
    print(f"Downloading {url}...", file=sys.stderr)
    handle = urlopen(url)
    memzip = BytesIO(handle.read())
    archive = zipfile.ZipFile(memzip)
    files = [f for f in archive.namelist() if "/workflow/" in f]
    base_dir = files[0].split("/")[0]
    parent = os.path.dirname(target_dir)
    extr_dir = os.path.join(parent, base_dir)
    os.makedirs(parent, exist_ok=True)
    archive.extractall(parent, files)
    shutil.copytree(extr_dir, target_dir, dirs_exist_ok=True)
    shutil.rmtree(extr_dir)


# paths
uvsnake_workdir = "processing/{workflow}/uvsnake/{run}_{layout}"
uvsnake_path, uvsnake_url, uvsnake_id = get_uvsnake_location()
if uvsnake_path is None:
    uvsnake_path = f"processing/_uvsnake/{uvsnake_id}"
    if not exists(uvsnake_path):
        download_uvsnake(uvsnake_url, uvsnake_path)
uvsnakefile_path = abspath(f"{uvsnake_path}/workflow/Snakefile")


rule uvsnake_gen_config:
    params:
        primer_config=lambda _: config["primers"],
        usearch_config=lambda wildcards: cfg[wildcards.workflow]["settings"]["usearch"],
        snakefile=uvsnakefile_path,
    input:
        sample_tab=lambda wildcards: "processing/{{workflow}}/input/sample_config/{technology}/{layout}/{run}/samples.tsv".format(
            **cfg.get_run_data(**wildcards)
        ),
    output:
        workdir=directory(uvsnake_workdir),
        config=uvsnake_workdir + "/config.yaml", 
        snakefile=uvsnake_workdir + "/Snakefile", 
    log:
        "logs/{workflow}/usearch/{run}_{layout}/gen_config.log",
    wildcard_constraints:
        workflow = r"[^/ ]+",
        technology = r"[^/ ]+",
        cluster_method = r"[^/ ]+",
        layout = r"(single|paired)",
        run = r"[^/ ]+",
    conda:
        "envs/basic.yaml"
    group:
        "denoise"
    script:
        "../scripts/uvsnake_gen_config.py"


rule uvsnake_prepare:
    params:
        command="prepare"
    input:
        snakefile=rules.uvsnake_gen_config.output.snakefile,
    output:
        trim_dir=directory(rules.uvsnake_gen_config.output.workdir + "/workdir/prepare_paired/2_trim"),
        stats=rules.uvsnake_gen_config.output.workdir + "/results/sample_report.tsv",
    log:
        "logs/{workflow}/usearch/{run}_{layout}/trim.log",
    conda:
        "snakemake"
    group:
        "denoise"
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
            "processing/{{workflow}}/uvsnake/{{run}}_{{layout}}/results/{primers}/{{cluster_method}}{rest}", 
            primers=cfg.primer_combinations_nomarker,
            rest=[".fasta", "_otutab.txt.gz", ".biom"],
        ),
    log:
        "logs/{workflow}/usearch/{run}_{layout}/{cluster_method}.log",
    conda:
        "snakemake"
    group:
        "denoise"
    threads:
        workflow.cores
    script:
        "../scripts/uvsnake_run.py"


rule uvsnake_copy_results:
    input:
        results=expand(
            "processing/{{workflow}}/uvsnake/{{run}}_paired/results/{primers}/{{cluster_method}}{rest}", 
            primers=cfg.primer_combinations_nomarker,
            rest=[".fasta", "_otutab.txt.gz", ".biom"],
        ),
        stats="processing/{workflow}/uvsnake/{run}_paired/results/sample_report.tsv",
    output:
        results=expand(
            "results/{{workflow}}/workflow_usearch_{{cluster_method}}/{{run}}_paired/{primers}/denoised{rest}", 
            primers=cfg.primer_combinations_flat,
            rest=[".fasta", "_otutab.txt.gz", ".biom"],
        ),
        stats="results/{workflow}/workflow_usearch_{cluster_method}/{run}_paired/sample_report.tsv",
    log:
        "logs/{workflow}/usearch/{run}_paired/copy_{cluster_method}.log",
    group:
        "denoise"
    shell:
        """
        outdir="$(dirname "{output.results[0]}")"
        for f in {input.results}; do
            name=$(basename "$f")
            out_name=denoised${{name#{wildcards.cluster_method}}}
            echo $name "$outdir"/$out_name
            cp -f "$f" "$outdir/$out_name"
        done
        cp "{input.stats}" "{output.stats}"
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
        "results/{workflow}/workflow_usearch_{cluster}/{run}_{layout}/_qc/multiqc_report.html",
    log:
        "logs/{workflow}/usearch/workflow_usearch_{cluster}/{run}_{layout}_paired_multiqc.log",
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        outdir="$(dirname "{output}")"
        multiqc -f -m fastqc -m cutadapt -o "$outdir" "{input.fastqc}" "{input.cutadapt_dir}" 2> {log}
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
        confidence=lambda w: cfg[w.workflow]["taxonomy"][w.marker][(w.db_name, w.tax_method)]["assign"]["confidence"],
        program=lambda w: cfg[w.workflow]["taxonomy"][w.marker][(w.db_name, w.tax_method)]["assign"]["program"],
        usearch_bin=config["software"]["usearch"]["binary"],
        maxaccepts=1,
        maxrejects=1,
    input:
        fa="results/{workflow}/workflow_{cluster}/{run}_{layout}/{marker}__{primers}/denoised.fasta",
        db=lambda w: "refdb/taxonomy/db_{source[preformatted]}_{source[source_id]}/flt_{filter_id}/utax.fasta".format(
            **cfg[w.workflow]["taxonomy"][w.marker][(w.db_name, w.tax_method)]
        ),
    output:
        taxtab="results/{workflow}/workflow_{cluster}/{run}_{layout}/{marker}__{primers}/taxonomy/{db_name}-sintax_usearch-{tax_method}.txt.gz",
        sintax="results/{workflow}/workflow_{cluster}/{run}_{layout}/{marker}__{primers}/taxonomy/sintax/{db_name}-sintax_usearch-{tax_method}.txt.gz",
    log:
        "logs/{workflow}/workflow_{cluster}/{run}_{layout}/{marker}__{primers}/taxonomy_sintax/{db_name}-{tax_method}.log",
    group:
        "taxonomy"
    conda:
        "envs/vsearch.yaml"
    threads:
        lambda w: workflow.cores \
        if cfg[w.workflow]["taxonomy"][w.marker][(w.db_name, w.tax_method)]["assign"]["program"] == "vsearch" \
        else 1
 