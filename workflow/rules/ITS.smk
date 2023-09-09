
localrules:
    clean_itsx,


rule itsx:
    params:
        par=config["ITSx"],
    input:
        fa="results/{workflow}/workflow_{cluster}/{run}/ITS__{primers}/clusters.fasta",
    output:
        pos="results/{workflow}/workflow_{cluster}/{run}/ITS__{primers}/ITSx/out.positions.txt",
    log:
        "logs/{workflow}/{run}/ITS__{primers}/{cluster}_ITSx.log",
    group:
        "ITS"
    conda:
        "envs/itsx.yaml"
    threads: min(workflow.cores, 8)  # too many does not make sense
    shell:
        """
        pos="{output.pos}"
        ITSx -i {input} -o "${{pos%.positions.txt}}" \
          -t {params.par[organism_groups]} \
          -E {params.par[e-value]} \
          --allow_single_domain {params.par[allow_single_domain]} \
          --complement {params.par[complement]} \
          --heuristics {params.par[heuristics]} \
          --graphical {params.par[graphical]} \
          --fasta {params.par[fasta]} \
          --preserve {params.par[preserve]} \
          --save_regions {params.par[save_regions]} \
          --partial {params.par[partial]} \
          --cpu {threads} \
          2> {log}
        """

rule clean_itsx:
    shell:
        "rm -Rf results/*/workflow_*/*/*/ITSx"
