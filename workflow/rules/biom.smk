
rule otutab_to_biom:
    input:
        otutab="results/{workflow}/workflow_{cluster}/{run}/{sub_path}/otutab.txt.gz",
    output:
        tmp_tab=temp("results/{workflow}/workflow_{cluster}/{run}/{sub_path}/otutab_tmp.txt"),
        biom="results/{workflow}/workflow_{cluster}/{run}/{sub_path}/otutab.biom",
    log:
        "logs/{workflow}/{run}/{sub_path}/{cluster}_otutab_to_biom.log",
    group:
        "cluster"
    priority: -100
    conda:
        "envs/biom.yaml"
    shell:
        """
        gzip -dc {input.otutab} > {output.tmp_tab}
        biom convert -i {input.otutab} \
          -o {output.biom} \
          --table-type 'OTU table' --to-json &> {log}
        """


rule biom_to_hdf5:
    input:
        biom="results/{workflow}/workflow_{cluster}/{run}/{sub_path}/otutab.biom",
    output:
        biom_hdf5="results/{workflow}/workflow_{cluster}/{run}/{sub_path}/otutab.hdf5.biom",
    log:
        "logs/{workflow}/{run}/{sub_path}/{cluster}_biom_to_hdf5.log",
    group:
        "cluster"
    priority: -100
    conda:
        "envs/biom.yaml"
    shell:
        """
        biom convert -i {input.biom}  \
          -o {output.biom_hdf5} \
          --table-type "OTU table" --to-hdf5 &> {log}
        """
