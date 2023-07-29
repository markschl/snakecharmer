import os
from os.path import *
import shutil
import lib
from collections import OrderedDict


localrules:
    dump_samples_yaml,
    dump_samples_tsv,
    dump_config,
    link_samples,
    link_data_dirs,



#### Start ####

cfg = lib.Config(config)


# Assemble dict with paths to link:
#  dict keys are the unique samples
#    - source file path
#    - nested output dir (strategy -> sample [ -> "_multiple" if duplicate names]
#    - unique output name [with optional suffix to make name unique]
link_paths = []

for strategy, samples in cfg.samples.items():
    for sample_name, paths in samples.items():
        sample_dir = join(strategy, sample_name)
        if len(paths) == 1:
            # sample name is unique -> link files directly
            # **note**: R2 will become R1 here if only R2 was supplied
            link_paths.append((
                (sample_name, sample_name),
                [(p, sample_dir, "{}_R{}.fastq.gz".format(sample_name, i + 1))
                for i, p in enumerate(paths[0])]
            ))
        else:
            # sample name is duplicated -> add a suffix number
            # and put them into the "_multiple" directory, from where they can be
            # combined
            for sample_i, read_paths in enumerate(paths):
                unique_name = "{}_{}".format(sample_name, sample_i + 1)
                link_paths.append((
                    (sample_name, unique_name),
                    [(p, join(sample_dir, '_multiple'),
                      "{}_R{}.fastq.gz".format(unique_name, i + 1))
                    for i, p in enumerate(read_paths)]
                ))

# flat version of link_paths without grouping by sample
link_paths_flat = [p for _, paths in link_paths 
                   for p in paths]


if len(link_paths_flat) == 0:
    raise Exception('No sample files found! Are the input directories/patterns correctly specified?')


#### Configuration ####


rule dump_samples_yaml:
    params:
        samples = cfg.samples,
    output:
        yaml="results/samples.yaml",
    log:
        "logs/dump_sampes.log"
    script:
        "../scripts/dump_samples_yaml.py"


rule dump_samples_tsv:
    params:
        link_paths = link_paths,
    output:
        tsv="results/samples.tsv",
    log:
        "logs/dump_sampes.log"
    script:
        "../scripts/dump_samples_tsv.py"


rule dump_config:
    params:
        config = lambda wildcards: cfg.pipelines[wildcards.name]
    output:
        "results/{name}/config.yaml",
    log:
        "logs/{name}/dump_config.log"
    script:
        "../scripts/dump_config.py"


#### Sample handling ####


rule prepare_primers:
    params:
        primers = cfg.primers,
    output:
        yaml='processing/primers/primers.yaml',
    log:
        "logs/prepare_primers.log",
    conda:
        "envs/consensus.yaml"
    script:
        "../scripts/prepare_primers.py"


rule link_data_dirs:
    params:
        workflow_cfg = cfg.pipelines
    output:
        [
            directory("results/{name}/data".format(**p))
            for p in cfg.pipelines.values()
            if p["is_simple"]
        ],
    log:
        "logs/link_data_dirs.log",
    script:
        "../scripts/link_data_dirs.py"


rule link_samples:
    params:
        link_paths = link_paths_flat
    output:
        grouped = [os.path.join("input", "grouped", path, name)
                   for _, path, name in link_paths_flat],
        unique = [os.path.join("input", "unique_samples", name)
                  for _, path, name in link_paths_flat]
    log:
        "logs/link_samples.log",
    script:
        "../scripts/link_samples.py"



rule combine_multiple_samples:
    input:
        lambda w: expand(
            "input/grouped/{{strategy}}/{{sample}}/_multiple/{{sample}}_{file_i}_R{{read}}.fastq.gz",
            file_i=[i + 1 for i in range(len(cfg.samples[w.strategy][w.sample]))]
        ),
    output:
        "input/grouped/{strategy}/{sample}/{sample}_R{read}.fastq.gz",
    log:
        "logs/combine_collected/{strategy}/{sample}_R{read}.log",
    group:
        "sample"
    conda:
        "envs/basic.yaml"
    shell:
        """
        zcat {input} | 
          gzip -c > {output} 2> {log}
        """


#### Steps after clustering ####


rule make_biom:
    input:
        tab="results/{name}/pipeline_{cluster}/{primers}/{strategy}/denoised_otutab.txt.gz",
    output:
        biom="results/{name}/pipeline_{cluster}/{primers}/{strategy}/denoised.biom",
    log:
        "logs/{name}/other/{strategy}/pipeline_{cluster}/{primers}/make_biom.log",
    group:
        "denoise"
    conda:
        "envs/biom.yaml"
    shell:
        """
        biom convert -i {input.tab} -o {output.biom} --table-type 'OTU table' --to-json
        """
