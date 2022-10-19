import logging
import gzip
import shutil
from collections import OrderedDict
from sys import stderr
from os.path import splitext, dirname, join
import os
import lib

localrules:
    dump_samples,
    dump_config,
    setup_project,


#### Functions ####


def get_link_paths(samples, outdir):
    link_paths = []
    combine_samples = {}
    for sample_name, by_read in samples:
        sample_dir = join(outdir, sample_name)
        for read_num, files in by_read:
            if len(files) == 1:
                root, f = next(iter(files))
                f = join(root, f)
                out = join(sample_dir, "{}_R{}.fastq.gz".format(sample_name, read_num))
                link_paths.append((f, out))
            else:
                files = sorted(files)
                files = [join(root, f) for root, f in files]
                file_i = list(range(1, len(files) + 1))
                names = [
                    "{}_s{}_R{}.fastq.gz".format(sample_name, i, read_num)
                    for f, i in zip(files, file_i)
                ]
                link_paths.extend(
                    (f, join(sample_dir, "_multiple", name))
                    for f, name in zip(files, names)
                )
                assert not (sample_name, read_num) in combine_samples
                combine_samples[(sample_name, read_num)] = file_i
    return link_paths, combine_samples


def combine_primers(primers_by_marker):
    out = {'forward': {}, 'reverse': {}}
    for marker, primers_by_dir in primers_by_marker.items():
        for dir_, primers in primers_by_dir.items():
            for name, seqs in primers.items():
                if name in out[dir_]:
                    assert seqs == out[dir_][name], 'Sequences of primer {} differ between different markers'.format(name)
                else:
                    out[dir_][name] = seqs
    return out


#### Start ####

cfg = lib.Config(config)


# determine paths to link
sample_data = {
    strategy: get_link_paths(cfg.samples[strategy], strategy)
    for strategy in ("single", "paired")
}
link_paths = [p for paths, _ in sample_data.values() for p in paths]
combine_samples = {strategy: data[1] for strategy, data in sample_data.items()}


#### Configuration ####


rule dump_samples:
    output:
        samples="results/samples.yaml",
    run:
        import os
        import yaml

        samples = {
            strategy: {
                sample: {
                    "R{}".format(r): [
                        os.path.relpath(join(prefix, fname)) for prefix, fname in files
                    ]
                    if len(files) > 1
                    else os.path.relpath(join(*list(files)[0]))
                    for r, files in sdata
                }
                for sample, sdata in s
            }
            for strategy, s in cfg.samples.items()
        }
        with open(output.samples, "w") as o:
            yaml.dump(samples, o)


rule dump_config:
    output:
        c="results/{name}/config.yaml",
    run:
        import os
        import yaml
        from copy import deepcopy

        with open(output.c, "w") as o:
            c = deepcopy(cfg[wildcards.name])
            del c["settings"]["input"]
            del c["settings"]["primers"]
            del c["settings"]["taxonomy_db_sources"]
            del c["settings"]["taxonomy_dbs"]
            del c["settings"]["taxonomy_methods"]
            c["taxonomy"] = {
                marker: {
                    "-".join(name): config for name, config in tax.items()
                }
                for marker, tax in c["taxonomy"].items()
            }
            yaml.dump(c, o)


#### Sample handling ####


rule setup_project:
    output:
        fprimers="processing/primers/forward.fasta",
        rprimers="processing/primers/reverse.fasta",
        rprimers_rev="processing/primers/reverse_rev.fasta",
        simple_sym=[
            directory("results/{name}/data".format(**p))
            for p in cfg.pipelines.values()
            if p["is_simple"]
        ],
    log:
        "logs/setup_project.log",
    run:
        with lib.file_logging(log) as out:
            # primer FASTA
            primers = combine_primers(cfg.primers)
            primers_rev = combine_primers(cfg.primers_rev)
            lib.make_primer_fasta(primers['forward'], output.fprimers)
            lib.make_primer_fasta(primers['reverse'], output.rprimers)
            lib.make_primer_fasta(primers_rev['reverse'], output.rprimers_rev)
            # symlink dirs
            for sym_dir in output.simple_sym:
                p = cfg.pipelines[os.path.basename(dirname(sym_dir))]
                res_dir = join(
                    dirname(sym_dir),
                    "pipeline_" + p["cluster"],
                    p["single_primercomb"],
                    p["single_strategy"],
                )
                if not os.path.exists(res_dir):
                    os.makedirs(res_dir)
                if os.path.exists(sym_dir):
                    os.remove(sym_dir)
                os.symlink(os.path.relpath(res_dir, os.path.dirname(sym_dir)), sym_dir)


rule link_samples:
    output:
        [join("input", p) for _, p in link_paths],
    log:
        "logs/link_samples.log",
    run:
        with lib.file_logging(log) as out:
            for source, target in link_paths:
                target = join("input", target)
                if os.path.exists(target):
                    os.remove(target)
                outdir = os.path.dirname(target)
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                os.symlink(os.path.abspath(source), os.path.abspath(target))
                source = os.path.relpath(source, ".")
                target = os.path.relpath(target, ".")
                print("{} > {}".format(source, target), file=out)


rule combine_multiple_samples:
    input:
        lambda w: expand(
            "input/{{strategy}}/{{sample}}/_multiple/{{sample}}_s{file_i}_R{{read}}.fastq.gz",
            file_i=combine_samples[w.strategy][(w.sample, int(w.read))],
        ),
    output:
        "input/{strategy}/{sample}/{sample}_R{read}.fastq.gz",
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


#### QC ####


rule fastqc:
    input:
        "input/{strategy}/{sample}/{prefix}.fastq.gz",
    output:
        html="results/_validation/fastqc/{strategy,[^/]+}/{sample,[^/]+}/{prefix}_fastqc.html",
        zip="results/_validation/fastqc/{strategy,[^/]+}/{sample,[^/]+}/{prefix}_fastqc.zip",
    log:
        "logs/fastqc/{strategy}/{sample}/{prefix}.log",
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        outdir=$(dirname {output.html})
        mkdir -p $outdir
        fastqc -q -f fastq -t 1 -o "$outdir" {input} 2> {log}
        """


rule multiqc_fastqc:
    input:
        [
            join(
                "results",
                "_validation",
                "fastqc",
                splitext(splitext(p)[0])[0] + "_fastqc.html",
            )
            for _, p in link_paths
        ],
    output:
        "results/_validation/multiqc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    conda:
        "envs/multiqc.yaml"
    shell:
        "multiqc -fm fastqc -o $(dirname {output}) results/_validation/fastqc 2> {log}"


#### Steps after clustering ####


rule make_biom:
    input:
        tab="results/{name}/pipeline_{cluster}/{primers}/{strategy}/denoised_otutab.txt.gz",
    output:
        biom="results/{name}/pipeline_{cluster}/{primers}/{strategy}/denoised.biom",
    log:
        "logs/{name}/other/{strategy}/pipeline_{cluster}/{primers}/make_biom.log",
    conda:
        "envs/biom.yaml"
    shell:
        """
        biom convert -i {input.tab} -o {output.biom} --table-type 'OTU table' --to-json
        """
