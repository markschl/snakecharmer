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
    make_primer_fasta,
    link_data_dirs,


#### Functions ####


def combine_primers(primers_by_marker):
    """
    Combines primers from different markers (in a nested dictionary structure)
    together (still grouped by forward/reverse), making sure that if the same
    primer name occurs in different markers, the primer sequences are the same
    """
    out = {"forward": {}, "reverse": {}}
    for marker, primers_by_dir in primers_by_marker.items():
        for dir_, primers in primers_by_dir.items():
            for name, seqs in primers.items():
                if name in out[dir_]:
                    assert (seqs == out[dir_][name]), \
                          "Sequences of primer {} differ between different markers".format(name)
                else:
                    out[dir_][name] = seqs
    return out


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


#### Configuration ####


rule dump_samples_yaml:
    params:
        # make sure the command is rerun if settings change
        pool = cfg.config['input']['pool_duplicates'],
        samples = cfg.samples,
        link_paths = link_paths,
    output:
        yml="results/samples.yaml",
    log:
        "logs/dump_sampes_yaml.log"
    run:
        import yaml

        # dict representation of OrderedDict in YAML
        from collections import OrderedDict
        
        class OrderedDumper(yaml.SafeDumper):
            def __init__(self, *args, **kwargs):
                super(OrderedDumper, self).__init__(*args, **kwargs)
                r = lambda self, data:  self.represent_mapping(
                    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, data.items()
                )
                self.add_representer(OrderedDict, r)


        with lib.file_logging(log):
            # YAML file
            # Even though dicts have ordered keys since Python 3.7 we
            # make sure that data is always written in the input order
            # by using OrderedDict            
            samples = OrderedDict((
                strategy, OrderedDict((
                    sample, OrderedDict((
                        "R{}".format(i + 1), [
                            relpath(path) for path in read_paths
                        ]
                        if len(read_paths) > 1
                        else relpath(read_paths[0]))
                        for i, read_paths in enumerate(zip(*sample_paths))
                    ))
                    for sample, sample_paths in strategy_paths.items()
                ))
                for strategy, strategy_paths in params.samples.items()
            )
            with open(output.yml, "w") as o:
                yaml.dump(samples, o, Dumper=OrderedDumper)



rule dump_samples_tsv:
    params:
        # make sure the command is rerun if settings change
        pool = cfg.config['input']['pool_duplicates'],
        samples = cfg.samples,
        link_paths = link_paths,
    output:
        tsv="results/samples.tsv"
    log:
        "logs/dump_sampes_tsv.log"
    run:
        with lib.file_logging(log):
            import os
            import csv
            # TSV file
            with open(output.tsv, "w") as o:
                w = csv.writer(o, delimiter="\t")
                w.writerow(["sample", "unique_sample", "strategy", "read_1_orig", "read_2_orig", "read_1", "read_2", "read_1_md5", "read_2_md5"])
                for names, paths in params.link_paths:
                    sample_name, unique_name = names
                    orig_files = [relpath(f, os.getcwd()) for f, _, _ in paths]
                    files = [f for _, _, f in paths]
                    assert len(files) <= 2
                    md5 = [lib.file_md5(p) for p, _, _ in paths]
                    if len(files) == 1:
                        w.writerow([sample_name, unique_name, 'single', orig_files[0], '', files[0], '', md5[0], ''])
                    else:
                        w.writerow([sample_name, unique_name, 'paired'] + orig_files + files + md5)


rule dump_config:
    params:
        # make sure the command is rerun if any setting changes
        pipeline_cfg = cfg.pipelines
    output:
        c="results/{name}/config.yaml",
    log:
        "logs/{name}/dump_config.log"
    run:
        import os
        import yaml
        from copy import deepcopy

        with lib.file_logging(log):
            with open(output.c, "w") as o:
                c = deepcopy(params.pipeline_cfg[wildcards.name])
                del c["settings"]["input"]
                del c["settings"]["primers"]
                del c["settings"]["taxonomy_db_sources"]
                del c["settings"]["taxonomy_dbs"]
                del c["settings"]["taxonomy_methods"]
                c["taxonomy"] = {
                    marker: {"-".join(name): config for name, config in tax.items()}
                    for marker, tax in c["taxonomy"].items()
                }
                yaml.dump(c, o)


#### Sample handling ####


rule make_primer_fasta:
    params:
        # make sure the command is rerun if settings change
        primers = cfg.primers,
        primers_rev = cfg.primers_rev
    output:
        f_primer_file="processing/primers/forward.fasta",
        r_primer_file="processing/primers/reverse.fasta",
        rprimer_rev_f="processing/primers/reverse_rev.fasta"
    log:
        "logs/make_primer_fasta.log",
    run:
        with lib.file_logging(log) as out:
            primers = combine_primers(params.primers)
            primers_rev = combine_primers(params.primers_rev)
            lib.make_primer_fasta(primers["forward"].items(), output.f_primer_file)
            lib.make_primer_fasta(primers["reverse"].items(), output.r_primer_file)
            lib.make_primer_fasta(primers_rev["reverse"].items(), output.rprimer_rev_f)


rule link_data_dirs:
    output:
        [
            directory("results/{name}/data".format(**p))
            for p in cfg.pipelines.values()
            if p["is_simple"]
        ],
    log:
        "logs/link_data_dirs.log",
    run:
        with lib.file_logging(log) as out:
            for sym_dir in output:
                p = cfg.pipelines[os.path.basename(dirname(sym_dir))]
                res_dir = join(
                    dirname(sym_dir),
                    "pipeline_" + p["cluster"],
                    p["single_primercomb"],
                    p["single_strategy"],
                )
                if not exists(res_dir):
                    os.makedirs(res_dir)
                if exists(sym_dir):
                    os.remove(sym_dir)
                os.symlink(relpath(res_dir, os.path.dirname(sym_dir)), sym_dir)


rule link_samples:
    output:
        grouped = [join("input", "grouped", path, name)
                   for _, path, name in link_paths_flat],
        unique = [join("input", "unique_samples", name)
                  for _, path, name in link_paths_flat]
    log:
        "logs/link_samples.log",
    run:
        with lib.file_logging(log) as out:
            # remove "input" dir if already present
            if exists("input"):
                shutil.rmtree("input")
            grouped = [(orig, join("input", "grouped", path, name))
                        for orig, path, name in link_paths_flat]
            unique = [(orig, join("input", "unique_samples", name))
                    for orig, path, name in link_paths_flat]
            for source, target in grouped + unique:
                if exists(target):
                    os.remove(target)
                outdir = dirname(target)
                if not exists(outdir):
                    os.makedirs(outdir)
                os.symlink(abspath(source), abspath(target))
                source = relpath(source, ".")
                target = relpath(target, ".")
                print("{} > {}".format(source, target), file=out)


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


#### QC ####


rule fastqc:
    input:
        "input/grouped/{strategy}/{sample}/{prefix}.fastq.gz",
    output:
        html="results/_validation/fastqc/{strategy,[^/]+}/{sample,[^/]+}/{prefix}_fastqc.html",
        zip="results/_validation/fastqc/{strategy,[^/]+}/{sample,[^/]+}/{prefix}_fastqc.zip",
    log:
        "logs/fastqc/{strategy}/{sample}/{prefix}.log",
    group:
        "qc"
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
        [join("results", "_validation", "fastqc", splitext(splitext(join(path, name))[0])[0] + "_fastqc.html")
        for _, path, name in link_paths_flat],
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
    group:
        "denoise"
    conda:
        "envs/biom.yaml"
    shell:
        """
        biom convert -i {input.tab} -o {output.biom} --table-type 'OTU table' --to-json
        """
