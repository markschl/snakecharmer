from collections import defaultdict
import csv
import hashlib
import os
from os.path import exists, abspath, relpath
import shutil
import sys
from typing import *

from utils import file_logging
from utils.sample_list import SampleList


def get_samples(run_meta, path_template):
    # get all necessary metadata
    _run_meta = [
        (
            (d["technology"], d["layout"], d["run"]),
            list(SampleList(path_template.format(**d).samples()))
        )
        for d in run_meta
    ]
    # get list of runs per sample
    sample2run = defaultdict(set)
    for run, samples in _run_meta:
        for sample in samples:
            sample2run[sample].add(run)
    # Which layout/run/technology combinations have duplicate samples?
    # Runs/layouts without duplicates are set to None
    dupes = sorted(set(
        run
        for runs in sample2run.values() if len(runs) > 1
        for run in runs
    ))
    # Define unique suffixes for run/layout combination that have duplicate samples.
    # For consistency, all sample names in these runs will receive a suffix, 
    # even if some of them are not duplicated.
    suffixes = {r: dupes.index(dupes[r]) + 1 for r in dupes}
    sample_dict = {}
    for run, samples in _run_meta:
        try:
            suffix = suffixes[run]
        except KeyError:
            suffix = None
        for sample, reads in samples:
            unique_sample = sample if suffix is None else f"{sample}_{suffix}"
            assert not unique_sample in sample_dict
            paths = [(path, f"{unique_sample}_R{i+1}.fastq.gz")
                     for i, path in enumerate(reads)]
            sample_dict[unique_sample] = (sample, run, paths)
    # convert to flat sorted list
    samples = [(unique_sample, *other) for unique_sample, other in sample_dict.items()]
    return sorted(samples)


def file_md5(filename):
    md5 = hashlib.md5()
    with open(filename, "rb") as f:
        # Read and update hash in chunks of 4K
        for chunk in iter(lambda: f.read(4096), b""):
            md5.update(chunk)
        return md5.hexdigest()


def do_link(path_iter, outdir):
    for paths in path_iter:
        for source, target in paths:
            target = os.path.join(outdir, target)
            assert not os.path.exists(target)
            os.symlink(abspath(source), abspath(target))
            # report
            source = relpath(source, ".")
            target = relpath(target, ".")
            print("{} > {}".format(source, target), file=sys.stdout)


def dump_tsv(samples, outfile):
    with open(outfile, "w") as o:
        w = csv.writer(o, delimiter="\t")
        header = ["technology", "run", "layout", 
                  "sample", "unique_sample",
                  "source_read_1", "source_read_2",
                  "read_1", "read_2",
                  "read_1_md5", "read_2_md5"]
        w.writerow(header)
        for unique_sample, sample, meta, paths in samples:
            orig_reads, target_files = zip(*paths)
            orig_reads = list(orig_reads)
            target_files = list(target_files)
            md5 = [file_md5(p) for p in orig_reads]
            if len(orig_reads) == 1:
                orig_reads.append("")
                target_files.append("")
                md5.append("")
            row = list(meta) + [sample, unique_sample] + orig_reads + target_files + md5
            w.writerow(row)


def link_samples(run_meta, read_dir, sample_file):
    # delete output dir to make sure there are no orphan files
    if exists(read_dir):
        shutil.rmtree(read_dir)
    if not exists(read_dir):
        os.makedirs(read_dir)
    # get the sample paths
    samples = get_samples(run_meta)
    # do the work
    do_link((paths for _, _, _, paths in samples), read_dir)
    dump_tsv(samples, sample_file)


with file_logging(snakemake.log[0]):
    link_samples(snakemake.params.run_meta, snakemake.input.path_template, 
                 snakemake.output.read_dir, snakemake.output.tsv)
