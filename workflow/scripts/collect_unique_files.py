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


def get_samples(run_cfg_dirs):
    # read sample files
    # TODO: indexes.tsv for demultiplexing
    sample_files = [os.path.join(d, "samples.tsv") for d in run_cfg_dirs]
    for f in sample_files:
        if not os.path.exists(f):
            shutil.rmtree(os.path.dirname(f))
            raise Exception(f"Sample file {sample_files[0]} does not exist")
    sample_lists = [SampleList(f) for f in sample_files]
    # get list of runs per sample
    sample2run = defaultdict(set)
    for path, l in zip(sample_files, sample_lists):
        # annotate run and technology
        p = path.split(os.sep)
        l.run = p[-2]
        l.technology = p[-4]
        assert l.layout == p[-3]  # layout should be correct already
        for sample, _ in l.samples():
            sample2run[sample].add((l.layout, l.run))
    # Which layout/run combinations have duplicate names?
    # Runs/layouts without dupes are set to None
    dupes = set(comb
                for r in sample2run.values() if len(r) > 1
                for comb in r)
    dupes = {r: r for r in dupes}  # convert to dict, where only values are further modified
    if len(set(layout for layout, _ in dupes.values())) == 1:
        dupes = {k: (None, v[1]) for k, v in dupes.items()}
    if len(set(run for _, run in dupes.values())) == 1:
        dupes = {k: (v[0], None) for k, v in dupes.items()}
    combinations = sorted(set(dupes.values()))
    assert not (None, None) in combinations
    # assign unique suffixes
    suffixes = {k: combinations.index(dupes[k]) + 1 for k in dupes}
    sample_dict = {}
    for l in sample_lists:
        try:
            suffix = suffixes[(l.layout, l.run)]
        except KeyError:
            suffix = None
        meta = (l.technology, l.layout, l.run)
        for sample, reads in l.samples():
            unique_sample = sample if suffix is None else f"{sample}_{suffix}"
            assert not unique_sample in sample_dict
            paths = [(path, f"{unique_sample}_R{i+1}.fastq.gz")
                     for i, path in enumerate(reads)]
            sample_dict[unique_sample] = (sample, meta, paths)
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


def link_samples(run_cfg_dirs, read_dir, sample_file):
    # delete output dir to make sure there are no orphan files
    if exists(read_dir):
        shutil.rmtree(read_dir)
    if not exists(read_dir):
        os.makedirs(read_dir)
    # get the sample paths
    samples = get_samples(run_cfg_dirs)
    # do the work
    do_link((paths for _, _, _, paths in samples), read_dir)
    dump_tsv(samples, sample_file)


with file_logging(snakemake.log[0]):
    link_samples(snakemake.input.run_cfg_dirs,  snakemake.output.read_dir, snakemake.output.tsv)
