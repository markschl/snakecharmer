import sys
import os
from os.path import exists, dirname, join, relpath, abspath
import shutil

from utils import file_logging


def link_samples(link_paths):
    # remove "input" dir if already present
    if exists("input"):
        shutil.rmtree("input")
    grouped = [(orig, join("input", "grouped", path, name))
                for orig, path, name in link_paths]
    unique = [(orig, join("input", "unique_samples", name))
            for orig, _path, name in link_paths]
    for source, target in grouped + unique:
        if exists(target):
            os.remove(target)
        outdir = dirname(target)
        if not exists(outdir):
            os.makedirs(outdir)
        os.symlink(abspath(source), abspath(target))
        source = relpath(source, ".")
        target = relpath(target, ".")
        print("{} > {}".format(source, target), file=sys.stdout)


with file_logging(snakemake.log[0]):
    link_samples(snakemake.params.link_paths)
