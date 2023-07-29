import sys
import os
from os.path import basename

from utils import file_logging


def link_samples(input, output):
    for source, target in zip(input, output):
        # IMPORTANT: source and target are assumed to be in the same order
        # (due to snakemake expand() using the same arguments)
        assert basename(source) == basename(target), "bug: files unordered"
        if os.path.exists(target):
            os.remove(target)
        outdir = os.path.dirname(target)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        os.symlink(os.path.abspath(source), os.path.abspath(target))
        source = os.path.relpath(source, ".")
        target = os.path.relpath(target, ".")
        print("{} > {}".format(source, target), file=sys.stdout)



with file_logging(snakemake.log[0]):
    link_samples(snakemake.input, snakemake.output)
