import os
from os.path import exists, abspath, relpath
import shutil
import sys
from typing import *

from utils import file_logging
from utils.sample_list import SampleList


def link_samples(sample_file, outdir):
    sample_list = SampleList(sample_file)
    # delete output dir to make sure there are no orphan files
    if exists(outdir):
        shutil.rmtree(outdir)
    if not exists(outdir):
        os.makedirs(outdir)
    # link
    for sample, read_paths in sample_list.samples():
        # R1 and R2 files
        for i, source in enumerate(read_paths):
            target = os.path.join(outdir, f"{sample}_R{i+1}.fastq.gz")
            os.symlink(abspath(source), abspath(target))
            # report
            source = relpath(source, ".")
            target = relpath(target, ".")
            print("{} > {}".format(source, target), file=sys.stdout)


with file_logging(snakemake.log[0]):
    link_samples(snakemake.input.samples, snakemake.output.fq)
