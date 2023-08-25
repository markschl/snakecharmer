import os
from os.path import abspath
from multiprocessing.dummy import Pool   # only threaded version necessary
from subprocess import check_call
import sys
# import gzip
# import shutil

import yaml

from utils import file_logging


def do_pooling(data):
    source_paths, target_path = data
    # make sure that all file names are the same
    assert len(set(os.path.basename(p) for p in list(source_paths) + [target_path])) == 1
    try:
        if len(source_paths[0]) == 1:
            # we only need to link
            assert len(source_paths) == 1
            os.symlink(abspath(source_paths[0]), abspath(target_path))
        else:
            d = os.path.dirname(target_path)
            if not os.path.exists(d):
                os.makedirs(d)
            check_call(["zcat {} | gzip -c > {}".format(" ".join(source_paths), target_path)], shell=True)
            # with gzip.open(target_path, 'w') as out:
            #     for path in source_paths:
            #         print("source", path, "to", target_path)
            #         with gzip.open(path) as f:
            #             shutil.copyfileobj(f, out)
            return data
    except Exception as e:
        return e


def pool_raw(pooling_info, outdir, processes=1):
    # list of input -> output paths
    with open(pooling_info) as f:
        pool_files = yaml.safe_load(f)
    # send these lists to the workers
    p = Pool(processes)
    # Create a flat sequence of read files to combine
    # and the corresponding output file name
    args = ((f, os.path.join(outdir, f"{sample}_{read}.fastq.gz"))
            for sample, sample_files in pool_files.items()
            for read, f in sample_files.items())
    for res in p.imap_unordered(do_pooling, args):
        if isinstance(res, Exception):
            raise res
        source_paths, target_path = res
        print("{} > {}".format(" ".join(source_paths), target_path),
              file=sys.stderr)


with file_logging(snakemake.log[0]):
    pool_raw(pooling_info=snakemake.input.yml,
                outdir=snakemake.output.fq,
                processes=snakemake.threads)
