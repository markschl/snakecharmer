import os
import csv
import hashlib

from utils import file_logging


def file_md5(filename):
    md5 = hashlib.md5()
    with open(filename, "rb") as f:
        # Read and update hash in chunks of 4K
        for chunk in iter(lambda: f.read(4096), b""):
            md5.update(chunk)
        return md5.hexdigest()


def dump_tsv(link_paths, outfile):
    with open(outfile, "w") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["sample", "unique_sample", "strategy", "read_1_orig", "read_2_orig", "read_1", "read_2", "read_1_md5", "read_2_md5"])
        for names, paths in link_paths:
            sample_name, unique_name = names
            orig_files = [os.path.relpath(f, os.getcwd()) for f, _, _ in paths]
            files = [f for _, _, f in paths]
            assert len(files) <= 2
            md5 = [file_md5(p) for p, _, _ in paths]
            if len(files) == 1:
                w.writerow([sample_name, unique_name, 'single', orig_files[0], '', files[0], '', md5[0], ''])
            else:
                w.writerow([sample_name, unique_name, 'paired'] + orig_files + files + md5)


with file_logging(snakemake.log[0]):
    dump_tsv(snakemake.params.link_paths, snakemake.output.tsv)
