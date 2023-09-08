import os
from os.path import dirname

from utils import file_logging


def link_data_dir(outdir_list, clust_files, data_dir):
    # write output dirs
    outdirs = [dirname(f) for f in clust_files]
    with open(outdir_list, "w") as o:
        o.writelines((l + "\n" for l in outdirs))
    if len(outdirs) == 1:
        source_dir = outdirs[0]
        if os.path.exists(data_dir):
            os.remove(data_dir)
        os.symlink(os.path.relpath(source_dir, dirname(data_dir)), data_dir)


with file_logging(snakemake.log[0]):
    link_data_dir(outdir_list=snakemake.output.out_list,
                  clust_files=snakemake.input.clust,
                  data_dir=snakemake.params.data_dir)
