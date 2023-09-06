import os

from utils import file_logging


def link_data_dir(clust_files, data_dir):
    if len(clust_files) == 1:
        source_dir = os.path.dirname(clust_files[0])
        if os.path.exists(data_dir):
            os.remove(data_dir)
        os.symlink(os.path.relpath(source_dir, os.path.dirname(data_dir)), data_dir)


with file_logging(snakemake.log[0]):
    link_data_dir(clust_files=snakemake.input.clust,
                  data_dir=snakemake.params.data_dir)
