import os

from utils import file_logging


def link_data_dir(source_dirs, data_dir, list_file):
    with open(list_file, "w") as out:
        for d in source_dirs:
            print(d, file=out)
    if len(source_dirs) == 1:
        source_dir = source_dirs[0]
        if os.path.exists(data_dir):
            os.remove(data_dir)
        os.symlink(os.path.relpath(source_dir, os.path.dirname(data_dir)), data_dir)


with file_logging(snakemake.log[0]):
    link_data_dir(source_dirs=snakemake.input.outdirs,
                data_dir=snakemake.params.data_dir,
                list_file=snakemake.output.list)
