import os
import yaml
from copy import deepcopy

from utils import file_logging


def dump_config(config, outfile):
    with open(outfile, "w") as o:
        c = deepcopy(config)  # make copy since we modify it
        del c["settings"]["input"]
        del c["settings"]["primers"]
        del c["settings"]["taxonomy_db_sources"]
        del c["settings"]["taxonomy_dbs"]
        del c["settings"]["taxonomy_methods"]
        c["taxonomy"] = {
            marker: {"-".join(name): config for name, config in tax.items()}
            for marker, tax in c["taxonomy"].items()
        }
        yaml.dump(c, o)


with file_logging(snakemake.log[0]):
    dump_config(snakemake.params.config, snakemake.output[0])
