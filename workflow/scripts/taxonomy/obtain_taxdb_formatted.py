import os
import shutil
import yaml

from tax_helpers import bin_file, zstd_bin_writer, fail_on_invalid
from utils import file_logging


def obtain_formatted(param_file, outfile):
    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    with open(param_file) as f:
        par = yaml.safe_load(f)
    
    fmt = par.pop("format")
    try:
        file = par.pop("file")
    except KeyError:
        raise Exception(f"The 'file' parameter must be present with databases of type '{fmt}'")
    fail_on_invalid(par)

    with bin_file(file, gz=False) as f, zstd_bin_writer(outfile) as o:
        shutil.copyfileobj(f, o)


with file_logging(snakemake.log[0]):
    obtain_formatted(
        snakemake.input.params,
        snakemake.output.db
    )
