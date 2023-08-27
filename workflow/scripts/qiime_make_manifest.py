import os

from utils import file_logging
from utils.sample_list import SampleList



def make_manifest(run_config, manifest, outdir):
    # read sample files
    # TODO: demultiplexing not yet implemented
    sample_file = os.path.join(run_config, "samples.tsv")
    sl = SampleList(sample_file)
    with open(manifest, "w") as f:
        sl.write(f, qiime_style=True, absolute_paths=True, outdir=outdir)


with file_logging(snakemake.log[0]):
    make_manifest(run_config=snakemake.input.run_config,
                manifest=snakemake.output.manifest,
                  outdir=snakemake.input.fq)
