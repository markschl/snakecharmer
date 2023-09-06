from os.path import join

from utils import file_logging
from utils.sample_list import SampleList


def make_sample_tab(sample_file, sample_file_out, sample_dir, subdir=None, qiime_style=False):
    if subdir is not None:
        sample_dir = join(sample_dir, subdir)
    sl = SampleList(sample_file)
    with open(sample_file_out, "w") as f:
        sl.write(f, absolute_paths=True, qiime_style=qiime_style,
                 out_pattern=sample_dir+"/{sample}_R{read}.fastq.gz")


with file_logging(snakemake.log[0]):
    make_sample_tab(sample_file=snakemake.input.tab,
                    sample_file_out=snakemake.output.tab,
                    sample_dir=snakemake.input.sample_dir,
                    **snakemake.params)
