import csv
from os.path import abspath

from utils import file_logging


def manifest_paired(output, sample_names):
    with open(output, "w") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(
            ["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"]
        )
        for sample in sample_names:
            paths = [
                abspath(
                    "input/grouped/paired/{sample}/{sample}_R{read}.fastq.gz".format(
                        sample=sample, read=read
                    )
                )
                for read in [1, 2]
            ]
            w.writerow([sample] + paths)


with file_logging(snakemake.log[0]):
    manifest_paired(snakemake.output[0], snakemake.params.samples)
