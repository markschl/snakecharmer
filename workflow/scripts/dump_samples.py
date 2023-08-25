import csv
import hashlib
from pathlib import Path

from utils import file_logging
from utils.sample_list import SampleList


def list_samples(sample_dirs, outfile):
    with open(outfile, "w") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["technology", "run", "layout", "sample", "read_1", "read_2"])
        for sample_dir in sample_dirs:
            if os.path.exists("")
            sl = SampleList(sample_file)
            # we can infer some values from the path
            path = Path(sample_file).parts
            technology = path[-4]
            layout = path[-3]
            run_name = path[-2]
            assert layout == sl.layout
            for sample, read_paths in sl.samples():
                md5 = [file_md5(p) for p in read_paths]
                if len(read_paths) == 1:
                    read_paths.append("")
                    md5.append("")
                w.writerow([technology, run_name, layout, sample] + read_paths + md5)


with file_logging(snakemake.log[0]):
    list_samples(snakemake.input.sample_files, snakemake.output.tsv)
