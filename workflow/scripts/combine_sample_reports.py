from copy import copy
import csv

from utils import file_logging


class SampleFile(object):
    def __init__(self, path):
        # parse path
        p = path.split("/")
        self.run, self.layout = p[-2].split("_", 1)
        # read data
        with open(path) as f:
            rdr = csv.reader(f, delimiter="\t")
            try:
                self.header = next(rdr)
            except StopIteration:
                self.header = []
            self.rows = [row for row in rdr if row]


def combine_reports(sample_reports, outfile):
    sample_files = [SampleFile(path) for path in sample_reports]
    # obtain header
    header = []
    for f in sample_files:
        f.field_idx = []
        for field in f.header:
            try:
                f.field_idx.append(header.index(field))
            except ValueError:
                f.field_idx.append(len(header))
                header.append(field)
    # write data
    with open(outfile, "w") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["run", "layout"] + header)
        empty_row = [""] * (2 + len(header))
        for f in sample_files:
            for row in f.rows:
                out = copy(empty_row)
                out[0] = f.run
                out[1] = f.layout
                for i, field in zip(f.field_idx, row):
                    out[2 + i] = field
                w.writerow(out)


with file_logging(snakemake.log[0]):
    combine_reports(snakemake.input.reports, snakemake.output.report)
