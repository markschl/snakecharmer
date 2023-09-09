from collections import OrderedDict
from copy import copy
import csv
import re

from utils import file_logging


class SampleReport(object):
    def __init__(self, path, pattern):
        # read data
        params = pattern.parse(path)
        self.group_header = list(params)
        self.groups = list(params.values())
        with open(path) as f:
            rdr = csv.reader(f, delimiter="\t")
            try:
                self.header = next(rdr)
            except StopIteration:
                self.header = []
            self.rows = [row for row in rdr if row]

class PathPattern(object):
    def __init__(self, pattern):
        self._pattern = re.compile(pattern)
    
    def parse(self, path):
        out = OrderedDict()
        for m in re.finditer(self._pattern, path):
            out.update(m.groupdict())
        return out


def combine_reports(sample_reports, path_pattern, outfile):
    pattern = PathPattern(path_pattern)
    sample_files = [SampleReport(path, pattern) for path in sample_reports]
    assert all(f.group_header == sample_files[0].group_header for f in sample_files[1:])
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
        w.writerow(f.group_header + header)
        empty_row = [""] * len(header)
        for f in sample_files:
            for row in f.rows:
                out = copy(empty_row)
                for i, field in zip(f.field_idx, row):
                    out[i] = field
                w.writerow(f.groups + out)


with file_logging(snakemake.log[0]):
    combine_reports(
        snakemake.input.reports,
        snakemake.params.path_pattern,
        snakemake.output.report
    )
