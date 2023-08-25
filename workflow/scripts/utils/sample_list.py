import csv
from os.path import abspath, join
import re
import sys

import yaml


class SampleList(object):
    default_header = {
        "paired": ["id", "R1", "R2"],
        "single": ["id", "R1"]
    }
    qiime_header = {
        "paired": ["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"],
        "single": ["sample-id", "absolute-filepath"]
    }

    def __init__(self, sample_file=None, layout=None, reserved_chars=None):
        self._samples = []
        if reserved_chars is None:
            self.reserved_re = None
        else:
            self.reserved_re = re.compile(f"[{reserved_chars}]")
        if sample_file is not None:
            assert layout is None
            self._read_samples(sample_file)
        else:
            assert layout in ("single", "paired")
            self.layout = layout
            self.ncol = 3 if layout == "paired" else 2
            self.header = None
        self.n_reads = self.ncol - 1

    def _read_samples(self, sample_file):
        with open(sample_file) as f:
            rdr = csv.reader(f, delimiter="\t")
            self.header = next(rdr)
            self.ncol = len(self.header)
            if self.ncol == 3:
                assert self.header == self.default_header["paired"] or self.header == self.qiime_header["paired"], (
                    "Unknown paired-end sample file header: {}. "
                    "Valid are {} or {}").format(self.header, self.header["single"], self.qiime_header["single"])
                self.layout = "paired"
            else:
                assert self.ncol == 2, (
                    "Invalid number of columns in sample file. "
                    "Either two (single-end) or three (paired-end) are expected")
                assert self.header == self.header["single"] or self.header == self.qiime_header["single"], (
                    "Unknown paired-end sample file header: {}. "
                    "Valid are {} or {}").format(self.header, self.header["single"], self.qiime_header["single"])
                self.layout = "single"

            for row in rdr:
                self.add(row[0], row[1:])
    
    def add(self, sample, reads):
        row = [sample] + list(reads)
        assert len(row) == self.ncol
        if self.reserved_re is not None:
            (row[0], n) = self.reserved_re.subn("_", row[0])
            if n > 0:
                print(f"Reserved characters replaced in sample name: {row[0]}",
                      file=sys.stderr)
        self._samples.append(row)

    def samples(self):
        for row in self._samples:
            yield row[0], row[1:]

    def write(self, handle, qiime_style=False, outdir=None, absolute_paths=False):
        wtr = csv.writer(handle, delimiter="\t")
        header = self.qiime_header if qiime_style else self.default_header
        wtr.writerow(header[self.layout])
        for row in self._samples:
            if outdir is not None:
                row[1:] = [join(outdir, f) for f in row[1:]]
            if absolute_paths:
                row[1:] = [abspath(f) for f in row[1:]]
            wtr.writerow(row)

    def write_yaml(self, handle):
        out = {}
        for sample, reads in self.samples():
            if len(reads) == 1:
                out[sample] = reads[0]
            else:
                out[sample] = {f"R{i+1}": f for i, f in enumerate(reads)}
        yaml.safe_dump(out, handle)
