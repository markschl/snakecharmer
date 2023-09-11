"""
QIIME-style to UTAX-style database conversion.
Reference: https://www.drive5.com/usearch/manual/tax_annot.html
"""

import re

import yaml

from tax_helpers import zstd_fasta_reader, zstd_fasta_writer, fail_on_invalid
from utils import file_logging



def convert_taxdb_utax(input, param_file, output):
    # we don't allow any parameters
    with open(param_file) as f:
        params = yaml.safe_load(f)
        fail_on_invalid(params)

    reserved_chars = re.compile("[ ,:]")
    rank_pat = re.compile('\s*?([a-z]+)__(.*?)\s*')
    with zstd_fasta_reader(input) as records, zstd_fasta_writer(output) as out:
        for rec in records:
            lineage = rec.description
            # some characters have a special meaning in the UTAX format, convert to '_'
            # (including spaces in names)
            lineage = reserved_chars.sub('_', lineage)
            # split into components
            try:
                lineage = [rank_pat.fullmatch(r).groups() for r in lineage.split(';')]
            except AttributeError:
                raise Exception("Not a valid QIIME-formatted lineage: {}".format(lineage))
            # remove empty ranks, since the UTAX format does not require every rank
            # in every lineage
            lineage_out = [(rank, name) for rank, name in lineage if name.strip() != ""]
            # if lineage != lineage_out:
            #     print("before", lineage)
            #     print("filter", lineage_out)
            # final format
            rec.id = "{id};tax={tax};".format(
                id=rec.id,
                tax=",".join(rank + ':' + name for rank, name in lineage_out)
            )
            rec.description = ""
            out.write_record(rec)


with file_logging(snakemake.log[0]):
    convert_taxdb_utax(snakemake.input.db, snakemake.input.params, snakemake.output.db)
