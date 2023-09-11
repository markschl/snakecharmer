import os
import re
from tempfile import mkstemp
import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tax_helpers import zstd_fasta_reader, zstd_fasta_writer, RANK_TRANS
from utils import file_logging

# translates rank codes to a shortened rank name (first three characters, lowercase)
rank_codes = dict(zip(*(
    RANK_TRANS.values(),
    (k[:3].lower() for k in RANK_TRANS.keys())
)))


def rank_propagate(input, output, unknown_root="Unknown", unknown_suffix="incertae_sedis"):
    """
    Does rank propagation of undefined *internal* names, while leaving the *terminal*
    undefined (empty) names empty. Thus, the output is a consistent lineage with
    all consecutive names being defined.
    If the first rank (usually kingdom/domain) is undefined, but some lower ranks
    are defined, then we fill in "Unknown".
    All lower ranks names are composed as:
    <lowest defined rank>_<rank code>_<unknown_suffix>.
    """
    prefix_re = re.compile('([a-z]+)__')
    with zstd_fasta_reader(input) as records, zstd_fasta_writer(output) as out:
        for rec in records:
            lineage = rec.description
            lineage = [r.strip() for r in lineage.split(';')]
            # first, determine the last defined rank
            n = len(lineage)
            for n in range(n - 1, -1, -1):
                is_defined = prefix_re.fullmatch(lineage[n]) is None
                if is_defined:
                    break
            # then, rank propagate all internal ranks
            # (we don't have to re-check the last defined item, so lineage[:n] is fine)
            modified = False
            last_defined_taxon = None
            last_defined_name = None
            for i, taxon in enumerate(lineage[:n]):
                m = prefix_re.fullmatch(taxon)
                if m is None:
                    # we have a name
                    last_defined_taxon = taxon
                else:
                    # empty -> construct a name
                    if last_defined_name is None and last_defined_taxon is not None:
                        # we do this splitting only once, if required
                        last_defined_name = last_defined_taxon.split("__", 1)[1]
                    rank_char = m.group(1)
                    try:
                        rank_code = rank_codes[rank_char]
                    except IndexError:
                        raise Exception(f"Unknown rank code: '{rank_char}'")
                    if last_defined_name is None:
                        name = last_defined_name = unknown_root
                    else:
                        name = f"{last_defined_name}_{rank_code}_{unknown_suffix}"
                    lineage[i] = taxon + name
                    modified = True
            if modified:
                rec.description = ";".join(lineage)
            out.write_record(rec)


### Tests ####

def _do_propagate(records):
    infile = mkstemp()[1]
    outfile = mkstemp()[1]
    with zstd_fasta_writer(infile) as w:
        for i, r in enumerate(records):
            w.write_record(SeqRecord(Seq(r[1]), id=str(i), description=r[0]))
    rank_propagate(infile, outfile)
    with zstd_fasta_reader(outfile) as recs:
        out = [(rec.description, str(rec.seq)) for rec in recs]
    os.remove(infile)
    os.remove(outfile)
    return out


class Tester(unittest.TestCase):
    # input -> output of rank propagation
    lineages = [
        (
            "k__Kingdom;f__Family;g__Genus;s__Genus species",
            "k__Kingdom;f__Family;g__Genus;s__Genus species",
        ),
        (
            "k__Kingdom;f__Family;g__;s__",
            "k__Kingdom;f__Family;g__;s__",
        ),
        (
            "k__Kingdom;f__;g__;s__",
            "k__Kingdom;f__;g__;s__",
        ),
        # we don't fill in the 'Unknown' rank if the whole lineage is empty
        (
            "k__;s__",
            "k__;s__",
        ),
        (
            "k__;p__ ;c__ ; g__Genus",
            "k__Unknown;p__Unknown_phy_incertae_sedis;c__Unknown_cla_incertae_sedis;g__Genus",
        ),
        (
            "k__Kingdom;p__; f__; g__Genus;s__Species",
            "k__Kingdom;p__Kingdom_phy_incertae_sedis;f__Kingdom_fam_incertae_sedis;g__Genus;s__Species",
        ),
    ]

    def test_propagation(self):
        modified = _do_propagate((orig, "") for orig, _ in self.lineages)
        modified = (lineage for lineage, seq in modified)
        for l, modified in zip(self.lineages, modified):
            _, expected = l
            self.assertEqual(modified, expected)


#### Entry point ####


if __name__ == '__main__':
    try:
        with file_logging(snakemake.log[0]):
            rank_propagate(snakemake.input.db, snakemake.output.db)
    except NameError:
        unittest.main()
