import os
import re
import sys
from tempfile import mkstemp
import unittest

import yaml
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tax_helpers import *
from utils import file_logging


class Filters(object):
    """
    Class handling the filters
    """
    # QIIME-style rank prefixes
    _prefix_re = re.compile('([a-z]+)__')

    # for ambiguous filtering
    DNA = set(b"ACGT")
    Ambigs = set(b"MRWSYKVHDBN")
    
    # list of filter keywords
    keywords = ['defined_rank', 'min_len', 'max_len', 'max_n', 'max_ambig']

    def __init__(
            self,
            defined_rank = None,
            min_len = None,
            max_len = None,
            max_n = None,
            max_ambig = None
    ):
        self.defined_rank = defined_rank
        self.min_len = min_len
        self.max_len = max_len
        self.max_n = max_n
        self.max_ambig = max_ambig
        self.initialized = False

    def init(self, lineage):
        """
        Set up filters given the first record
        """
        # any ambiguity filter present?
        self._ambig_filter = self.max_n is not None or self.max_ambig is not None
        self._no_ambigs = self._ambig_filter and self.max_n in (None, 0) and self.max_ambig in (None, 0)

        # validate defined_rank
        self.defined_rank_prefix = None
        if self.defined_rank is not None:
            try:
                ranks_short = [self._prefix_re.match(r.strip()).group(1) for r in lineage.split(';')]
            except AttributeError:
                raise Exception("Not a valid QIIME-formatted lineage: {}".format(lineage))
            assert self.defined_rank in RANK_TRANS, "Unknown rank name for filtering: {}".format(self.defined_rank)
            rank_short = RANK_TRANS[self.defined_rank]
            self.defined_rank_prefix = rank_short + '__'
            if rank_short not in ranks_short:
                rev_dict = dict(zip(RANK_TRANS.values(), RANK_TRANS.keys()))
                raise Exception((
                    "Minimum unambiguous rank name not found in lineage: {}. "
                    "Even though this is a valid rank, it should actually exist in the database."
                    "Available ranks: {}").format(
                        self.defined_rank, ", ".join(rev_dict[r] for r in ranks_short)))
            # ensure ranks are in order
            possible_ranks = list(RANK_TRANS.values())
            try:
                rank_i = [possible_ranks.index(r) for r in ranks_short]
            except IndexError:
                raise Exception("Invalid rank found among the ranks of the first lineage: {}".format(", ".join(ranks_short)))
            assert rank_i == sorted(rank_i), "Ranks in the first lineage are not in the correct order"


    def check(self, id: str, lineage: str, seq: Seq) -> bool:
        # initialize if necessary
        if not self.initialized:
            self.init(lineage)
            self.initialized = True

        # ambiguous rank filter
        if self.defined_rank_prefix is not None:
            lineage = [l.strip() for l in lineage.split(";")]
            try:
                # if the empty prefix is found, this means that the rank
                # is not defined.
                undef_i = lineage.index(self.defined_rank_prefix)
                # This hit could be an intermediate "blank" rank, therefore we have to further
                # check if all lower ranks are empty. This procedure should still be faster
                # than iterating backwards on every lineage, intermediate blank
                # ranks are usually not so frequent.
                if all(self._prefix_re.fullmatch(r) is not None for r in lineage[(undef_i+1):]):
                    return False
            except ValueError:
                pass
        
        # length filters
        if self.min_len is not None and len(seq) < self.min_len:
            return False
        if self.max_len is not None and len(seq) > self.max_len:
            return False
        
        # sequence ambiguity filters
        if self._ambig_filter:
            # first, check if there are any ambiguities at all
            seq_bytes = bytes(seq)
            bases = set(seq_bytes)
            ambigs = bases.difference(self.DNA)
            # validate
            invalid = ambigs.difference(self.Ambigs)
            assert len(invalid) == 0,\
                "Invalid characters found in sequence: {}".format(
                    ", ".join(chr(b) for b in sorted(invalid))
                )
            if len(ambigs) > 0:
                # found some ambigs: in case we don't allow any, we can return
                if self._no_ambigs:
                    return False
                # otherwise, we count DNA bases to obtain the number of ambigs
                # This solution seems to be fastest (faster than using collections.Counter)
                dna_count = sum(seq_bytes.count(b) for b in self.DNA)
                ambig_count = len(seq_bytes) - dna_count
                if self.max_ambig is not None and ambig_count > self.max_ambig:
                    return False
                if self.max_n is not None and seq_bytes.count(b"N") > self.max_n:
                    return False
        # if none of the checks failed, the sequence may be included
        return True                


def filter_taxdb(input, filtered_out, cfg_file):

    with open(cfg_file) as f:
        cfg = yaml.safe_load(f)
    
    invalid = [k for k in cfg if not k in Filters.keywords]
    assert len(invalid) == 0, \
        "Invalid taxonomy database filter keyword(s) found: " + " ".join(invalid)

    if len(cfg) > 0:
        # there seems to be something to filter
        with zstd_fasta_reader(input) as records, zstd_fasta_writer(filtered_out) as out:
            filters = Filters(**cfg)
            written = 0
            i = -1
            for i, rec in enumerate(records):
                if filters.check(rec.id, rec.description, rec.seq) is True:
                    out.write_record(rec)
                    written += 1
            assert i is not None, "No records found"
            print(f"{written} of {i+1} records written to output", file=sys.stderr)
    else:
        # nothing to do, simply link to the output file
        if os.path.exists(filtered_out):
            os.remove(filtered_out)
        print(f"No filtering to be done, linking {input} to {filtered_out}", file=sys.stderr)
        os.symlink(os.path.relpath(input, os.path.dirname(filtered_out)), filtered_out)



### Tests ####

def _do_filter(records, **cfg):
    infile = mkstemp()[1]
    outfile = mkstemp()[1]
    cfg_file = mkstemp()[1]
    with open(cfg_file, "w") as f:
        yaml.safe_dump(cfg, f)
    with zstd_fasta_writer(infile) as w:
        for i, r in enumerate(records):
            w.write_record(SeqRecord(Seq(r[1]), id=str(i), description=r[0]))
    filter_taxdb(infile, outfile, cfg_file)
    with zstd_fasta_reader(outfile) as recs:
        out = [(rec.description, str(rec.seq)) for rec in recs]
    os.remove(infile)
    os.remove(outfile)
    os.remove(cfg_file)
    return out


class Tester(unittest.TestCase):
    tax_records = [
        ("k__Kingdom;f__Family;g__Genus;s__Genus species", "TGTATAGCAATAG"),
        ("k__Kingdom;f__Family;g__;s__", "TGTATAGCAATAG"),
        ("k__Kingdom;f__;g__;s__", "TGTATAGCAATAG"),
        # lineages with intermediate blanks
        ("k__;f__Family;g__;s__", "A"),
        ("k__;f__;g__Genus;s__Genus_species", "A"),
    ]
    tax_with_blanks = [
    ]
    invalid_tax = [
        ("f__Family;k__;g__;s__", ""),
        ("", ""),
        ("k_", ""),
        ("k__;s__;", ""),
    ]
    ambig_records = [
        ("k__Kingdom;s__Species", "TGTATAGCAATAG"),
        ("k__Kingdom;s__Species", "NNAGCAA"),  # length = 7
        ("k__Kingdom;s__Species", "BRYAGCAATAG"),
        ("k__Kingdom;s__Species", "RNNAGCAAT"),  # length = 9
    ]
    invalid_seq = [
        ("k__", "ATGZX")
    ]

    def test_no_filter(self):
        rec = self.tax_records + self.invalid_tax + self.ambig_records + self.invalid_seq
        flt = _do_filter(rec)
        self.assertEqual(flt, rec)

    def test_rank_filter(self):
        flt = _do_filter(self.tax_records, defined_rank="species")
        self.assertEqual(flt, [self.tax_records[i] for i in [0, 4]])
        flt = _do_filter(self.tax_records, defined_rank="genus")
        self.assertEqual(flt, [self.tax_records[i] for i in [0, 4]])
        flt = _do_filter(self.tax_records, defined_rank="family")
        self.assertEqual(flt, [self.tax_records[i] for i in [0, 1, 3, 4]])
        with self.assertRaises(Exception) as ctx:
            _do_filter(self.tax_records, defined_rank="subfamily")
        self.assertTrue('Minimum unambiguous rank name not found in lineage: subfamily.' in str(ctx.exception))
        with self.assertRaises(Exception) as ctx:
            _do_filter(self.tax_records, defined_rank="noname")
        self.assertTrue('Unknown rank name for filtering: noname' in str(ctx.exception))

    def test_ambig(self):
        flt = _do_filter(self.ambig_records, max_ambig=3)
        self.assertEqual(flt, [self.ambig_records[i] for i in [0, 1, 2, 3]])
        flt = _do_filter(self.ambig_records, max_ambig=2)
        self.assertEqual(flt, [self.ambig_records[i] for i in [0, 1]])
        flt = _do_filter(self.ambig_records, max_ambig=2, max_n=1)
        self.assertEqual(flt, [self.ambig_records[i] for i in [0]])
        flt = _do_filter(self.ambig_records, max_n=2)
        self.assertEqual(flt, [self.ambig_records[i] for i in [0, 1, 2, 3]])
        flt = _do_filter(self.ambig_records, max_n=1)
        self.assertEqual(flt, [self.ambig_records[i] for i in [0, 2]])
        flt = _do_filter(self.ambig_records, max_n=0)
        self.assertEqual(flt, [self.ambig_records[i] for i in [0]])

    def test_lenfilter(self):
        flt = _do_filter(self.ambig_records, min_len=7)
        self.assertEqual(flt, [self.ambig_records[i] for i in [0, 1, 2, 3]])
        flt = _do_filter(self.ambig_records, min_len=8)
        self.assertEqual(flt, [self.ambig_records[i] for i in [0, 2, 3]])
        flt = _do_filter(self.ambig_records, min_len=10)
        self.assertEqual(flt, [self.ambig_records[i] for i in [0, 2]])
        flt = _do_filter(self.ambig_records, min_len=8, max_len=9)
        self.assertEqual(flt, [self.ambig_records[i] for i in [3]])

    def test_comb(self):
        rec = self.tax_records + self.ambig_records
        flt = _do_filter(rec, defined_rank="genus", max_ambig=3, max_n=0)
        self.assertEqual(flt, [rec[i] for i in [0, 4, 5, 7]])
        flt = _do_filter(rec, defined_rank="kingdom", max_ambig=2, max_n=1)
        self.assertEqual(flt, [rec[i] for i in [0, 1, 2, 3, 4, 5]])
        flt = _do_filter(rec, defined_rank="species", max_ambig=2, min_len=8)
        self.assertEqual(flt, [rec[i] for i in [0, 5]])

    def test_invalid_tax(self):
        # errors are only checked with rank filter
        with self.assertRaises(Exception) as ctx:
            _do_filter([self.invalid_tax[0]], defined_rank="kingdom")
        self.assertTrue('Ranks in the first lineage are not in the correct order' in str(ctx.exception))
        for rec in self.invalid_tax[1:]:
            with self.assertRaises(Exception) as ctx:
                _do_filter([rec], defined_rank="kingdom")
            self.assertTrue('Not a valid QIIME-formatted lineage' in str(ctx.exception))

    def test_invalid_seq(self):
        # errors are only checked with ambig. filter
        with self.assertRaises(Exception) as ctx:
            _do_filter(self.invalid_seq, max_ambig=1)
        self.assertTrue('Invalid characters found in sequence: X, Z' in str(ctx.exception))


#### Entry point ####


if __name__ == '__main__':
    try:
        with file_logging(snakemake.log[0]):
            filter_taxdb(snakemake.input.all, snakemake.output.filtered, 
                        snakemake.input.cfg)
    except NameError:
        unittest.main()
