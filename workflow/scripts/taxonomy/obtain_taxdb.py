"""
This script obtains taxonomic databases from different sources
in the QIIME format.

All methods are expected to output the same ranks for all records 
(consistent number and order).

Ambiguous/missing information is represented by an empty string:
k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Sordariales;f__;g__;s__
# TODO: they may also be completely empty without rank char?
https://forum.qiime2.org/t/meaning-of-ambiguous-taxonomic-labels-in-greengenes/14428

Missing intermediate ranks may be empty. Rank name propagation to fill these
will be done automatically for classifiers that require it.
"""
from copy import copy
import csv
from io import TextIOWrapper
from itertools import zip_longest
import json
import os
from os.path import join, dirname, basename
import re
import sys
from tempfile import mkstemp
import unittest
from urllib.request import urlopen
import yaml

from tax_helpers import *
from utils import file_logging



def report_problem(permissive, assertion, text):
    """
    Warns or raises AssertionError upon problem depending on on permissive setting
    """
    if not assertion:
        if permissive:
            print(text, file=sys.stderr)
        else:
            assert assertion, text


class IdMismatch(Exception):
    pass


def read_aligned(fa_handle, tax_handle, full_desc=False):
    with fasta_reader(fa_handle, full_desc=full_desc) as seqs:
        lineages = csv.reader(tax_handle, delimiter="\t")
        for rec, tax in zip_longest(seqs, lineages):
            if rec is None or tax is None or rec.id != tax[0]:
                raise IdMismatch()
            yield rec, tax[1]


def read_inconsistent(fa_handle, tax_handle, full_desc=False):
    lineages = {l[0]: l[1] for l in csv.reader(tax_handle, delimiter="\t") if l}
    with fasta_reader(fa_handle, full_desc=full_desc) as seqs:
        for rec in seqs:
            try:
                yield rec, lineages[rec.id]
            except IndexError:
                raise Exception("Sequence ID '{}' not found in taxonomy".format(rec.id))


def report_count(recs):
    i = -1
    for i, r in enumerate(recs):
        yield r
    print("Obtained {} records".format(i + 1), file=sys.stderr)



class AmbigMatcher(object):
    """
    Helper class for matching ambiguous keywords/patterns (custom supplied) in files
    from unknown sources.
    """
    def __init__(self, patterns, sp_patterns=None):
        self.patterns = self.parse_config(patterns)
        self._sp_patterns_raw = sp_patterns
        self.sp_rank_i = None
        self.initialized = False

    def initialize(self, lineage):
        ranks = [r.strip().split("__", 1)[0] for r in lineage]
        self.empty_lineage = [rank + "__" for rank in ranks]
        self.sp_patterns = None
        if self._sp_patterns_raw is not None:
            self.sp_patterns = self.parse_config(self._sp_patterns_raw)
            try:
                self.sp_rank_i = ranks.index(RANK_TRANS['species'])
            except ValueError:
                print("Species ambiguity patterns specified, but no species in lineages",
                    file=sys.stderr)
        self.has_patterns = len(self.patterns) > 0 or self.sp_rank_i is not None

    @staticmethod
    def parse_config(cfg):
        assert isinstance(cfg, list)
        fixed = []
        regex = []
        for patt in cfg:
            if isinstance(patt, dict):
                assert len(patt) == 1 and next(iter(patt.keys())) == 'regex',\
                    "Invalid ambiguity pattern supplied"
                regex.append(re.compile(next(iter(patt.values()))))
            else:
                fixed.append(patt)
        return fixed, regex

    def is_ambig(self, rank_name, fixed, regex):
        # TODO: should the rank be removed from the name? only relevant for regex patterns matching '^'
        for pat in fixed:
            if pat in rank_name:
                return True
        for pat in regex:
            if pat.search(rank_name) is not None:
                return True
        return False

    def clean(self, lineage):
        # we need the first lineage to know which ranks to expect
        # and being able to fully initialize
        if not self.initialized:
            self.initialize(lineage)
            self.initialized = True
        # clean lineage
        if self.has_patterns:
            # regular patterns
            assert len(lineage) == len(self.empty_lineage), "lineage length mismatch"
            for i in range(len(lineage)):
                empty = self.empty_lineage[i]
                rank_name = lineage[i]
                assert rank_name.startswith(empty), "rank mismatch"
                if self.is_ambig(rank_name, *self.patterns):
                    lineage[i] = empty
            # species patterns
            if self.sp_rank_i is not None \
                and self.is_ambig(lineage[self.sp_rank_i], *self.sp_patterns):
                lineage[self.sp_rank_i] = self.empty_lineage[self.sp_rank_i]
        return lineage


def obtain_qiime_qza(target, sequences=None, taxonomy=None, ambig=[], ambig_sp=None, **other):
    """
    This simple implementation opens QIIME artifacts without requiring the QIIME2
    software (which would have to be loaded only for that specific method).
    see https://dev.qiime2.org/latest/storing-data/archive
    """
    fail_on_invalid(other)

    @contextmanager
    def _open_qza(path, type, filename, encoding="ascii", errors="replace"):
        with archive(path, "zip") as a:
            uuid = os.path.split(a.namelist()[0])[0]
            with a.open(uuid + '/metadata.yaml') as f, TextIOWrapper(f, encoding="ascii") as f:
                meta = yaml.safe_load(f)
                f = meta.get("format")
                assert f == type, "Unexpected format for {}: {}".format(path, f)
                assert meta["uuid"] == uuid
            with a.open("{}/data/{}".format(uuid, filename)) as f:
                yield TextIOWrapper(f, encoding=encoding, errors=errors)

    # usually, sequences and taxonomy should be aligned, but
    # we fall back to unordered parsing if the other fails
    with local_file(sequences) as seqfile, local_file(taxonomy) as taxfile:
        ambig_matcher = AmbigMatcher(ambig, ambig_sp)
        for i, read_fn in enumerate([read_aligned, read_inconsistent]):
            with _open_qza(seqfile, "DNASequencesDirectoryFormat", "dna-sequences.fasta") as seqs,\
                _open_qza(taxfile, "TSVTaxonomyDirectoryFormat", "taxonomy.tsv") as tax:
                try:
                    with zstd_fasta_writer(target) as out:
                        for rec, lineage in report_count(read_fn(seqs, tax)):
                            lineage = [r.strip() for r in lineage.split(";")]
                            ambig_matcher.clean(lineage)
                            rec.description = ";".join(lineage)
                            out.write_record(rec)
                except IdMismatch:
                    if i == 0:
                        print("Sequences and taxonomy not aligned, retrying using taxonomy dictionary...",
                              file=sys.stderr)
                    else:
                        raise Exception("Sequence ID mismatch")



def obtain_qiime(target, file=None, ambig=[], ambig_sp=None, **other):
    """
    Converter for the flat-file FASTA format with QIIME-style lineages
    in the header.
    This essentially only cleans the lineages (if any ambiguity keywords were supplied)
    """
    fail_on_invalid(other)
    with zstd_fasta_writer(target) as out, fasta_reader(file) as recs:
        filter_ambig = ambig or ambig_sp
        ambig_matcher = AmbigMatcher(ambig, ambig_sp)
        for rec in report_count(recs):
            if filter_ambig:
                lineage = [r.strip() for r in rec.description.split(";")]
                ambig_matcher.clean(lineage)
                rec.description = ";".join(lineage)
            # now we can collapse and write the output
            out.write_record(rec)


def obtain_utax(target, file=None, ambig=[], ambig_sp=None, ignore_problems=False, **other):
    """
    UTAX format converter
    Since the UTAX format allows for missing ranks in part of the lineages,
    we have to infer these first, which makes the code somewhat more
    complicated.
    """
    fail_on_invalid(other)

    # helper functions
    def set_remove(s, key):
        try:
            s.remove(key)
        except KeyError:
            return False
        return True
    
    def split_rank(rank):
        s = rank.strip().split(':', 1)
        assert len(s) == 2, "Invalid rank: {}".format(rank)
        return s
    
    lineage_pat = re.compile(r"(.+?);tax=(.+?)\s*;?\s*")
    def get_id_lineage(lineage):
        try:
            return lineage_pat.fullmatch(lineage).groups()
        except AttributeError:
            raise Exception("Not an UTAX-formatted lineage: {}".format(lineage))

    with local_file(file) as path:
        # obtain a list of ranks
        print("Looking for taxonomic ranks...", file=sys.stderr)
        with fasta_reader(path, full_desc=True) as recs:
            ranks = set(r.strip().split(':', 1)[0]
                        for rec in recs
                        for r in get_id_lineage(rec.description)[1].split(','))
            ranks.discard('')
            assert len(ranks) > 0, "The UTAX lineages are empty"
            all_ranks = [r for r in RANK_TRANS.values() if set_remove(ranks, r)]
            assert len(ranks) == 0, \
                ("Some ranks are not known and this method is not smart enough "
                 "to know where to put them: {}").format(", ".join(ranks))

        # now we can read the file again, converting the lineages
        print("Converting lineages...", file=sys.stderr)
        ambig_matcher = AmbigMatcher(ambig, ambig_sp)
        with zstd_fasta_writer(target) as out, fasta_reader(path, full_desc=True) as recs:
            reserved = re.compile(r"(__|;)")
            empty_lineage = [r + '__' for r in all_ranks]
            for rec in report_count(recs):
                rec.id, lineage = get_id_lineage(rec.description)
                # replace reserved chars
                lineage = reserved.sub('_', lineage)
                # fill in ranks
                lineage_out = copy(empty_lineage)
                rank, name = None, None
                rank_iter = iter(lineage.split(','))
                item = next(rank_iter)
                rank, name = split_rank(item)
                for i, exp_rank in enumerate(all_ranks):
                    if rank == exp_rank:
                        # if rank matches -> insert it in the lineage and
                        # proceed to next
                        lineage_out[i] = '{}__{}'.format(rank, name)
                        try:
                            item = next(rank_iter)
                            rank, name = split_rank(item)
                        except StopIteration:
                            rank, name = None, None
                report_problem(ignore_problems, rank is None, (
                    "The order of ranks was not as expected: {} "
                    "(found {}:{})").format(rec.description, rank, name))

                # then we clean ambiguous names (if configured so)
                ambig_matcher.clean(lineage_out)

                # now we can collapse and write the output
                rec.description = ";".join(lineage_out)
                out.write_record(rec)


def obtain_gtdb(target, file=None, **other):
    fail_on_invalid(other)
    assert file is not None and isinstance(file, (str, list)),\
        "'file' in GTDB config should be list or single string"
    if isinstance(file, str):
        file = [file]
    with zstd_fasta_writer(target) as out:
        desc_pat = re.compile(r" \[\w+=.+?\].*")
        for f in file:
            print(f"Obtaining {file}", file=sys.stderr)
            with archive(f, format="tar") as t:
                # parse the content and write to output
                f = t.getnames()
                assert len(f) == 1
                with TextIOWrapper(t.extractfile(f[0]), "ascii", errors="replace") as fa,\
                    fasta_reader(fa) as seqs:
                    for rec in report_count(seqs):
                        # remove [...=...] annotations from end
                        rec.description = desc_pat.sub('', rec.description)
                        out.write_record(rec)


def obtain_unite_otus(target, doi=None, file=None, date=None, threshold="dynamic", kind="regular", **other):
    fail_on_invalid(other)
    if file is None:
        # Query the UNITE API using the DOI in order to obtain the file URL
        assert doi is not None, "Reference databases of format 'unite' need a 'doi' or an 'url' specified"
        assert not doi.startswith("http"), "Please specify the UNITE DOI without preceding https://doi.org/"
        res = urlopen("https://api.plutof.ut.ee/v1/public/dois/?format=vnd.api%2Bjson&identifier=" + doi).read()
        d = json.loads(res)
        files = [(f["name"], f["url"]) for f in d["data"][0]["attributes"]["media"]]
        if date is not None:
            date = str(date)
            filtered = [f for f in files if date in f[0]]
            assert len(filtered) == 1, \
                "UNITE file 'date' setting does not identify one file. Available: {}".format(", ".join(f for f, _ in files))
            file = filtered[0][1]
        else:
            if len(files) > 1:
                print("WARNING: selected the last file in the list: {}. Available are: {}. Set 'date' to select another file.".format(
                    files[-1][0], ", ".join(f for f, _ in files)
                ), file=sys.stderr)
            file = files[-1][1]
    
    # obtain the file
    print(f"Obtaining {file}", file=sys.stderr)
    with archive(file, format="tar") as t:
        # select the correct files
        if kind == "regular":
            dname = ""
        else:
            assert kind == "developer", "Invalid database 'kind' (set to 'regular' or 'developer')"
            dname = "developer"
        files = [f for f in t.getnames() if dirname(f) == dname]
        file_pat = re.compile(r"sh_[a-z]+_qiime_ver\d+_([a-z0-9]+)_.+?\.(fasta|txt)")
        matches = ((file_pat.fullmatch(basename(f)), f) for f in files)
        file_info = {m.groups(): f for m, f in matches if m is not None}
        threshold = str(threshold)
        try:
            fa = file_info[(threshold, "fasta")]
            tax = file_info[(threshold, "txt")]
        except KeyError:
            raise Exception((
                "Correct UNITE FASTA/taxonomy files not found, check values for 'threshold' and 'kind'. "
                "If correct, there may be a bug to report. Available files are: {}".format(", ".join(files))))
        print(f"Reading from {fa} and {tax}", file=sys.stderr)

        # precompile regex patterns to replace
        sp_pattern = re.compile(r"s__\w+?[_ ]sp\.?")
        unknown_pattern = re.compile(r"([a-z])__(\w+?_[a-z]{3}_Incertae_sedis|unidentified|)")

        # parse FASTA and add lineages
        with zstd_fasta_writer(target) as out:
            with TextIOWrapper(t.extractfile(tax), "ascii", errors="replace") as lineages, \
                TextIOWrapper(t.extractfile(fa), "ascii", errors="replace") as seqs:
                for rec, lineage in report_count(read_inconsistent(seqs, lineages)):
                    lineage = [n.strip() for n in lineage.split(";")]
                    # make undefined ranks empty
                    if sp_pattern.fullmatch(lineage[-1]) is not None:
                        lineage[-1] = "s__"
                    for i in range(len(lineage)):
                        m = unknown_pattern.fullmatch(lineage[i])
                        if m is not None:
                            lineage[i] = "{}__".format(m.group(1))
                    rec.description = ";".join(lineage)
                    out.write_record(rec)


def obtain_midori(target, prefix=None, version=None, marker=None, kind=None, include_ambig=None, remove_num=None, **other):
    fail_on_invalid(other)
    if prefix is None:
        assert version is not None and marker is not None and kind is not None, "Reference databases of format 'midori' need 'version', 'marker' and 'kind' defined, or alternatively an 'url_prefix'"
        assert kind in {'longest', 'uniq'}, "Reference databases of format 'midori' need kind=uniq or kind=longest"
        base_url = "https://www.reference-midori.info/forceDownload.php?fName=download/Databases/GenBank{}/".format(version)
        if include_ambig is True:
            prefix = base_url + "QIIME_sp/{}/MIDORI2_{}_NUC_SP_GB{}_{}_QIIME".format(kind, kind.upper(), version, marker)
        else:
            prefix = base_url + "QIIME/{}/MIDORI2_{}_NUC_GB{}_{}_QIIME".format(kind, kind.upper(), version, marker)

    # obtain the files
    print(f"Obtaining from: {prefix}", file=sys.stderr)
    seqfile = prefix + ".fasta.gz"
    taxfile = prefix + ".taxon.gz"

    # compile regex patterns
    num_pat = re.compile(r"_\d+(;|$)")
    undef_pat = {
        name[0]: re.compile(r"{}__{}_.+".format(name[0], name))
        for name in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus']
    }
    # Non-filtered version including sp. and others defined only at higher ranks
    # -> we have to convert undefined names to empty strings.
    # The patterns follow information from https://doi.org/10.1002/edn3.303 and comparisons with
    # filtered files from Midori. Still, it was not possible to obtain the exact
    # same result, differences are however very small (see validation in test_data/taxonomy/midori).
    sp_pattern=r"""
    [_\s]
    (
        (
            cf\.|aff\.|sp\.|environment|undescribed|uncultured|complex|unclassified
            |nom\.|_nomen\s.*|_nom\.\s.*|nud\.|unidentif\.|indet\.|gen\.|nr\.
            |taxon\s\w+
        )
        [_\s\d]
        |(sp|cf)[\._][A-Z0-9]
        |sp\.$
    )
    """
    sp_pattern = re.compile(sp_pattern, re.VERBOSE)

    with textfile(seqfile, gz=True) as fa, textfile(taxfile, gz=True) as tax:
        # then filter and write to output
        with zstd_fasta_writer(target) as out:
            for rec, lineage in report_count(read_aligned(fa, tax)):
                if remove_num is True:
                    lineage = num_pat.sub(r"\1", lineage)
                lineage = lineage.split(";")
                # Set undefined names to empty strings
                # Set species names matching one of the above patterns to undefined
                # The species identified as undefined here vs. species filtered
                # from the QIIME_SP file (=QIIME version) do not match perfectly,
                # but almost (see validation code above)
                assert lineage[-1].startswith("s__")
                if sp_pattern.search(lineage[-1]) is not None:
                    lineage[-1] = "s__"
                # remove _family_..., _phylum_..., etc.
                for i in range(len(lineage)-1):
                    taxon = lineage[i]
                    rank_char = taxon[0]
                    if undef_pat[rank_char].fullmatch(taxon) is not None:                    
                        lineage[i] = taxon[:3]
                if lineage[-2] == "g__" and lineage[-1] != "s__" and not lineage[-1].startswith("s__["):
                    print(f"WARNING: species was not recognized as ambiguous, "
                          "but the genus is: {lineage[-1]}. "
                          "The species was set to undefined.",
                           file=sys.stderr)
                    lineage[-1] = "s__"
                rec.description = ";".join(lineage)
                out.write_record(rec)


def obtain_taxdb(param_file, outfile):
    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    with open(param_file) as f:
        d = yaml.safe_load(f)

    try:
        fmt = d.pop("format")
    except:
        raise Exception("The database format ('format' keyword) must be defined for {}".format(d["name"]))

    if fmt in ("unite_otus", "midori", "gtdb", "utax", "qiime", "qiime_qza"):
        func = getattr(sys.modules[__name__], "obtain_" + fmt)
        func(target=outfile, **d)
    else:
        raise Exception("Unknown taxonomy database format: {}".format(fmt))


#### Tests ####


def _obtain_fasta(expected_file=None, **cfg):
    outfile = mkstemp()[1]
    cfg_file = mkstemp()[1]
    with open(cfg_file, "w") as f:
        yaml.safe_dump(cfg, f)
    obtain_taxdb(cfg_file, outfile)
    expected = None
    if expected_file:
        with open(expected_file) as f:
            expected = f.read().strip('\n')
    with zstd_reader(outfile) as f:
        found = f.read().strip('\n')
    os.remove(outfile)
    os.remove(cfg_file)
    return expected, found


class Tester(unittest.TestCase):
    data = join("test_data", "taxonomy")

    def test_unite(self):
        """
        The file contains a small subset of sequences downloaded here:
        https://unite.ut.ee/repository.php
        (licensed CC BY-SA 4.0)
        """
        unite = join(self.data, "unite")
        exp, found = _obtain_fasta(
            join(unite, "v97.fasta"), 
            file=join(unite, "unite_test.tgz"),
            format="unite_otus",
            threshold=97
        )
        self.assertEqual(exp, found)
        unite = join(self.data, "unite")
        exp, found = _obtain_fasta(
            join(unite, "v97_dev.fasta"), 
            file=join(unite, "unite_test.tgz"),
            format="unite_otus",
            threshold=97,
            kind="developer"
        )
        self.assertEqual(exp, found)
        # code for checking download functionality:
        # exp, found = _obtain_fasta(
        #     doi="10.15156/BIO/2483918",
        #     format="unite_otus",
        #     date="29.11.2022",
        #     threshold=97
        # )

    def test_gtdb(self):
        """
        The file contains a small subset of sequences downloaded here:
        https://data.gtdb.ecogenomic.org/releases
        (licensed CC BY-SA 4.0)
        """
        path = join(self.data, "gtdb")
        exp, found = _obtain_fasta(
            join(path, "gtdb.fna"), 
            file=join(path, "gtdb.tar.gz"),
            format="gtdb",
        )
        self.assertEqual(exp, found)

    def test_midori(self):
        """
        The files contain a small subset of sequences downloaded here:
        https://www.reference-midori.info
        """
        path = join(self.data, "midori")
        exp, found = _obtain_fasta(
            join(path, "clean.fasta"),
            prefix=join(path, "MIDORI2_LONGEST_NUC_SP_GB256_lrRNA_QIIME"),
            format="midori",
        )
        self.assertEqual(exp, found)
        exp, found = _obtain_fasta(
            join(path, "clean_nonum.fasta"),
            prefix=join(path, "MIDORI2_LONGEST_NUC_SP_GB256_lrRNA_QIIME"),
            format="midori",
            remove_num=True
        )
        self.assertEqual(exp, found)
        # _obtain_fasta(
        #     format="midori",
        #     version=256,
        #     marker="lrRNA",
        #     kind="longest",
        #     remove_num=True,
        #     include_ambig=True
        # )

    def test_utax(self):
        """
        The file contains a small subset of sequences downloaded here:
        https://www.drive5.com/usearch/manual/sintax_downloads.html
        """
        path = join(self.data, "utax")
        exp, found = _obtain_fasta(
            join(path, "qiime_formatted.fasta"),
            file=join(path, "rdp_16s_v18.fa.gz"),
            format="utax",
        )
        self.assertEqual(exp, found)
        # _obtain_fasta(
        #     file="https://www.drive5.com/sintax/silva_18s_v123.fa.gz",
        #     format="utax",
        #     ignore_problems=True,
        #     ambig=["uncultured"]
        # )


    def test_qiime_qza(self):
        """
        The file contains a small subset of sequences downloaded here:
        https://docs.qiime2.org/2023.5/data-resources/#data-resources
        The QZA archives are stripped down and don't contain provenience/md5 hashes
        """
        path = join(self.data, "qiime_qza")
        exp, found = _obtain_fasta(
            join(path, "qiime_flat.fasta"),
            sequences=join(path, "silva-138-99-seqs-515-806.qza"),
            taxonomy=join(path, "silva-138-99-tax-515-806.qza"),
            format="qiime_qza",
        )
        self.assertEqual(exp, found)
        exp, found = _obtain_fasta(
            join(path, "qiime_flat.fasta"),
            sequences=join(path, "silva-138-99-seqs-515-806.qza"),
            taxonomy=join(path, "silva-138-99-tax-515-806_reordered.qza"),
            format="qiime_qza",
        )
        self.assertEqual(exp, found)
        exp, found = _obtain_fasta(
            join(path, "qiime_flat_noambig.fasta"),
            sequences=join(path, "silva-138-99-seqs-515-806.qza"),
            taxonomy=join(path, "silva-138-99-tax-515-806.qza"),
            format="qiime_qza",
            ambig=["uncultured"],
            ambig_sp=[{"regex": "[ _]sp\.?$"}]
        )
        self.assertEqual(exp, found)
        # _obtain_fasta(
        #     sequences="http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.fna.qza",
        #     taxonomy="http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.tax.qza",
        #     format="qiime_qza",
        #     ambig=["uncultured"],
        #     ambig_sp=[{"regex": "[ _]sp\.?$"}]
        # )

    def test_qiime(self):
        """
        Flat QIIME-style annotations in FASTA descriptions
        Same data as in test_qiime_qza.
        """
        path = join(self.data, "qiime")
        exp, found = _obtain_fasta(
            join(path, "qiime_out1.fasta"),
            file=join(path, "qiime_input.fasta"),
            format="qiime",
            ambig=["uncultured"]
        )
        self.assertEqual(exp, found)
        exp, found = _obtain_fasta(
            join(path, "qiime_out2.fasta"),
            file=join(path, "qiime_input.fasta"),
            format="qiime",
            ambig=["uncultured", "human_gut"],
            ambig_sp=[{"regex": "sp\.?$"}]
        )
        self.assertEqual(exp, found)

    def test_ambig_matcher(self):
        m = AmbigMatcher([], None)
        assert m.clean(["s__undefined"]) == ["s__undefined"]
        m = AmbigMatcher(["undefined"], None)
        assert m.clean(["k__Fungi", "s__undefined"]) == ["k__Fungi", "s__"]
        assert m.clean(["k__Fungi", "s__uncultured"]) == ["k__Fungi", "s__uncultured"]
        m = AmbigMatcher(["undefined"], [" sp."])
        assert m.clean(["k__undefined", "s__undefined sp."]) == ["k__", "s__"]
        assert m.clean(["k__Fungi", "s__Some sp.345"]) == ["k__Fungi", "s__"]
        m = AmbigMatcher([], [{"regex": "\ssp\.$"}])
        assert m.clean(["k__undefined", "g__Genusp."]) == ["k__undefined", "g__Genusp."]
        m = AmbigMatcher([], [{"regex": "\ssp\.$"}])
        assert m.clean(["k__undefined", "s__undefined sp."]) == ["k__undefined", "s__"]
        assert m.clean(["k__Fungi", "s__Some sp.345"]) == ["k__Fungi", "s__Some sp.345"]
        m = AmbigMatcher([], [{"regex": "\ssp\.$"}])
        m.clean(["k__"])
        with self.assertRaises(Exception) as ctx:
            m.clean(["k__", "s__"])
        self.assertTrue("length mismatch" in str(ctx.exception))
        m = AmbigMatcher([])
        m.clean(["k__Kingdom", "g__Genus"])
        with self.assertRaises(Exception) as ctx:
            m.clean(["k__Kingdom", "s__Genus"])
        self.assertTrue("rank mismatch" in str(ctx.exception))


#### Entry point ####


if __name__ == '__main__':
    try:
        with file_logging(snakemake.log[0]):
            obtain_taxdb(snakemake.input.yml, snakemake.output.db)
    except NameError:
        unittest.main()
