from collections import OrderedDict
from contextlib import contextmanager
import gzip
from io import IOBase, TextIOWrapper
import os
import re
import tarfile
from tempfile import mkstemp
from urllib.request import urlopen, urlretrieve
import zipfile

import zstandard as zstd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser, FastaWriter


# Possible taxonomic ranks with their corresponding QIIME prefix,
# following https://github.com/bokulich-lab/RESCRIPt
RANK_TRANS = OrderedDict([
    ('domain', 'd'),
    ('superkingdom', 'sk'),
    ('kingdom', 'k'),
    ('subkingdom', 'ks'),
    ('superphylum', 'sp'),
    ('phylum', 'p'),
    ('subphylum', 'ps'),
    ('infraphylum', 'pi'),
    ('superclass', 'sc'),
    ('class', 'c'),
    ('subclass', 'cs'),
    ('infraclass', 'ci'),
    ('cohort', 'co'),
    ('superorder', 'so'),
    ('order', 'o'),
    ('suborder', 'os'),
    ('infraorder', 'oi'),
    ('parvorder', 'op'),
    ('superfamily', 'sf'),
    ('family', 'f'),
    ('subfamily', 'fs'),
    ('tribe', 't'),
    ('subtribe', 'ts'),
    ('genus', 'g'),
    ('subgenus', 'gs'),
    ('species group', 'ss'),
    ('species subgroup', 'sgs'),
    ('species', 's'),
    ('subspecies', 'ssb'),
    ('forma', 'for')
])


def _fasta_reader(handle, full_desc=False):
    """
    Custom FASTA parser that doesn't include the IDs in the description
    """
    for title, seq in SimpleFastaParser(handle):
        try:
            id, desc = title.split(" ", 1)
        except ValueError:
            id, desc = title, ""
        if full_desc:
            desc = title
        yield SeqRecord(Seq(seq), id=id, description=desc)


@contextmanager
def zstd_writer(filename, **wrapper_args):
    with open(filename, "wb") as outh:
        cctx = zstd.ZstdCompressor()
        with cctx.stream_writer(outh) as out:
            with TextIOWrapper(out, **wrapper_args) as w:
                yield w


@contextmanager
def zstd_fasta_writer(filename, encoding="ascii", errors="replace", **kwargs):
    with zstd_writer(filename, encoding=encoding, errors=errors, **kwargs) as out:
        yield FastaWriter(out, wrap=None)


@contextmanager
def zstd_reader(file, encoding="ascii", errors="replace", **kwargs):
    with bin_file(file, gz='auto', **kwargs) as inh:
        cctx = zstd.ZstdDecompressor()
        with cctx.stream_reader(inh) as f:
            with TextIOWrapper(f, encoding=encoding, errors=errors) as r:
                yield r

@contextmanager
def zstd_fasta_reader(filename, full_desc=False, **kwargs):
    with zstd_reader(filename, encoding="ascii", errors="replace", **kwargs) as f:
        yield _fasta_reader(f, full_desc=full_desc)
        

def _get_url_path(file=None, path=None, url=None):
    if file is not None:
        if re.match("https?://", file):
            url = file
        else:
            path = file
    assert not (path is not None and url is not None), \
        "Both path and url are defined, which one should be used?"
    return path, url


@contextmanager
def bin_file(file=None, path=None, url=None, gz="auto"):
    path, url = _get_url_path(file, path, url)
    if url is not None:
        with urlopen(url) as resp:
            if gz is True or gz == "auto" and resp.info().get("Content-Encoding") == "gzip":
                with gzip.GzipFile(fileobj=resp) as gz:
                    yield gz
            else:
                yield resp
    else:
        assert path is not None, "Either path or url should be defined"
        with open(path, "rb") as f:
            if gz is True or gz == "auto" and f.read(2) == b"\x1f\x8b":
                f.seek(0)
                with gzip.open(f) as gz:
                    yield gz
            else:
                assert gz in (False, "auto")
                f.seek(0)
                yield f


@contextmanager
def textfile(file=None, encoding='ascii', errors='replace', **settings):
    with bin_file(file=file, **settings) as f:
        with TextIOWrapper(f, encoding=encoding, errors=errors) as w:
            yield w


@contextmanager
def fasta_reader(file=None, full_desc=False, **settings):
    if isinstance(file, IOBase):
        yield _fasta_reader(file, full_desc=full_desc)
    else:
        with textfile(file=file, **settings) as f:
            yield _fasta_reader(f, full_desc=full_desc)


@contextmanager
def local_file(file=None, path=None, url=None):
    path, url = _get_url_path(file, path, url)
    if path is not None:
        yield path
    if url is not None:
        tmp = mkstemp()[1]
        urlretrieve(url, tmp)
        yield tmp
        os.remove(tmp)


@contextmanager
def archive(path=None, format="tar", **kwargs):
    with local_file(file=path, **kwargs) as f:
        if format == "tar":
            _open = tarfile.open
        else:
            assert format == "zip"
            _open = zipfile.ZipFile
        with _open(f) as out:
            yield out
