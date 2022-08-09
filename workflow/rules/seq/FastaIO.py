"""
Minimal FASTA reader and writer.
"""

from .seq import SeqRecord


def from_file(filepath, encoding='utf-8'):
    f = open(filepath, encoding=encoding)  # TODO: file not closed; cannot use with statement
    return parse(f)


def parse(handle):
    header = None
    for line in handle:
        line = line.rstrip('\r\n')
        if not line.startswith(";"):
            if line.startswith('>'):
                if header is not None:
                    yield _make_record(header, seq)
                header = line[1:].strip()
                seq = []
            elif header is not None:
                seq.append(line.rstrip())

    if header is not None:
        yield _make_record(header, seq)


def _make_record(header, seq):
    header = header.split(" ", 1)
    return SeqRecord(
        header[0],
        "".join(seq).replace(" ", ""),
        header[1].strip() if len(header) == 2 else ""
    )


def write(records, handle, **kwarg):
    if isinstance(records, SeqRecord):
        _write_fasta(records, handle, **kwarg)
    else:
        for r in records:
            _write_fasta(r, handle, **kwarg)


def to_file(records, filename, encoding='utf-8'):
    with open(filename, 'w', encoding=encoding) as out:
        write(records, out)


def _write_fasta(r, handle, wrap=None):
    seq = r.seq
    if wrap is not None:
        seq = "\n".join([seq[i:i+wrap] for i in range(0, len(seq), wrap)])
    desc = r.description if r.description else ""
    handle.write(">%s  %s\n%s\n" % (r.id, desc, seq))
