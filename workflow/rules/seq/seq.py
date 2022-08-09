"""
Sequence representation
"""


try:
    from string import maketrans # python 2
except ImportError:
    maketrans = bytes.maketrans

ambig = {
    'N': {'A', 'T', 'G', 'C'},
    'R': {'A', 'G'},
    'S': {'G', 'C'},
    'K': {'T', 'G'},
    'B': {'T', 'G', 'C'},
    'H': {'A', 'T', 'C'},
    'Y': {'T', 'C'},
    'W': {'A', 'T'},
    'M': {'A', 'C'},
    'V': {'A', 'G', 'C'},
    'D': {'A', 'T', 'G'}
}
#ambig = {ord(char): {ord(char) for char in equals} for char, equals in _ambig.items()}

complements = maketrans(b'ATCGRYKMBVDHNSWatcgrykmbvdhnsw -',
                        b'TAGCYRMKVBHDNSWtagcyrmkvbhdnsw -')


def reverse_complement(seq):
    return seq.translate(complements)[::-1]


def variants(seq):
    """
    Returns all variants of a degenerated sequence (useful with primers)
    :param seq: DNA sequence with IUPAC codes for ambiguity
    :return: list of sequences without ambiguity
    """
    var = []
    for i, char in enumerate(seq):
        if char in ambig:
            for char2 in ambig[char]:
                var += variants(seq[:i] + char2 + seq[i+1:])
    if not var:
        var.append(seq)
    return var


class SeqRecord(object):
    #__slots__ = ['id', 'description', 'seq', 'qual']
    def __init__(self, id, seq, description='', qual=None):
        self.id = id
        self.description = description
        self.seq = seq
        self.qual = qual


