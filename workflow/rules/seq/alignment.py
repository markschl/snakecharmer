#!/usr/bin/python



from .seq import ambig


rev_ambig = {}
for k, v in ambig.items():
    rev_ambig[''.join(sorted(v))] = k


def alignment_map(records, endgaps=False, positions=None):
    """
    Parses an alignment and returns information about character occurrences
    at each position.
    
    :param records: iterable of :class:`seq.SeqRecord` objects with *aligned* 
        sequences.
    :param endgaps: End gaps will not be counted as gaps unless this option is True. 
        Thus, it is possible that empty dictionaries are returned if alignments have 
        positions with gaps only at the beginning/end
    :returns: list of dictionaries, where each dict contains the character counts
        at the given position in the alignment.

    Example::
        >> seqs = [
            seq.SeqRecord('r1', 'ATTGC'),
            seq.SeqRecord('r2', 'AT-CC'),            
        ]
        >> alignment_map(seqs)
        [{'A': 2}, {'T': 2}, {'-': 1, 'T': 1}, {'C': 1, 'G': 1}, {'C': 2}]
    """
    occurrences = None

    try:
        for r, record in enumerate(records):
            if r == 0:
                seqlen = len(record.seq)
                rng = positions or range(seqlen)
                occurrences = [{} for _ in rng]
            seq = (record.seq if endgaps else _remove_endgaps(record.seq)).upper()
            for column in rng:
                char = seq[column]
                if char != ' ' or endgaps:
                    try:
                        occurrences[column][char] += 1
                    except KeyError:
                        occurrences[column][char] = 1
                    
    except IndexError:
        raise Exception(
            "Sequence length of '%s' does not match other sequences (%d != %d)! Aligned sequences need to have the same length." % (
             record.id, len(record.seq), len(occurrences))
        )
    return occurrences


def compare(seq1, seq2, compare_endgaps=False, ignore_gaps=False, case_insensitive=True):
    if not compare_endgaps:
        seq1, seq2 = _remove_endgaps(seq1), _remove_endgaps(seq2)
    if case_insensitive:
        seq1 = seq1.upper()
        seq2 = seq2.upper()
    n = min(len(seq1), len(seq2))
    matching, counted = bytearray(n), bytearray(n)
    for i in range(n):
        if not seq1[i] == seq2[i] == '-':
            if matches(seq1[i], seq2[i], ignore_gaps=ignore_gaps):
                matching[i] += 1
            counted[i] += 1
    return list(matching), list(counted)


def n_mismatches(seq1, seq2, **param):
    matching, counted = compare(seq1, seq2, **param)
    return len(matching) - sum(matching)


def fast_match(seq1, seq2, compare_endgaps=False, ignore_gaps=False):
    if seq1 == seq2:
        return True
    if not compare_endgaps:
        seq1, seq2 = _remove_endgaps(seq1), _remove_endgaps(seq2)
    n = min(len(seq1), len(seq2))
    for i in range(n):
        if not seq1[i] == seq2[i] == '-' and not matches(seq1[i], seq2[i], ignore_gaps=ignore_gaps):
                return False
    return True


def matches(char1, char2, ignore_gaps=False):
    """

    Args:
        char1:
        char2: must be contained in char2 (A matches R, but not other way round)
        ignore_gaps:

    Returns:

    """
    return char1 == char2 or \
           char2 in ambig and char1 in ambig[char2] or \
           char1 in ambig and char2 in ambig and ambig[char1].issubset(ambig[char2]) or \
           char1 == ' ' or char2 == ' ' or \
           ignore_gaps and (char1 == '-' or char2 == '-')


def consensus(records, threshold=0, endgaps=False, ignore_gaps=False):
    occurrences = alignment_map(records, endgaps)
    if occurrences is None:
        raise Exception('Alignment is empty')
    return ''.join(get_consensus_char(o, threshold, ignore_gaps=ignore_gaps)
                   for o in occurrences)


def _remove_endgaps(seq):
    l = len(seq)
    s = seq.lstrip('-')
    if len(s) < l:
        seq = ' ' * (l - len(s)) + s
    s = seq.rstrip('-')
    if len(s) < l:
        seq = s + ' ' * (l - len(s))
    return seq


def get_consensus_char(occurrences, threshold=0, ignore_gaps=False):
    """
    :param occurrences: dict with base codes as keys and counts as values
    :param threshold: number between 0 and 1 denoting the frequency a base needs
        in order to be chosen as consensus. If the frequency is below the threshold,
        the consensus will be an ambiguity code representing a combination of bases
        whose cumulative frequency is above the threshold. With the default threshold
        (0), the most frequent base will always be chosen as consensus.
    :param ignore_gaps: Do not care about gaps at all, only return bases if found
    """
    # replace ambiguous
    for char in list(occurrences.keys()):
        try:
            for char2 in ambig[char]:
                occurrences[char2] += occurrences[char]
            del occurrences[char]
        except KeyError:
            pass

    threshold = threshold * sum(occurrences.values())
    charlist = list(zip(list(occurrences.values()), list(occurrences.keys())))
    charlist = iter(sorted(charlist, reverse=True))
    nn = 0
    chars = []
    for n, char in charlist:
        nn += n
        chars.append(char)
        if nn >= threshold:
            # add ties
            for next_n, next_char in charlist:
                if n != next_n:
                    break
                chars.append(next_char)                
            if len(chars) == 1:
                return char
            if '-' in chars:
                if ignore_gaps:
                    chars.remove('-')
                    if len(chars) == 1:
                        return chars[0]
                else:
                    return 'N'

            return _get_ambig(chars) if chars else '-'
    
    return _get_ambig(chars) if chars else '-'


def _get_ambig(chars):
    """
    :param chars: list/tuple of base codes (A, T, G or C)
    :returns: IUPAC ambiguous code corresponding to the given base combination
    """
    if not chars:
        raise Exception('Empty character sequence!')
    try:
        return rev_ambig[''.join(sorted(chars))]
    except KeyError:
        return 'N'

