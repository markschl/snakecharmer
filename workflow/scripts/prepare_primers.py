import yaml
from seq_consensus import consensus
from utils import file_logging


__complements = bytes.maketrans(
    b'ATCGRYKMBVDHNSWatcgrykmbvdhnsw -',
    b'TAGCYRMKVBHDNSWtagcyrmkvbhdnsw -'
)


def reverse_complement(seq):
    return seq.translate(__complements)[::-1]


def parse_combine_primers(primers_by_marker):
    """
    Combines primers from different markers (in a nested dictionary structure)
    together (still grouped by forward/reverse), making sure that if the same
    primer name occurs in different markers, the primer sequences are the same
    (even though this may be a rare case).
    """
    out = {"forward": {}, "reverse": {}}
    for _marker, primers_by_dir in primers_by_marker.items():
        for dir_, primers in primers_by_dir.items():
            for name, seqs in primers.items():
                # remove anchors and split comma-delimited oligo lists
                anchor = False
                if seqs.startswith('^'):
                    anchor = True
                    seqs = seqs[1:]
                seqs = [s.strip() for s in seqs.split(',')]
                # check if seq already present
                if name in out[dir_]:
                    assert (seqs == out[dir_][name]['seq']), \
                          "Sequences of primer {} differ between different markers".format(name)
                else:
                    out[dir_][name] = {'seq': seqs, 'anchor': anchor}
    return out


def oligo_consensus(seqs, threshold=0.8):
    seqs = list(seqs)
    # we allow for different lengths by adding 3' terminal gaps
    # (oligos are expected to be aligned to 5', global alignment is a requirement
    # for many clustering algorithms)
    maxlen = max(len(s) for s in seqs)
    return consensus((s + '.' * (maxlen - len(s)) for s in seqs),
                     threshold=threshold,
                     free_endgaps=True,
                     end_gap_char='.')


def process_primers(primers_by_marker, single_method='consensus:50'):
    out = parse_combine_primers(primers_by_marker)
    # reverse complement of reverse sequence
    # (often needed for trimming reverse primer of merged amplicon)
    for dir_, seq_dict in list(out.items()):
        out[dir_ + '_rev'] = {
            name: {
                'seq': [reverse_complement(s) for s in d['seq']],
                'anchor': d['anchor']
            }
            for name, d in seq_dict.items()
        }
    # consensus
    if single_method.startswith('consensus:'):
        t = float(single_method.replace('consensus:', '').strip()) / 100.
        fn = lambda seqs: oligo_consensus(seqs, threshold=t)
    else:
        raise Exception('Currently no method other than "consensus:<threshold>" available to obtain a single oligo')
    for key, seq_dict in list(out.items()):
        out[key + '_consensus'] = {
            name: {
                'seq': fn(d['seq']),
                'anchor': d['anchor']
            }
            for name, d in seq_dict.items()
        }
    return out


with file_logging(snakemake.log[0]):
    out = process_primers(snakemake.params.primers)
    with open(snakemake.output.yaml, 'w') as f:
        yaml.safe_dump(out, f, sort_keys=False)
