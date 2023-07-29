from os.path import dirname, basename
import csv
import sys
from collections import defaultdict


def write_stats(merge_stats, trim_stats, filter_stats, primer_combinations, outfile):
    # read data
    merge = read_stats(merge_stats, key=lambda f: basename(dirname(f)))
    trim = read_stats(trim_stats, key=lambda f: basename(dirname(dirname(f))))
    flt = read_stats(
        filter_stats,
        key=lambda f: (basename(dirname(dirname(f))), basename(dirname(f))),
    )
    flt2 = defaultdict(dict)
    for k, n in flt.items():
        flt2[k[1]][k[0]] = n

    # combine and write to output
    with open(outfile, 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        header = ['sample', 'raw', 'merged', '(%)', 'fwd-primer', '(% of merged)', 'rev-primer', '(% of merged)', 'long enough', '(% of merged)']
        for p in primer_combinations:
            header += [p, "filtered", "(% of trimmed)"]
        writer.writerow(header)
        for sample in merge:
            n_raw, n_merged = merge[sample]
            m2, n_fwd, n_rev, n_long_enough = trim[sample]
            assert (n_merged == m2), "Number of sequences in logfiles from read merging and primer trimming does not match: {} vs. {}".format(n_merged, m2)
            row = [
                sample, 
                n_raw,
                n_merged, percent(n_merged, n_raw),
                n_fwd, percent(n_fwd, n_merged),
                n_rev, percent(n_rev, n_merged),
                n_long_enough, percent(n_long_enough, n_merged)
            ]
            f = flt2[sample]
            for p in primer_combinations:
                kept, removed = f[p]
                total = kept + removed
                row += [total, kept, percent(kept, total)]
            writer.writerow(row)


def read_stats(files, key):
    out = {}
    for f in files:
        with open(f, "r") as h:
            out[key(f)] = [int(n) for n in next(csv.reader(h, delimiter="\t"))]
    return out


def percent(x, y):
    return round(100 * x / y, 2) if y > 0 else 0.0



with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    write_stats(snakemake.input.merge, snakemake.input.trim, snakemake.input.filter,
                snakemake.params.primer_combinations,
                snakemake.output[0])
