import sys
import yaml
from subprocess import check_call

from utils import file_logging


def trim_paired(primer_file, input, output, f_primer, r_primer, err_rate=None, threads=None, min_length=None):
    with open(primer_file) as f:
        primers = yaml.safe_load(f)
        
    anchored = lambda s: '^' + s['seq'] if s['anchor'] else s['seq']
    cmd = [
        "qiime", "cutadapt", "trim-paired",
        "--i-demultiplexed-sequences", input,
        "--p-cores", str(threads),
        "--p-adapter-f", "{fwd}...{r_rev};optional".format(
            fwd=anchored(primers['forward_consensus'][f_primer]),
            r_rev=primers['reverse_rev_consensus'][r_primer]['seq']
        ),
        "--p-adapter-r", "{rev}...{f_rev};optional".format(
            rev=anchored(primers['reverse_consensus'][r_primer]),
            f_rev=primers['forward_rev_consensus'][f_primer]['seq']
        ),
        "--p-error-rate", str(err_rate),
        "--p-overlap", "10", # TODO: configure
        "--p-minimum-length", str(min_length),
        "--p-discard-untrimmed",
        "--verbose",
        "--o-trimmed-sequences", output
    ]
    print("Call: " + " ".join(cmd), file=sys.stderr)
    check_call(cmd, stdout=sys.stdout, stderr=sys.stderr)


with file_logging(snakemake.log[0]):
    trim_paired(
        snakemake.input.yaml, 
        snakemake.input.demux, 
        snakemake.output.qza,
        f_primer=snakemake.wildcards.f_primer, 
        r_primer=snakemake.wildcards.r_primer,
        err_rate=snakemake.params.err_rate,
        min_length=snakemake.params.min_length,
        threads=snakemake.threads,
    )
