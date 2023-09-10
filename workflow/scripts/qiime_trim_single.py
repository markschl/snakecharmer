import sys
import yaml
from subprocess import check_call

from utils import file_logging


def trim_paired(primer_file, input, output, f_primer, r_primer, rev_read=False, err_rate=None, threads=None, min_length=None):
    with open(primer_file) as f:
        primers = yaml.safe_load(f)
    
    anchored = lambda s: '^' + s['seq'] if s['anchor'] else s['seq']
    if not rev_read:
        f_seq = anchored(primers['forward_consensus'][f_primer])
        r_rev = primers['reverse_rev_consensus'][r_primer]['seq']
    else:
        f_seq = anchored(primers['reverse_consensus'][r_primer])
        r_rev = primers['forward_rev_consensus'][f_primer]['seq']

    # see https://docs.qiime2.org/2023.7/plugins/available/cutadapt/trim-single/
    cmd = [
        "qiime", "cutadapt", "trim-single",
        "--i-demultiplexed-sequences", input,
        "--p-cores", str(threads),
        "--p-adapter", f"{f_seq}...{r_rev};optional",
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
        rev_read=snakemake.params.rev_read,
        threads=snakemake.threads,
    )
