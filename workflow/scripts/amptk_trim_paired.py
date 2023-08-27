import os
import sys
import yaml
from subprocess import check_call

from utils import file_logging


def anchored(seq_d):
    if seq_d['anchor']:
        print(f"WARNING: The primer: {seq_d['seq']} should be anochored, but "
              "this is not possible with the Amptk workflow",
              file=sys.stderr)
    return seq_d['seq']


def trim_paired(input_files, demux_out,
                primer_file, f_primer_name, r_primer_name,
                program=None,
                usearch_bin=None,
                err_rate=None,
                min_len=None,
                threads=1):
    # read primers
    with open(primer_file) as f:
        primers = yaml.safe_load(f)
    f_primer = anchored(primers["forward_consensus"][f_primer_name])
    r_primer = anchored(primers["reverse_consensus"][r_primer_name])
    # determine mismatches (depending on average of forward/reverse primer lengths)
    # and subsequent rounding
    primer_mismatch = round((len(f_primer) + len(r_primer)) / 2 * err_rate)
    # NOTE: we limit the max. number of primer mismatches to 2.
    # The reason is that Amptk apparently does **another** primer trimming after 
    # merging the already trimmed reads.
    # Too liberal mismatch thresholds lead to many unspecific primer matches
    # and consequently to unwanted trimming of reads.
    if primer_mismatch > 2:
        print(f"WARNING: The maximum primer mismatches were limited to 2 with "
              "the Amptk workflow (would be {primer_mismatch} with current "
              "max_error_rate setting)",
              file=sys.stderr)
        primer_mismatch = 2
        
    # prepare input/output dirs
    # TODO: depends on this exact naming scheme
    prefix = demux_out.replace('.demux.fq.gz', '')
    outdir = os.path.dirname(prefix)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # Since we need to change the directory while executing, get the correct
    # relative input path
    indir = os.path.abspath(os.path.dirname(input_files[0]))
    # All files in the input dir are parsed by Amptk, therefore make sure that
    # there are no other files
    f1 = set(os.path.basename(f) for f in input_files)
    f2 = set(os.path.basename(f) for f in os.listdir(indir) if f.endswith('.fastq.gz'))
    assert f1 == f2, \
        "Amptk input dir has unknown files: {}".format(', '.join(f2.difference(f1)))

    # run
    cmd = [
        "amptk", "illumina",
        "-i", indir,
        "-o", os.path.basename(prefix),
        "-f", f_primer,
        "-r", r_primer,
        "--min_len", str(min_len),
        "--trim_len", "10000000000",  # high enough to never be longer
        "--cpus", str(threads),
        "--cleanup",
        "--require_primer=on",
        "--rescue_forward=off",
        "--primer_mismatch", str(primer_mismatch),
        "--merge_method", program
    ]
    if program == "usearch":
        assert usearch_bin is not None
        cmd += ["--usearch", usearch_bin]

    print("Call: " + " ".join(cmd), file=sys.stderr)
    check_call(cmd, cwd=outdir, stdout=sys.stdout, stderr=sys.stderr)


with file_logging(snakemake.log[0]):
    trim_paired(
        snakemake.input.fq,
        snakemake.output.demux,
        primer_file=snakemake.input.primers_yaml, 
        f_primer_name=snakemake.wildcards.f_primer,
        r_primer_name=snakemake.wildcards.r_primer,
        program=snakemake.params.program,
        usearch_bin=snakemake.params.usearch_bin,
        err_rate=snakemake.params.err_rate,
        min_len=snakemake.params.min_len,
        threads=snakemake.threads,
    )
