import gzip
import os
import shutil
import sys
from subprocess import check_call

from utils import file_logging


def cluster_paired(method,
                trimmed_in,
                clustered_out, otutab_out,
                usearch_par,
                dada2_par,
                unoise_program=None,
                usearch_bin=None,
                threads=1):
    
    cmd = [
        "amptk", method,
        "-i", os.path.basename(trimmed_in),
        "-o", method,
        "--cpus", str(threads)
    ]

    # add method specific arguments
    # TODO: inconsistent to have these settings under usearch even though used for DADA2 as well
    # (on the other hand, Amptk uses an USEACH-like procedure for DADA2 as well)
    maxee = usearch_par["merge"]["expected_length"] * usearch_par["filter"]["max_error_rate"]
    if method == "unoise3":
        cmd += [
                "--usearch",
                usearch_bin,
                "--maxee",
                str(maxee),
                "--method",
                unoise_program,
                "--minsize",
                str(usearch_par["unoise3"]["min_size"]),
            ]
        if unoise_program == "usearch":
            assert usearch_bin is not None
            cmd += ["--usearch", usearch_bin]
    else:
        assert method == "dada2", "Unknown / unimplemented Amptk command"
        cmd += [
            "--maxee",
            str(maxee),
            "--chimera_method",
            dada2_par["chimera_method"],
        ]
        p = dada2_par["pooling_method"]
        if p == "pooled":
            cmd.append("--pool")
        elif p == "pseudo":
            cmd.append("--pseudopool")

    outdir = os.path.dirname(trimmed_in)
    print("Call: " + " ".join(cmd), file=sys.stderr)
    check_call(cmd, cwd=outdir, stdout=sys.stdout, stderr=sys.stderr)

    # copy files
    results_dir = os.path.dirname(clustered_out)
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    shutil.copy2(os.path.join(outdir, method + '.ASVs.fa'), clustered_out)
    with open(os.path.join(outdir, method + '.otu_table.txt'), 'rb') as i:
        with gzip.open(otutab_out, 'wb') as o:
            shutil.copyfileobj(i, o)


with file_logging(snakemake.log[0]):
    cluster_paired(
        snakemake.params.method,
        snakemake.input.demux,
        clustered_out=snakemake.output.clustered,
        otutab_out=snakemake.output.tab,
        unoise_program=snakemake.params.unoise_program,
        usearch_bin=snakemake.params.usearch_bin,
        usearch_par=snakemake.params.usearch_par,
        dada2_par=snakemake.params.dada2_par,
        threads=snakemake.threads,
    )
