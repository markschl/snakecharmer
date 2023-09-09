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
    
    command = "cluster" if method == "uparse" else method

    cmd = [
        "amptk", command,
        "-i", os.path.basename(trimmed_in),
        "-o", method,
        "--cpus", str(threads)
    ]

    # add method specific arguments
    # TODO: inconsistent to have these settings under usearch even though used for DADA2 as well
    # (on the other hand, Amptk uses an USEACH-like procedure for DADA2 as well)
    maxee = usearch_par["merge"]["expected_length"] * usearch_par["filter"]["max_error_rate"]
    if method in ("unoise3", "uparse"):
        cmd += [
                "--usearch",
                usearch_bin,
                "--maxee",
                str(maxee),
            ]
        if method == "unoise3":
            cmd += [
                "--method",
                unoise_program,
                "--minsize",
                str(usearch_par["unoise3"]["min_size"]),
            ]
            if unoise_program == "usearch":
                assert usearch_bin is not None
                cmd += ["--usearch", usearch_bin]
            cluster_file = method + '.ASVs.fa'
        elif method == "uparse":
            assert usearch_bin is not None
            cmd += [
                "--minsize",
                str(usearch_par["uparse"]["min_size"]),
                "--usearch", usearch_bin,
            ]
            cluster_file = method + '.cluster.otus.fa'
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
        cluster_file = method + '.ASVs.fa'

    outdir = os.path.dirname(trimmed_in)
    print("Call: " + " ".join(cmd), file=sys.stderr)
    check_call(cmd, cwd=outdir, stdout=sys.stdout, stderr=sys.stderr)

    # copy files
    results_dir = os.path.dirname(clustered_out)
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    shutil.copy2(os.path.join(outdir, cluster_file), clustered_out)
    with open(os.path.join(outdir, method + '.otu_table.txt'), 'rb') as i:
        with gzip.open(otutab_out, 'wb') as o:
            shutil.copyfileobj(i, o)


with file_logging(snakemake.log[0]):
    cluster_paired(
        trimmed_in=snakemake.input.demux,
        clustered_out=snakemake.output.clustered,
        otutab_out=snakemake.output.tab,
        threads=snakemake.threads,
        **snakemake.params
    )
