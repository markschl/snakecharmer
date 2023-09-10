import re
import sys
from subprocess import check_call

from utils import file_logging


def denoise_paired(
        input, 
        layout,
        seq_out, 
        tab_out, 
        stats_out,
        trunc_qual=None,
        trunc_len=None,
        max_err=None,
        # the following settings have defaults (can be undefined without error)
        merge_maxdiffs=0,
        chimera_method=None,
        pooling_method=None,
        threads=1,
        **unused
):
    cmd = ["qiime", "dada2"]

    # https://docs.qiime2.org/2023.7/plugins/available/dada2/denoise-paired/
    if layout == "paired":
        cmd += [
            "denoise-paired",
            "--p-trunc-len-f", "{:.0f}".format(trunc_len["fwd"]),
            "--p-trunc-len-r", "{:.0f}".format(trunc_len["rev"]),
            "--p-max-ee-f", str(max_err["fwd"]),
            "--p-max-ee-r", str(max_err["rev"]),
        ]
    else:
        dir_ = "fwd"
        if layout.endswith(".rev"):
            layout = re.sub(r"\.rev$", "", layout)
            dir_ = "rev"
        assert layout == "single"
        # https://docs.qiime2.org/2023.7/plugins/available/dada2/denoise-single
        cmd += [
            "denoise-single",
            "--p-trunc-len", "{:.0f}".format(trunc_len[dir_]),
            "--p-max-ee", str(max_err[dir_]),
        ]

    # common options

    if chimera_method is None:
        chimera_method = "consensus"
    elif chimera_method == "per-sample":
        chimera_method = "none"
    if pooling_method is None:
        pooling_method = "independent"
    elif pooling_method == "pooled":
        print('Warning: "pooling_method = pooled" not possible for QIIME2',
              file=sys.stderr)
    cmd += [
        "--i-demultiplexed-seqs", input,
        "--p-trunc-q", "{:.0f}".format(trunc_qual),
        "--p-n-reads-learn", "1000000",  # TODO: not configurable
        "--p-chimera-method", chimera_method,
        "--verbose",
        "--p-n-threads", str(threads),
        "--o-representative-sequences", seq_out,
        "--o-table", tab_out,
        "--o-denoising-stats", stats_out
    ]

    if merge_maxdiffs > 0:
        if merge_maxdiffs > 1:
            print("WARNING: 'merge_maxdiffs' is > 1, but QIIME only allows up to one difference.",
                  file=sys.stderr)
        cmd.append("--p-allow-one-off")

    print("Call: " + " ".join(cmd), file=sys.stderr)
    check_call(cmd, stdout=sys.stdout, stderr=sys.stderr)


with file_logging(snakemake.log[0]):
    denoise_paired(
        snakemake.input.trim,
        layout=snakemake.wildcards.layout,
        seq_out=snakemake.output.asvs,
        tab_out=snakemake.output.tab,
        stats_out=snakemake.output.stats,
        threads=snakemake.threads,
        **snakemake.params.par
    )
