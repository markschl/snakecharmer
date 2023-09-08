from os.path import abspath
from subprocess import check_call
import sys

from utils import file_logging


def run_uvsnake(snakefile, command, threads=1):
    cmd = [
        "snakemake",
        "--use-conda",
        "--cache",
        "--cores", str(threads),
        "--snakefile", abspath(snakefile),
        f"uvsnake_{command}"
    ]
    print("Running {}...".format(" ".join(cmd)), file=sys.stderr)
    check_call(cmd, stdout=sys.stdout, stderr=sys.stderr)


with file_logging(snakemake.log[0]):
    run_uvsnake(
        snakefile=snakemake.input.snakefile,
        command=snakemake.params.command,
        threads=snakemake.threads,
    )
