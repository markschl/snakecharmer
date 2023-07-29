from utils import file_logging, FastaIO
import yaml


def _make_fasta(primers, outfile, add_anchor):
    with open(outfile, 'w') as f:
        for _id, d in primers:
            for s in d['seq']:
                if d['anchor']:
                    s = add_anchor(s)
                FastaIO.write(FastaIO.SeqRecord(_id, s), f)


def make_primer_fasta(primers, fwd_out, rev_out):
    _make_fasta(primers["forward"].items(), fwd_out, lambda s: '^' + s)
    _make_fasta(primers["reverse_rev"].items(), rev_out, lambda s: s + '$')


with file_logging(snakemake.log[0]):
    with open(snakemake.input.yaml) as f:
        primers = yaml.safe_load(f)
    make_primer_fasta(primers, snakemake.output.forward, snakemake.output.reverse_rev)
