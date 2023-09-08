import yaml

from utils import file_logging
from utils.sample_list import SampleList


def make_pooling_list(sample_files, sample_file_out, info_file_out):
    # read sample files
    # TODO: indexes.tsv for demultiplexing
    sample_lists = [SampleList(f) for f in sample_files]
    # there should be only one layout present (paired or single)
    unique_layout = set(l.layout for l in sample_lists)
    assert len(unique_layout) == 1
    layout = next(iter(unique_layout))
    # get a map of sample -> read files
    sample_dicts = [dict(l.samples()) for l in sample_lists]
    unique_names = sorted(set(s for d in sample_dicts for s in d))
    # write sample lists
    sample_list_out = SampleList(layout=layout)
    sample_pools = {}
    for sample in unique_names:
        # obtain a list of all run files for that sample
        run_files = [d[sample] for d in sample_dicts if sample in d]
        # sort the runs by path to have a consistent output
        run_files.sort()
        # group by R1/R2
        run_files = list(zip(*run_files))
        out_reads = [f"{sample}_R{i+1}.fastq.gz" for i in range(sample_list_out.n_reads)]
        sample_list_out.add(sample, out_reads)
        sample_pools[sample] = {f"R{i+1}": f for i, f in enumerate(run_files)}
    # write new sample list to files
    # TODO: we only list file names, not paths. This sample list is mostly needed
    # for the pipeline to work, but the paths are unimportant since already known
    with open(sample_file_out, "w") as f:
        sample_list_out.write(f)
    # write YAML file with all pooling info
    # note: sample keys are sorted (unless sort_keys=False), but the
    # but the order of read pooling is consistent (files are in a sorted list)
    with open(info_file_out, "w") as f:
        yaml.safe_dump(sample_pools, f)


with file_logging(snakemake.log[0]):
    make_pooling_list(sample_files=snakemake.input.sample_files,
                      sample_file_out=snakemake.output.sample_file,
                      info_file_out=snakemake.output.yml)
