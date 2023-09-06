from utils import file_logging
from utils.sample_list import SampleList


def collect_sample_lists(run_meta, reserved_chars, path_template):
    for d in run_meta:
        tsv_out = path_template.format(ext="tsv", **d)
        l = SampleList(d["sample_file"], reserved_chars=reserved_chars)
        with open(tsv_out, "w") as o:
            l.write(o)
        yaml_out = path_template.format(ext="yaml", **d)
        with open(yaml_out, "w") as o:
            l.write_yaml(o)
    

with file_logging(snakemake.log[0]):
    collect_sample_lists(**snakemake.params)
