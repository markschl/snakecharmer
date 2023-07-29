from os.path import relpath
import yaml
from collections import OrderedDict

from utils import file_logging, OrderedDumper



def dump_yaml(samples, outfile):
    # YAML file
    samples = OrderedDict((
        strategy, OrderedDict((
            sample, OrderedDict((
                "R{}".format(i + 1), [
                    relpath(path) for path in read_paths
                ]
                if len(read_paths) > 1
                else relpath(read_paths[0]))
                for i, read_paths in enumerate(zip(*sample_paths))
            ))
            for sample, sample_paths in strategy_paths.items()
        ))
        for strategy, strategy_paths in samples.items()
    )
    with open(outfile, "w") as o:
        yaml.dump(samples, o, Dumper=OrderedDumper)



with file_logging(snakemake.log[0]):
    dump_yaml(snakemake.params.samples, snakemake.output.yaml)
