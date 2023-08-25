import sys
from collections import OrderedDict, defaultdict
import traceback
from contextlib import contextmanager
import yaml
from typing import *


@contextmanager
def file_logging(f):
    with open(f, "w") as handle:
        sys.stderr = sys.stdout = handle
        try:
            yield handle
        except Exception:
            traceback.print_exc(file=handle)
            raise
        finally:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__


class OrderedDumper(yaml.SafeDumper):
    """
    OrderedDict representation in YAML
    Even though dicts have ordered keys since Python 3.7, using this makes
    sure that data is always written in the input order.
    """
    def __init__(self, *args, **kwargs):
        super(OrderedDumper, self).__init__(*args, **kwargs)
        r = lambda self, data:  self.represent_mapping(
            yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, data.items()
        )
        self.add_representer(OrderedDict, r)


def read_sample_file(path) -> Dict[
        str,  # layout: single / paired
        Dict[
            str,    # sample name
            List[  # list of paths with same sample name (may be > 1 with duplicate names and pool_duplicates=True)
                Tuple[str]  # read file paths: first = forward; [second = reverse]
            ]
        ]
    ]:
    """
    Read samples / paths to obtain back the same data structure as
    in group_samples from collect_sample_files.py
    """
    header_expected = ["sample", "R1", "R2"]    
    with open(path) as f:
        rdr = csv.reader(f, delimiter="\t")

        # read and validate header
        header = [field.strip().lower().replace(" ", "_") for field in next(rdr)]
        assert header == header_expected, \
            "Unexpected sample file header, expected columns: {}".format(
                ", ".join(header_expected))

        # read samples
        by_layout = defaultdict(lambda: defaultdict(list))
        for row in rdr:
            row = [entry.strip() for entry in row]
            if not row:
                continue
            sample = row[0]
            assert len(sample) > 0, "Empty sample for entry: {}".format(", ".join(row))
            assert len(row) == len(header), "Malformed sample file: unequal column numbers"
            paths = row[1:]
            if len(paths[1]) == 0:
                del paths[1]
                layout = "single"
            else:
                layout = "paired"
            assert len(paths[0]) > 0, "R1 path empty for entry: {}".format(", ".join(row))
            by_layout[layout][sample].append(paths)

    n_files = sum(len(paths) for paths in by_layout.values())
    assert n_files > 0, "Sample file was empty!"
    return by_layout
