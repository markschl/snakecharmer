import sys
from collections import OrderedDict
import traceback
from contextlib import contextmanager
import yaml


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
