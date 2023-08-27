import sys
import traceback
from contextlib import contextmanager


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
