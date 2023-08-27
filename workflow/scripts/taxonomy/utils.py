from contextlib import contextmanager
import sys
import traceback


# TODO: code repeated from scripts/utils/utils.py
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
