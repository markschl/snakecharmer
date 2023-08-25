#!/usr/bin/env python3

import csv
from glob import glob
import os
import re
import sys
from collections import defaultdict, OrderedDict
from typing import *


sample_patterns = {
    'illumina': r"(.+?)_S\d+_L\d+_R([12])_\d{3}\.fastq\.gz",
    'sample_read': r"(.+?)_R([12])\.fastq\.gz"
}


class SamplePatternMismatch(Exception):
    def __init__(self, sample, pattern, *args, **kwargs):
        msg = (
            'Sample name "{}" not matched by the given pattern "{}". Is '
            'the sample_pattern option correctly specified? '
            'Note that if specifying a directories list as input,'
            'all files need to be actual read files (e.g. Illumina index files present)'
            'in the directory as well will cause this error.\n'
            'Patterns provided by name are:\n{}\n'
            'Regex patterns can be tested e.g. on https://regex101.com or https://regexr.com'
        ).format(
            sample, pattern,
            "\n".join(f'* {name} = r"{patt}"' for name,
                      patt in sample_patterns.items())
        )
        super().__init__(msg, *args, **kwargs)


def format_list(l, cutoff=10):
    out = ", ".join(l[:10])
    if len(l) > cutoff:
        out += '...'
    return out


class RunWildcard(object):
    """
    Expands {run} wildcards in paths
    """

    def __init__(self, path, from_end=False):
        parts = list(path.split(os.sep))
        n_parts = len(parts)
        patterns = [(i, re.compile(re.escape(p).replace(
            "\\{run\\}", '(.+?)'))) for i, p in enumerate(parts) if "{run}" in p]
        assert len(
            patterns) == 1, "Invalid number of {run} wildcards (only one allowed)"
        self.pattern_index, self.pattern = patterns[0]
        if from_end:
            self.pattern_index = self.pattern_index - n_parts
        parts[self.pattern_index] = parts[self.pattern_index].replace(
            "{run}", "*")
        self.glob_pattern = os.path.join(*parts)

    def get_run(self, path):
        parts = path.split(os.sep)
        return self.pattern.fullmatch(parts[self.pattern_index]).group(1)


__run = "{run}"


def expand_file_pattern(pattern: str, recursive: bool = False) -> Generator[str, None, None]:
    # parse run wildcard
    run_wildcard = None
    n = pattern.count(__run)
    assert n <= 1, "Too many {run} wildcards (only one allowed)"
    if n == 1:
        if recursive:
            parts = pattern.split("**")
            if __run in parts[0]:
                run_wildcard = RunWildcard(parts[0])
                parts[0] = run_wildcard.glob_pattern
            else:
                assert __run in parts[-1], \
                    "The {run} wildcard cannot be in the middle between two recursive /**/ expressions, it must be in the first or last part of the pattern."
                run_wildcard = RunWildcard(parts[-1], from_end=True)
                parts[-1] = run_wildcard.glob_pattern
            pattern = "**".join(parts)
        else:
            run_wildcard = RunWildcard(pattern)
            pattern = run_wildcard.glob_pattern

    # do the glob search
    for path in glob(pattern, recursive=recursive):
        run = None if run_wildcard is None else run_wildcard.get_run(path)
        yield run, path


def list_dir(directory: str, recursive: bool = False) -> Generator[str, None, None]:
    # expand run directory
    if __run in directory:
        w = RunWildcard(directory)
        dirs = [(w.get_run(d), d)
                for d in glob(w.glob_pattern) if os.path.isdir(d)]
    else:
        dirs = [(None, directory)]

    for run, _dir in dirs:
        if recursive is True:
            for root, _, fnames in os.walk(_dir):
                for path in fnames:
                    yield run, os.path.join(root, path)
        else:
            for path in os.listdir(_dir):
                yield run, os.path.join(_dir, path)


def collect_files(
        directories: Optional[Iterable[str]] = None,
        patterns: Optional[Iterable[str]] = None,
        recursive: bool = False,
        sample_ext: Optional[Container[str]] = None,
) -> Generator[str, None, None]:
    """
    Collects all input files, given a sequence of 
    directories and/or a sequence of glob patterns.
    Relative paths will be interpreted relative to
    `base_dir`.
    The generator yields the individual file paths.
    """
    if directories is None and patterns is None:
        raise Exception(
            'At least one of "directories" and "patterns" must be defined in "input"')
    if patterns is not None:
        for pattern in patterns:
            found = False
            for run, path in expand_file_pattern(pattern, recursive):
                if os.path.isfile(path):
                    yield run, path
                    found = True
            assert found, "Pattern had no matches: '{}'".format(pattern)

    if directories is not None:
        for _dir in directories:
            found = False
            for run, f in list_dir(_dir, recursive):
                if os.path.isfile(f) and (sample_ext is None or any(f.endswith(ext) for ext in sample_ext)):
                    found = True
                    yield run, f
            assert found, (
                "No file found in {}. Maybe the allowed file extensions are "
                "too restrictive, or you meant to enable recursive search?"
                .format(_dir))


def get_sample_pattern(sample_pattern: str) -> str:
    """
    Prepares name patterns, returning the corresponding Regex
    for often used keywords.
    """
    try:
        sample_pattern = sample_patterns[sample_pattern.lower()]
    except KeyError:
        pass
    return sample_pattern


def parse_sample_pattern(f, pattern, reserved_pattern=None) -> Tuple[str, int, bool]:
    """
    Matches sample name and read number in a file name, given a
    name pattern.
    """
    m = pattern.fullmatch(f)
    if m is None:
        raise SamplePatternMismatch(f, pattern.pattern)
    try:
        sample_name = m.group(1)
    except IndexError:
        sample_name = m.group("sample")
    try:
        read = m.group(2)
    except IndexError:
        try:
            read = m.group("read")
        except AttributeError:
            # assuming no read group
            read = '1'
    assert (sample_name is not None and read is not None), (
        "Regular expression in 'sample_pattern' needs at least one group "
        "matching the sample name, and an optional second group matching "
        "the read number. They can also be named: (?P<sample>...) and (?P<read>...)")
    assert (read in ('1', '2')), \
        'Read number in file name must be 1 or 2, found instead "{}". ' \
        'Is the Regex pattern (sample_pattern) correct?'.format(read)
    # rename if necessary
    renamed = False
    if reserved_pattern is not None:
        sample_name, n = reserved_pattern.subn("_", sample_name)
        renamed = n > 0
    return sample_name, int(read), renamed


def collect_samples(sample_pattern: str, reserved_chars=None, **param) -> Tuple[str, Optional[str], str, int]:
    """
    This function collects sample files from `directories` and `patterns`
    (see `_collect_files`) and parses the sample names given using
    a defined pattern (see `parse_pattern`), optionally normalizing problematic
    characters.
    The generator yields a tuple of (sample name,  run name, file name, read number),
    whereby read number is 1 or 2.
    """
    _name_pat = re.compile(get_sample_pattern(sample_pattern))
    _reserved_pat = None if reserved_chars is None else re.compile(
        "[{}]".format(reserved_chars))
    renamed_samples = []
    for run, f in collect_files(**param):
        sample_name, read, renamed = parse_sample_pattern(
            os.path.basename(f), _name_pat, reserved_pattern=_reserved_pat)
        if renamed:
            renamed_samples.append(sample_name)
        yield sample_name, run, f, read
    if renamed_samples:
        renamed_samples = list(OrderedDict(
            zip(renamed_samples, renamed_samples)).keys())
        print(("Problematic characters in sample names were replaced by "
               "underscores: {}").format(format_list(renamed_samples)),
              file=sys.stderr)


def group_samples(
        paired_filter=None,
        **collect_args
) -> Dict[
    Optional[str],  # run name or None
    Dict[
        str,  # layout: single / paired
        Dict[
            str,        # sample name
            Tuple[str]  # read file paths: first = forward; [second = reverse]
        ]
    ]
]:
    # group files by run, then by sample name
    # TODO: defaultdict is only sorted in recent Python versions
    by_run = defaultdict(lambda: defaultdict(set))
    for sample_name, run, filename, read in collect_samples(**collect_args):
        by_run[run][sample_name].add((read, filename))

    # Now we can determine the read layout
    # Group by run -> read layout -> sample
    out = OrderedDict()
    r2_only = []
    for run, by_sample in by_run.items():
        by_layout = out[run] = OrderedDict(
            (('single', OrderedDict()), ('paired', OrderedDict())))
        for sample_name, read_files in by_sample.items():
            # sort the unique files (which were stored in a set)
            read_files = sorted(read_files, key=lambda x: x[0])
            # split into read indices and read paths
            read_idx, paths = zip(*read_files)
            # check for duplicates
            assert sorted(set(read_idx)) == sorted(read_idx), (
                "Several files with the same sample name found, are they from "
                "different runs? In this case, "
                "use the {{run}} wildcard.\n{}").format(format_list(paths))
            # now we can obtain the read layout (single/paired)s
            n_reads = len(read_idx)
            # other than 1 or 2 should be impossible (error caught earlier)
            assert n_reads in (1, 2)
            if n_reads == 2:
                # paired-end
                assert read_idx == (1, 2), (
                    "Two read files present for sample {}, but read numbers "
                    "are {} and {} instead of 1 and 2").format(sample_name, *read_idx)
                if paired_filter is None:
                    layout = 'paired'
                else:
                    if layout == "fwd-only":
                        layout = 'single'
                        paths = paths[0]
                    else:
                        layout = 'single_rev'
                        assert layout == "rev-only", \
                             f"Invalid value for 'paired_filter': {paired_filter}"
                        paths = paths[1]
            elif n_reads == 1:
                # single-end
                # other than 1 or 2 should be impossible (error caught earlier)
                assert read_idx[0] in (1, 2)
                if read_idx[0] == 2:
                    r2_only.append(sample_name)
                layout = 'single'
            by_layout[layout][sample_name] = paths
    if r2_only:
        print(("Only reverse read file (No. 2) present for some samples, "
               "is this correct?\n{}").format(format_list(r2_only)),
              file=sys.stderr
              )
    return out


def make_manifest(outdir, 
                  header_single, header_paired,
                  path_template="{layout}_{run}.tsv",
                  default_run="run1", 
                  relative_paths=True, 
                  **param):
    by_run = group_samples(**param)
    runs = []
    header_map = {"single": header_single, "paired": header_paired}
    unique_paths = set()  # for detecting duplicates
    for run, by_layout in by_run.items():
        if run is None:
            run = default_run
        for layout, samples in by_layout.items():
            if samples:
                outfile = path_template \
                    .format(outdir=outdir, layout=layout, run=run) \
                    .replace('/', os.sep)
                assert not outfile in unique_paths, (
                    "make_manifest: file {} already exists, is 'path_template' "
                    "correct?".format(outfile))
                unique_paths.add(outfile)
                if not os.path.exists(os.path.dirname(outfile)):
                    os.makedirs(os.path.dirname(outfile))
                with open(outfile, "w") as out:
                    w = csv.writer(out, delimiter="\t")
                    w.writerow(header_map[layout])
                    for sample, paths in samples.items():
                        if not relative_paths:
                            paths = [os.path.abspath(p) for p in paths]
                        w.writerow([sample] + list(paths))
                runs.append((run, layout))
    return runs


if __name__ == "__main__":
    import argparse
    from functools import partial
    comma_list = partial(str.split, sep=',')

    p = argparse.ArgumentParser(
        description="Script for making a manifest file as input for amplicon pipelines such as QIIME",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("outdir",
                   help="Output directory for manifest files."
                        "The files will be placed in the following path: "
                        "{outdir}/single/{run}.tsv or {outdir}/paired/{run}.tsv")
    p.add_argument("-d", "--directory", dest="directories", action="append",
                   help="Directory in which to search for FASTQ files. "
                        "Multiple patterns can be added by using this argument "
                        "several times")
    p.add_argument("-p", "--pattern", dest="patterns", action="append",
                   help="Glob pattern for finding FASTQ files. "
                        "Multiple patterns can be added by using this argument "
                        "several times")
    p.add_argument("-r", "--recursive", action="store_true",
                   help="Search directories recursively and match glob patterns "
                       "recursively (unlimited directory depth, specified using /**/)")
    p.add_argument("--rel", "--relative-paths", dest="relative_paths",
                   action="store_true",
                   help="Return paths relative to the current directory instead "
                        "of absolute paths (which is the default).")
    p.add_argument("--reserved", dest="reserved_chars", default="-. ",
                   help="List of reserved characters in sample names, "
                        "which will be converted to underscores. Default: '-. '")
    p.add_argument("-e", "--sample-ext", type=list,
                   default=[".fastq.gz", ".fq.gz"],
                   help="File extension(s) to match in directories")
    p.add_argument("-s", "--sample-pattern", default="illumina",
                   help="Sample pattern: a regular expression matching the whole"
                        "sample file name. It needs at least one group matching the "
                        "sample name, and an optional second group matching "
                        "the read number. They can also be named: (?P<sample>...) "
                        "and (?P<read>...)")
    p.add_argument("--hs", "--header-single", dest="header_single", type=comma_list,
                   default="sample-id,absolute-filepath",
                   help="Single-end manifest header (comma delimited list). "
                       "The default is the QIIME naming.")
    p.add_argument("--hd", "--header-paired", dest="header_paired", type=comma_list,
                   default="sample-id,forward-absolute-filepath,reverse-absolute-filepath",
                   help="Single-end manifest header (comma delimited list). "
                       "The default is the QIIME naming.")
    p.add_argument("--path-template", default="{outdir}/{layout}_{run}.tsv",
                   help="Template for creating the sample file(s). "
                        "Subdirectories are automatically created.")
    p.add_argument("-f", "--paired-filter", choices={None, "fwd-only", "rev-only"},
                   help="Keep only forward/R1 (-f fwd-only) or reverse/R2 "
                        "(-f rev-only) read files in a paired layout, resulting "
                        "in a single-end layout file. This setting does not affect "
                        "samples that are already in single-end mode. "
                        "If only keeping reverse/R2 reads, the resulting layout name will "
                        "actually be 'single_rev', indicating that the reads do contain "
                        "the reverse primer and not the forward primer")

    args = vars(p.parse_args())
    try:
        assert args["directories"] is not None or args["patterns"] is not None, \
            "Please specify at least one -d/--directory or -p/--pattern"

        make_manifest(**args)

    except Exception as e:
        print(f"Error: {e}\n\n", file=sys.stderr)
        p.print_help()
        exit(1)
