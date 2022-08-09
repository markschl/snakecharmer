from collections import defaultdict
from copy import deepcopy
import glob
from itertools import product
import re
import sys
import os
import glob
from seq import FastaIO, seq
from seq.alignment import consensus

from snakemake.workflow import srcdir
from os.path import dirname

# Set environment variable of pipeline root dir to allow post-deploy 
# scripts to run other scripts stored in that directory
# TODO: may be kind of a hack, not sure if it works in all cases
os.environ['PIPELINE_DIR'] = dirname(dirname(dirname(srcdir('.'))))


#### Configuration ####


def collect_samples(name_pattern, *args, **kwargs):
    _name_pat = re.compile(parse_pattern(name_pattern))
    for f in collect_sample_files(*args, **kwargs):
        root, fname = os.path.split(os.path.abspath(f))
        sample_name, read = parse_sample(fname, _name_pat)
        yield sample_name, root, fname, read


def collect_sample_files(directories=None, patterns=None, recursive=False):
    if directories is None and patterns is None:
        raise Exception('At least one of "directories" and "patterns" must be defined in "input"')
    if patterns is not None:
        for pattern in patterns:
            for f in glob.glob(os.path.expanduser(pattern), recursive=recursive):
                yield f
    if directories is not None:
        if recursive is True:
            for rd in directories:
                for root, _, fnames in os.walk(rd):
                    for f in fnames:
                        yield os.path.join(root, f)
        else:
            for root in directories:
                for f in os.listdirs(root):
                    yield os.path.join(root, f)


def parse_pattern(name_pattern):
    name_pattern = name_pattern.lower()
    if name_pattern == 'illumina':
        name_pattern = r"(.+?)_S\d+_L\d+_R([12])_\d{3}\.fastq\.gz"
    elif name_pattern == 'sample_read':
        name_pattern = r"(.+?)_R([12])\.fastq\.gz"        
    return name_pattern


def parse_sample(f, pattern):
    m = pattern.match(f)
    if m is None:
        #return None, None
        raise Exception('Invalid sample name: {}'.format(f))
    sample_name = m.group(1)
    read = m.group(2)
    assert read in ('1', '2')
    if '-' in sample_name:
        print('- in sample name: {}. Does not work with USEARCH pipeline'.format(sample_name), file=sys.stderr)
    return sample_name, int(read)


def group_samples(**collect_args):
    # group by sample -> read
    d = defaultdict(lambda: defaultdict(set))
    for sample_name, root, fname, read in collect_samples(**collect_args):
        d[sample_name][read].add((root, fname))
    # group by forward only / paired end
    by_strategy = {'single': [], 'paired': []}
    for sample_name, by_read in d.items():
        by_read = sorted(by_read.items())
        n = len(by_read)
        assert n > 0
        if n > 2:
            raise Exception('More than 2 read files for sample {}'.format(sample_name))
        elif n == 2:
            if len(by_read[0][1]) != len(by_read[1][1]):
                raise Exception('Read file number mismatch for sample {}'.format(sample_name))
            by_strategy['paired'].append((sample_name, by_read))
        elif n == 1:
            by_strategy['single'].append((sample_name, by_read))
    return by_strategy



def parse_primers(primers):
    for p in primers:
        assert isinstance(p, dict)
        id = next(iter(p.keys()))
        seqs = next(iter(p.values()))
        seqs = [s.strip() for s in seqs.split(',')]
        yield id, seqs


def make_primer_fasta(primers, outfile):
    with open(outfile, 'w') as f:
        for id, seqs in primers.items():
            for s in seqs:
                FastaIO.write(seq.SeqRecord(id, s), f)


def recursive_update(target, other):
    if hasattr(target, 'items') and hasattr(target, '__getitem__') and hasattr(other, 'items'):
        for name, values in other.items():
            target[name] = recursive_update(target[name], values) if name in target else values
        return target
    return other


#### Setup ####

class Config(object):
    #strategy_names = ['merged', 'notmerged_R1', 'notmerged_R2']
    workding_dir = 'processing'
    read_num_map = {'single': [1], 'paired': [1, 2]}

    def __init__(self, config):
        self.config = config
        self.pipeline = config['pipelines']
        self._init_dirs()
        self._read_samples()
        self._read_primers()
        self._read_cmp_files()
        self._init_config()
        # from pprint import pprint; pprint(vars(self))
    
    def _read_samples(self):
        _cfg = self.config['input']
        self.samples = group_samples(
            directories=_cfg.get('directories', None),
            patterns=_cfg.get('patterns', None),
            name_pattern=_cfg['name_pattern'],
            recursive=_cfg.get('recursive', False)
        )
        self.sample_names = {strategy: [s for s, _ in samples] 
                            for strategy, samples in self.samples.items()}
        self.sequencing_strategies = [name for name, samples in self.samples.items() if samples]
        self.strategy_sample_read = [
            (strategy, sample, read)
            for strategy, samples in self.samples.items()
            for sample, _ in samples
            for read in self.read_num_map[strategy]
        ]

    def _read_cmp_files(self):
        # reads files that may be compared with ASVs/OTUS
        if 'compare' in self.config:
            self.cmp_files = deepcopy(self.config['compare'])
            cfg = self.cmp_files.pop('default_settings', {})
            for d in self.cmp_files.values():
                d.update(cfg)
        else:
            self.cmp_files = {}

    def _init_dirs(self):
        # c_dir because the primers are placed there
        if not os.path.isdir(self.workding_dir):
            os.makedirs(self.workding_dir)

    def _read_primers(self):
        # if self.config['primers']['process_unmerged'] is True:
        #     self.unmerged_read_idx = self.read_idx

        self.fwd_primers = dict(parse_primers(self.config['primers']['forward']))
        self.rev_primers = dict(parse_primers(self.config['primers']['reverse']))
        self.fwd_primers_rev = {p: [seq.reverse_complement(s) for s in seqs] for p, seqs in self.fwd_primers.items()}
        self.rev_primers_rev = {p: [seq.reverse_complement(s) for s in seqs] for p, seqs in self.rev_primers.items()}

        self.fwd_primers_consensus = {name: consensus(seq.SeqRecord(name, s) for s in primers)
                                      for name, primers in self.fwd_primers.items()}
        self.rev_primers_consensus = {name: consensus(seq.SeqRecord(name, s) for s in primers)
                                      for name, primers in self.rev_primers.items()}
        self.rev_primers_consensus_rev = {name: seq.reverse_complement(s) for name, s in self.rev_primers_consensus.items()}
        self.fwd_primers_consensus_rev = {name: seq.reverse_complement(s) for name, s in self.fwd_primers_consensus.items()}

        # obtain primer combinations
        self.primer_combinations = []
        comb = self.config['primers'].get('combinations', 'default')
        if comb == 'default':
            assert len(self.fwd_primers) == 1 and len(self.rev_primers) == 1, \
                "With multiple forward/reverse primers, a list of combinations needs to be specified (forward...reverse)"
            self.primer_combinations.append(list(self.fwd_primers)[0] + '...' + list(self.rev_primers)[0])
        else:
            assert isinstance(comb, list)
            for c in comb:
                s = c.split('...')
                assert len(s) == 2
                assert s[0] in self.fwd_primers, 'Unknown forward primer: {}'.format(s[0])
                assert s[1] in self.rev_primers, 'Unknown reverse primer: {}'.format(s[1])
                self.primer_combinations.append(c)

    def _init_config(self):
        # prepare taxonomy databases
        for db in self.config['taxonomy_dbs'].values():
            # complete optional values
            if not 'defined' in db:
                db['defined'] = '_all'
            db.update(self.config['taxonomy_db_sources'][db['db']])

        # parse pipeline definitions
        self.pipelines = self.config['pipelines']
        del self.config['pipelines']
        for name, p in self.pipelines.items():
            p['name'] = name
            # clustering pipeline
            # copy settings over, add extra settings overriding the defaults
            if 'settings' in p:
                settings = p['settings']
                p['settings'] = deepcopy(self.config)
                recursive_update(p['settings'], settings)
            else:
                p['settings'] = deepcopy(self.config)
            
            # determine/adjust sequencing strategies 
            # TODO: some more work needed to activate this option
            # if p.get('forward_only', None) is True:
            #     p['sequencing_strategies'] = ['single']
            # else:
            p['sequencing_strategies'] = self.sequencing_strategies

            # prepare list of taxonomy assignment methods (with corresponding databases)
            settings = p['settings']
            dbs = settings['taxonomy_dbs']
            db_names = list(dbs.keys())
            methods = settings['taxonomy_methods']
            method_names = list(methods.keys())
            if p['taxonomy'] != 'default':
                _dbs = p['taxonomy']['dbs']
                _methods = p['taxonomy']['methods']
                assert not set(_dbs).difference(db_names)
                assert not set(_methods).difference(method_names)
                db_names = _dbs
                method_names = _methods
            p['taxonomy'] = {
                (db, method): {'db_name': db, 'method_name': method, **dbs[db], **methods[method]}
                for db, method in product(db_names, method_names)
            }

            # is it a 'simple' pipeline (with only one results dir, see setup_project())?
            p['is_simple'] = len(p["sequencing_strategies"]) == 1 and len(self.primer_combinations) == 1
            # simplify path generation
            p['single_primercomb'] = self.primer_combinations[0] if p['is_simple'] else None
            p['single_strategy'] = p['sequencing_strategies'][0] if p['is_simple'] else None

    # allow access to pipeline settings in a dict-like way
    def __getitem__(self, pipeline):
        return self.pipelines[pipeline]


#### Helpers ####


from snakemake.io import Log
from contextlib import contextmanager
import traceback


@contextmanager
def file_logging(f):
    if isinstance(f, Log):
        f = f[0]
    assert isinstance(f, str)
    with open(f, "w") as handle:
        try:
            yield handle
        except Exception as e:
            traceback.print_exc(file=handle)
            raise

