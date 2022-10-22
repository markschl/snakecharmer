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
# TODO: kind of a hack, not sure if it works in all cases; possibly replace with standard Snakemake scripts
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
    pat = name_pattern.lower()
    if pat == 'illumina':
        name_pattern = r"(.+?)_S\d+_L\d+_R([12])_\d{3}\.fastq\.gz"
    elif pat == 'sample_read':
        name_pattern = r"(.+?)_R([12])\.fastq\.gz"
    return name_pattern


def parse_sample(f, pattern):
    m = pattern.match(f)
    if m is None:
        raise Exception('Sample name "{}" not matched by Regex pattern "{}". Is the name_pattern option correctly specified? Regex patterns can be debugged e.g. on https://regexr.com'.format(f, pattern.pattern))
    sample_name = m.group(1)
    read = m.group(2)
    assert sample_name is not None and read is not None, 'Regular expression in "name_pattern" needs to have exactly two groups, the first for the sample name and the second for the read number.'
    assert read in ('1', '2'), 'Read number in file name must be 1 or 2, found instead "{}". Is the Regex pattern (name_pattern) correct?'.format(read)
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
        assert isinstance(p, dict), "Primers must be defined in the form: 'name: sequence1,sequence2,...'"
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
        self.primers = {}
        self.primers_rev = {}
        self.primers_consensus = {}
        self.primers_consensus_rev = {}
        self.primer_combinations = {}
        self.primer_combinations_flat = []
        self.markers = list(self.config['primers'])
        for marker, primers in self.config['primers'].items():
            if marker == 'trim_settings':
                # ignore settings
                continue
            assert isinstance(primers, dict), 'Invalid primer settings for marker {}'.format(marker)
            # parse primers
            pr = self.primers[marker] = {dir_: dict(parse_primers(primers[dir_])) for dir_ in ['forward', 'reverse']}
            # reverse complement versions
            self.primers_rev[marker] = {
                direction: {p: [seq.reverse_complement(s) for s in seqs] for p, seqs in pseqs.items()}
                for direction, pseqs in pr.items()
            }
            # primer consensus
            cns = self.primers_consensus[marker] = {
                direction: {name: consensus(seq.SeqRecord(name, s) for s in seqs)
                                      for name, seqs in pseqs.items()}
                for direction, pseqs in pr.items()
            }
            # reverse complement consensus
            self.primers_consensus_rev[marker] = {
                direction: {name: seq.reverse_complement(s) for name, s in cons.items()}
                for direction, cons in cns.items()
            }
            # obtain primer combinations
            combinations = primers.get('combinations', 'default')
            if combinations == 'default':
                self.primer_combinations[marker] = []
                for fwd, rev in product(pr['forward'], pr['reverse']):
                    self.primer_combinations[marker].append('{}...{}'.format(fwd, rev))
                    self.primer_combinations_flat.append('{}__{}...{}'.format(marker, fwd, rev))
            else:
                assert isinstance(combinations, list)
                self.primer_combinations[marker] = combinations
                for c in combinations:
                    s = c.split('...')
                    assert len(s) == 2
                    assert s[0] in self.primers['forward'], 'Unknown forward primer: {}'.format(s[0])
                    assert s[1] in self.primers['reverse'], 'Unknown reverse primer: {}'.format(s[1])
                    self.primer_combinations_flat.append('{}__{}'.format(marker, c))
        # make sure the same primer combinations don't occur in different markers
        comb = [c for _, comb in self.primer_combinations.items() for c in comb]
        assert len(set(comb)) == len(comb), 'Primer combinations cannot be identical across markers'

    def _init_config(self):
        # prepare taxonomy databases
        for dbs in self.config['taxonomy_dbs'].values():
            for db in dbs.values():
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
            # first, replace 'default'
            settings = p['settings']
            if p['taxonomy'] == 'default':
                p['taxonomy'] = {
                    'dbs': {marker: list(dbs) for marker, dbs in settings['taxonomy_dbs'].items()},
                    'methods': list(settings['taxonomy_methods'])
                }
            else:
                assert 'dbs' in p['taxonomy']
                assert 'methods' in p['taxonomy']
            # validate method names
            method_names = p['taxonomy']['methods']
            method_cfg = settings['taxonomy_methods']
            assert not set(method_names).difference(method_cfg)
            # then, for every marker, assemble the combinations
            tax = {}
            for marker, db_names in p['taxonomy']['dbs'].items():
                db_cfg = settings['taxonomy_dbs'][marker]
                assert not set(db_names).difference(db_cfg)
                tax[marker] = {
                    (db, method): {'marker': marker, 'db_name': db, 'method_name': method, **db_cfg[db], **method_cfg[method]}
                    for db, method in product(db_names, method_names)
                }
            # finally, replace the initial taxonomy configuration by the parsed one
            p['taxonomy'] = tax

            # is it a 'simple' pipeline (with only one results dir, see setup_project())?
            p['is_simple'] = len(p["sequencing_strategies"]) == 1 and len(self.primer_combinations_flat) == 1
            # simplify path generation
            p['single_primercomb'] = self.primer_combinations_flat[0] if p['is_simple'] else None
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

