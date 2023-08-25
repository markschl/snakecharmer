from hashlib import sha256
from itertools import product
import os
from copy import deepcopy
from os.path import dirname
from typing import *

from snakemake.workflow import srcdir

# Set environment variable of workflow root dir to allow post-deploy
# scripts to run other scripts stored in that directory
os.environ['PIPELINE_DIR'] = dirname(dirname(dirname(srcdir('.'))))


#### Setup ####

class Config(object):
    # strategy_names = ['merged', 'notmerged_R1', 'notmerged_R2']
    # working_dir = 'processing'
    # read_num_map = {'single': [1], 'paired': [1, 2]}
    
    # The following database types can be used directly
    # (specifically formatted or pre-trained), but 
    # cannot be imported to the internal taxonomy FASTA file format and then
    # filtered and converted to other formats.
    formatted_dbs = {'qiime_nb', 'idtaxa'}    

    def __init__(self, config):
        self.config = config
        self.workflow = config['workflows']
        self._get_primer_combinations()
        self._get_cmp_files()
        self._init_taxonomy()
        self._init_config()
        self._assemble_taxonomy()
        # from pprint import pprint; pprint(vars(self))

    def _get_cmp_files(self):
        # reads files that may be compared with ASVs/OTUS
        if 'compare' in self.config:
            self.cmp_files = deepcopy(self.config['compare'])
            cfg = self.cmp_files.pop('default_settings', {
                # TODO: investigate, how default settings can be taken from config.schema.yaml
                'maxaccepts': 64,
                'maxrejects': 64,
                'maxhits': 1,
            })
            for d in self.cmp_files.values():
                d.update(cfg)
        else:
            self.cmp_files = {}
    
    # def _init_dirs(self):
    #     # c_dir because the primers are placed there
    #     if not os.path.isdir(self.working_dir):
    #         os.makedirs(self.working_dir)

    def _get_primer_combinations(self):
        """
        Assemble the primer settings
        """
        # if self.config['primers']['process_unmerged'] is True:
        #     self.unmerged_read_idx = self.read_idx
        self.primers = {}
        self.primer_combinations = {}
        self.primer_combinations_flat = []
        self.markers = list(self.config['primers'])
        for marker, primers in self.config['primers'].items():
            if marker == 'trim_settings':
                # ignore settings
                continue
            assert isinstance(primers, dict), \
                'Invalid primer settings for marker {}'.format(marker)
            # parse primers
            pr = self.primers[marker] = {dir_: dict(parse_primers(
                primers[dir_])) for dir_ in ['forward', 'reverse']}
            # obtain primer combinations
            combinations = primers.get('combinations', 'default')
            if combinations == 'default':
                self.primer_combinations[marker] = []
                for fwd, rev in product(pr['forward'], pr['reverse']):
                    self.primer_combinations[marker].append(
                        '{}...{}'.format(fwd, rev))
                    self.primer_combinations_flat.append(
                        '{}__{}...{}'.format(marker, fwd, rev))
            else:
                assert isinstance(combinations, list), \
                    'Primer combinations of marker {} are not in list form'.format(marker)
                self.primer_combinations[marker] = combinations
                for c in combinations:
                    s = c.split('...')
                    assert (len(s) == 2), \
                        "Primer combinations must be in the form 'forward...reverse'. " \
                        "Encountered '{}' instead.".format(c)
                    assert (s[0] in self.primers['forward']), \
                        'Unknown forward primer: {}'.format(s[0])
                    assert (s[1] in self.primers['reverse']), \
                        'Unknown reverse primer: {}'.format(s[1])
                    self.primer_combinations_flat.append('{}__{}'.format(marker, c))
        # make sure the same primer combinations don't occur in different markers
        comb = [c for _, comb in self.primer_combinations.items()
                for c in comb]
        assert (len(set(comb)) == len(comb)), \
            'Identical primer combinations found across markers.' \
            'This is not allowed.'

    def _init_taxonomy(self):
        """"
        Sets up the taxonomy database sources and filtering options
        """
        # (1) parse contents of taxonomy.yaml:
        
        # Add "ids" (SHA-256 hashes of the database configuration),
        # which are used for caching database downloads.
        # We don't check for invalid/superfluous entries even though
        # they will alter the hash.
        self.taxdb_sources_by_hash = {}
        for name, dbconfig in self.config['taxonomy_db_sources'].items():
            _id = config_hash(dbconfig.items())
            dbconfig['source_id'] = _id
            dbconfig['name'] = name
            dbconfig['preformatted'] = 'preformatted' if dbconfig['format'] in self.formatted_dbs else 'regular'
            self.taxdb_sources_by_hash[_id] = dbconfig
        
        # (2) parse contents of config.yaml

        # Prepare taxonomy databases in config file, adding in the database
        # source configuration (such as type, URL, etc.) to have one dict
        # with all information.
        # We also hash the filtering options, the hashes are part of the database
        # file paths to allow caching of filtered+trained databases.
        # Snakemake uses nested paths for the caching:
        # db_hash / flt_hash / dbname.fasta.zst.
        self.taxdb_filter_by_hash = {}
        for dbs in self.config['taxonomy_dbs'].values():  # per-marker
            for name, db in dbs.items():  # dbs within marker
                _id = config_hash(((k, v) for k, v in db.items() if k != 'db'), empty_str="unchanged")
                db['filter_id'] = _id
                db['name'] = name
                self.taxdb_filter_by_hash[_id] = db
                db['source'] = self.config['taxonomy_db_sources'][db['db']]

    def _init_config(self):
        """"
        Initialize workflow settings, overriding defaults by workflow-specific
        configuration if present.
        """
        # parse workflow definitions
        self.workflows = self.config['workflows']
        del self.config['workflows']
        for name, p in self.workflows.items():
            p['name'] = name
            # clustering workflow
            # copy settings over, add extra settings overriding the defaults
            if 'settings' in p:
                settings = p['settings']
                p['settings'] = deepcopy(self.config)
                recursive_update(p['settings'], settings)
            else:
                p['settings'] = deepcopy(self.config)

    def _assemble_taxonomy(self):
        """
        Assemble the workflow-specific taxonomic assignment settings
        """
        # parse workflow definitions
        for name, p in self.workflows.items():
            # prepare list of taxonomy assignment methods (with corresponding databases)
            # first, replace 'default' by all db/method combinatinos
            settings = p['settings']
            if p['taxonomy'] == 'default':
                p['taxonomy'] = {
                    'dbs': {marker: list(dbs) for marker, dbs in settings['taxonomy_dbs'].items()},
                    'methods': list(settings['taxonomy_methods'])
                }
            else:
                assert isinstance(p['taxonomy'], dict)
                assert 'dbs' in p['taxonomy'], \
                    "'dbs' option missing in 'taxonomy' definition of workflow {}.".format(name)
                assert 'methods' in p['taxonomy'], \
                    "'methods' option missing in 'taxonomy' definition of workflow {}.".format(name)
            # validate method names
            method_names = p['taxonomy']['methods']
            method_cfg = settings['taxonomy_methods']
            strange_methods = set(method_names).difference(method_cfg)
            assert not strange_methods, \
                "Workflow {} has taxonomy method names that are not listed in 'taxonomy_methods'" \
                "below in the config file: {}".format(name, ', '.join(strange_methods))
            # set VSEARCH as default program if not already present
            # TODO: Jsonschema should do this, but the default keyword has no effect
            for cfg in method_cfg.values():
                if 'program' not in cfg:
                    cfg['program'] = 'vsearch'
            # then, for every marker, combine the assignment method and the
            # source database configurations into a single dictionary
            tax = {}
            for marker, db_names in p['taxonomy']['dbs'].items():
                db_cfg = settings['taxonomy_dbs'][marker]
                assert not set(db_names).difference(db_cfg)
                tax[marker] = {
                    (db, method): {
                        'marker': marker,
                        'db_name': db,
                        'method_name': method,
                        'assign': method_cfg[method],
                        **db_cfg[db]
                    }
                    for db, method in product(db_names, method_names)
                }
            # finally, replace the initial taxonomy configuration by the parsed one
            p['taxonomy'] = tax


    # allow access to workflow settings in a dict-like way
    def __getitem__(self, workflow):
        return self.workflows[workflow]


#### Helpers ####


def parse_primers(primers: Iterable[Dict[str, str]]) -> Generator[Tuple[str, List[str]], None, None]:
    for p in primers:
        assert isinstance(p, dict) and len(p) > 0, \
            "Primers must be defined in the form: 'name: sequence1,sequence2,...'"
        yield next(iter(p.items()))


def recursive_update(target, other):
    """
    Helper that recursively updates a nested dictionary with the contents of
    another dictionary
    """
    if hasattr(target, 'items') and hasattr(target, '__getitem__') and hasattr(other, 'items'):
        for name, values in other.items():
            target[name] = recursive_update(target[name], values) \
                if name in target else values
        return target
    return other


def config_hash(items: Iterable[Tuple[str, str]], empty_str: Optional[str] = None) -> str:
    """
    Simple function that takes an iterable of key/value pairs,
    which are joined, sorted and merged into a string, for which
    the SHA-256 digest is returned.
    """
    sorted_cfg = sorted(str(k) + str(v) for k, v in items)
    if empty_str is not None and sorted_cfg == "":
        return empty_str
    s = "".join(sorted_cfg).encode('utf-8')
    return sha256(s).hexdigest()


# def get_nested(d, *keys):
#     """
#     Get a nested dict entry or None if non-existent
#     """
#     if len(keys) == 0:
#         return d
#     try:
#         sub = d[keys[0]]
#         return get_nested(sub, *keys[1:])
#     except KeyError:
#         return None
