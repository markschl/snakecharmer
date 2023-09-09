from collections import defaultdict
from hashlib import sha256
from itertools import product
import copy
from typing import *


from scripts.utils.sample_list import SampleList


def list_runs(input_cfg):
    for technology, data in input_cfg.items():
        for run, data in data.items():
            if isinstance(data, str):
                sample_file = data
                reverse = False
            else:
                sample_file = data["samples"]
                reverse = data.get("reverse", False)
            yield {
                "technology": technology,
                "run": run,
                "sample_file": sample_file,
                "orientation": "reverse" if reverse else "forward"
            }


class Config(object):
    # strategy_names = ['merged', 'notmerged_R1', 'notmerged_R2']
    # working_dir = 'workdir'
    # read_num_map = {'single': [1], 'paired': [1, 2]}
    
    # The following database types can be used directly
    # (specifically formatted or pre-trained), but 
    # cannot be imported to the internal taxonomy FASTA file format and then
    # filtered and converted to other formats.
    formatted_dbs = {'qiime_nb', 'idtaxa'}    

    def __init__(self, config):
        self.config = config
        self.workflow = config['workflows']
        self._init_input()
        self._get_primer_combinations()
        self._init_taxonomy()
        self._init_config()
        self._assemble_taxonomy()
        # from pprint import pprint; pprint(vars(self))

    def _init_input(self):
        grouped = defaultdict(lambda: defaultdict(dict))  # technology -> layout -> metadata
        self.runs = self.run_pools = OrderedDict()  # (run, layout) -> metadata
        for d in list_runs(self.config["input"]):
            run = d["run"]
            # read sample list and add some more metadata
            assert not run in self.runs, (
                f"Duplicate run name found: {run}. This is not allowed, "
                "even acros different sequencing technologies."
            )
            d["layout"] = layout = SampleList.infer_layout(d["sample_file"])
            if layout == "single" and d["orientation"] == "reverse":
                d["layout"] = layout = "single_rev"
            grouped[d["technology"]][layout][run] = d
            self.runs[(run, layout)] = d

        # the 'run_pools' dict should contain pooled runs if pool_raw: true
        if self.config["pool_raw"] is True:
            self.run_pools = OrderedDict()
            for technology, data in sorted(grouped.items()):
                for layout, runs in sorted(data.items()):
                    run_concat = "_".join(sorted(runs))
                    run = run_concat + "_pool"
                    self.run_pools[(run, layout)] = {
                        "technology": technology,
                        "run_list": run_concat,
                        "run": run,
                        "layout": layout,
                        "sample_files": sorted(r["sample_file"] for r in runs.values())
                    }

    def _get_primer_combinations(self):
        """
        Assemble the primer settings
        """
        # if self.config['primers']['process_unmerged'] is True:
        #     self.unmerged_read_idx = self.read_idx
        self.primers = {}
        self.primer_combinations = {}
        self.primer_combinations_flat = []
        self.primer_combinations_nomarker = set()
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
                    comb = f'{fwd}...{rev}'
                    self.primer_combinations[marker].append(comb)
                    self.primer_combinations_nomarker.add(comb)
                    self.primer_combinations_flat.append(f'{marker}__{comb}')
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
                    self.primer_combinations_nomarker.add(c)
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
                p['settings'] = copy.deepcopy(self.config)
                recursive_update(p['settings'], settings)
            else:
                p['settings'] = copy.deepcopy(self.config)

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

    def _get_runs(self, workflow=None, pooled=True):
        if pooled and (workflow is None or self[workflow]["settings"]["pool_raw"]):
            return self.run_pools
        return self.runs

    def get_runs(self, workflow=None, pooled=True):
        for d in self._get_runs(workflow, pooled).values():
            is_pool = d["run"].endswith("_pool")
            # print(workflow, pooled, d, is_pool)
            if pooled or not is_pool:
                yield d

    def get_run_data(self, workflow=None, run=None, layout=None, pooled=True, **unused):
        return self._get_runs(workflow, pooled=pooled)[(run, layout)]

    def read_samples(self, path, workflow=None, run=None, layout=None, pooled=True, **unused):
        d = self.get_run_data(workflow, run, layout, pooled=pooled)
        path = path.format(workflow=workflow, **d)
        l = SampleList(path)
        return {"sample": [s for s, _ in l.samples()], "read": [str(i+1) for i in range(l.n_reads)]}

    def tax_config(self, workflow, marker, db_name, tax_method, **ignored):
        return self.workflows[workflow]["taxonomy"][marker][(db_name, tax_method)]


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
