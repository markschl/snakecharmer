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

    def __init__(self, config):
        self.config = config
        self.workflow = config['workflows']
        self._init_input()
        self._get_primer_combinations()
        self._init_config()
        self._init_taxonomy()
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
                d["layout"] = layout = "single.rev"
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

    def _init_config(self):
        """"
        Initialize workflow settings, overriding defaults by workflow-specific
        configuration if present.
        """
        # parse workflow definitions
        self.workflows = self.config['workflows']
        reduced_cfg = copy.deepcopy(self.config)
        del reduced_cfg['workflows']
        del reduced_cfg['taxonomy']
        del reduced_cfg['taxonomy_db_sources']
        for name, p in self.workflows.items():
            p['name'] = name
            p['pipeline'] = p['cluster'].split('_', 1)[0]
            # clustering workflow
            # copy settings over, add extra settings overriding the defaults
            if 'settings' in p:
                settings = p['settings']
                p['settings'] = copy.deepcopy(reduced_cfg)
                recursive_update(p['settings'], settings)
            else:
                p['settings'] = copy.deepcopy(reduced_cfg)
        # possible technology/layout combinations that work
        # (populated in corresponding snakefiles)
        self.cluster_capabilities = {}

    def _init_taxonomy(self):
        """"
        Sets up the taxonomy database sources and filtering options
        """
        # first deep-copy the relevant settings, since empty dicts seem to be
        # shared and thus unwanted modifications can happen
        tax_cfg = self.config['taxonomy'] = copy.deepcopy(self.config['taxonomy'])
        tax_source_cfg = self.config['taxonomy_db_sources'] = copy.deepcopy(self.config['taxonomy_db_sources'])
        # First, parse contents of taxonomy.yaml:
        # Add "ids" (SHA-256 hashes of the database configuration),
        # which are used for caching database downloads.
        self.taxdb_sources_by_hash = {}
        for dbconfig in tax_source_cfg.values():
            d = copy.deepcopy(dbconfig)
            dbconfig['source_id'] = _id = config_hash(d.items())
            self.taxdb_sources_by_hash[_id] = d
        
        # Prepare taxonomy databases listed in config.yaml:
        # We hash the filtering options, the hashes are part of the database
        # file paths to allow caching of filtered+trained databases.
        # Snakemake uses nested paths for the caching:
        # db_hash / flt_hash / dbname.fasta.zst.
        self.taxdb_filter_by_hash = {}
        for marker_dbs in tax_cfg['dbs'].values():  # per-marker
            for flt_config in marker_dbs.values():  # dbs within marker
                clean_db = copy.deepcopy(flt_config)
                del clean_db['db']
                flt_config['filter_id'] = _id = config_hash(clean_db.items(), empty_str="unfiltered")
                if _id not in self.taxdb_filter_by_hash:
                    self.taxdb_filter_by_hash[_id] = clean_db

        # Similarly as above, get database training settings (if any),
        # referenced by their hash
        self.taxdb_training_cfg_by_hash = {}
        for mcfg in tax_cfg["methods"].values():
            tcfg = mcfg["train"] = copy.deepcopy(mcfg["train"])
            d = copy.deepcopy(tcfg)
            tcfg['conversion_id'] = _id = config_hash(d.items(), empty_str="standard")
            if _id not in self.taxdb_training_cfg_by_hash:
                self.taxdb_training_cfg_by_hash[_id] = d
            print("added", tcfg, d, self.taxdb_training_cfg_by_hash)

        # Assemble the workflow-specific taxonomic assignment settings
        # First, initialize global method/db combinations
        tax_combinations = tax_cfg.get('combinations', 'all')
        if tax_combinations == 'all':
            tax_combinations = {
                marker: list(product(db_names, tax_cfg['methods']))
                for marker, db_names in tax_cfg['dbs'].items()
            }
        # individual workflow definitions
        for wcfg in self.workflows.values():
            workflow_comb = wcfg.get('taxonomy', 'all')
            if workflow_comb == 'all':
                workflow_comb = tax_combinations
            tax = {}
            method_config = tax_cfg["methods"]
            for marker, comb in workflow_comb.items():
                marker_cfg = tax[marker] = {}
                marker_db_cfg = tax_cfg['dbs'][marker]
                for db_name, method_name in comb:
                    try:
                        filter_cfg = marker_db_cfg[db_name]
                    except KeyError:
                        raise Exception(f"Unknown database: {db_name}")
                    try:
                        d = method_config[method_name]
                    except KeyError:
                        raise Exception(f"Unknown taxonomy assignment method: {method_name}")
                    try:
                        source = tax_source_cfg[filter_cfg['db']]
                    except KeyError:
                        raise Exception(f"Unknown taxonomy source database: {filter_cfg['db']}")
                    marker_cfg[(db_name, method_name)] = {
                        'marker': marker,
                        'db_name': db_name,
                        'method_name': method_name,
                        # make sure to copy everything, so later modifications
                        # (outside of the Config object in common.smk) will not
                        # cause confusion
                        'assign': copy.deepcopy(d),
                        'source': copy.deepcopy(source),
                        'filter': copy.deepcopy(filter_cfg)
                    }
            wcfg["taxonomy"] = tax

        # The following dict contains format requirements for database assignment
        # methods, it is populated in individual snakefiles.
        # This is needed in order to know which method/database combinations
        # are possible and which aren't
        self.taxonomy_formats = {}

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
    if empty_str is not None and len(sorted_cfg) == 0:
        return empty_str
    s = "".join(sorted_cfg).encode('utf-8')
    return sha256(s).hexdigest()



def get_repo_location(**repo_config):
    path = repo_config.get("path", None)
    tag = repo_config.get("tag", None)
    commit = repo_config.get("commit", None)
    # remote
    base_url = "https://github.com/{github}/archive".format(**repo_config)
    if tag is not None:
        url = f"{base_url}/refs/tags/{tag}.zip"
        id_ = tag
    elif commit is not None:
        url = f"{base_url}/{commit}.zip"
        id_ = commit
    else:
        id_ = "local"
        url = None
    assert url is not None or path is not None, \
        "Either tag or commit or path must be defined with uvsnake source"
    return path, url, id_


def download_repo(url, target_dir):
    """
    Downloads the 'uvsnake' pipeline to the working directory.
    This solution was chosen over specifying a remote
    Snakefile with github(...) due to Python modules not being
    included (see also https://github.com/snakemake/snakemake/issues/1632)
    """
    from urllib.request import urlopen
    from io import BytesIO
    import zipfile
    import sys
    import os
    import shutil
    print(f"Downloading {url}...", file=sys.stderr)
    handle = urlopen(url)
    memzip = BytesIO(handle.read())
    archive = zipfile.ZipFile(memzip)
    files = [f for f in archive.namelist() if "/workflow/" in f]
    base_dir = files[0].split("/")[0]
    parent = os.path.dirname(target_dir)
    extr_dir = os.path.join(parent, base_dir)
    os.makedirs(parent, exist_ok=True)
    archive.extractall(parent, files)
    shutil.copytree(extr_dir, target_dir, dirs_exist_ok=True)
    shutil.rmtree(extr_dir)


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
