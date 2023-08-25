
import os
import yaml

from utils import file_logging


def write_db_config(dbconfig, outfile, exclude=None):
    if exclude is not None:
        for k in exclude:
            try:
                del dbconfig[k]
            except KeyError:
                pass
    yml = yaml.safe_dump(dbconfig)
    # check for the case of a hash collision after updating a ref. database
    # within the same project. Of course, the probability is close to zero
    # and this will not detect collisions across projects.
    if os.path.exists(outfile):
        with open(outfile) as f:
            assert f.read() == yml, """
                "Potential hash collision for database {}. Try deleting the "
                "contents of 'refdb', empty the Snakemake cache or rename the "
                "database in taxonomy.yaml
                """.format(dbconfig["name"])
    with open(outfile, "w") as out:
        out.write(yml)


with file_logging(snakemake.log[0]):
    write_db_config(snakemake.params.dbconfig, snakemake.output.yml, 
                    exclude=snakemake.params.get("exclude", None))
