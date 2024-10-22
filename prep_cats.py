
import os, sys, yaml
from src.cluster_challenge.utils import slurm_submit



config_file = sys.argv[1]
with open(config_file) as fstream :
	cfg = yaml.safe_load(fstream)


job_id = slurm_submit('prepare_catalogs', config_file)

