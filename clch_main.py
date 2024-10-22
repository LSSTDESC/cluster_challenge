
import yaml
from src.cluster_challenge.utils import slurm_submit, run_argparse

args = run_argparse()

if 'matching' in args.processes :
	job_id = slurm_submit('matching', args.cfgFile)
else :
	job_id = None

if 'performance' in args.processes :
	job_id = slurm_submit('performance', args.cfgFile, dep=job_id)
