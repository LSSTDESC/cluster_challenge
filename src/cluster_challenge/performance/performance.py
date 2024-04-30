
import time
import subprocess
import os, sys, yaml
import src.saving_figures as sfigs

current_dir = os.path.dirname(__file__)


def slurm_submit(task, config, job_list=[]) :
	slurm_mem = 8
	time.sleep(5)

	with open(config) as fstream :
		cfg = yaml.safe_load(fstream)

	cat1 = cfg['cats']['cat1']
	cat2 = cfg['cats']['cat2']
	
	mt_method = cfg['matching']['method']
	mt_pref   = cfg['matching']['pref']
	if mt_method == 'member' :
		mt_params = [cfg['matching']['minimum_share_fraction']]
	else :
		mt_params = [cfg['matching']['delta_z'], cfg['matching']['match_radius']]
	
	slogs_path = f"{current_dir}/slurm_outputs/{cat1}_{cat2}/{mt_method}/{mt_pref}/"
	if not os.path.exists(slogs_path) :
		os.makedirs(slogs_path)

	dep = ""
	if job_list != [] :
		dep = " --dependency=afterany"
		for job in job_list:
			dep += f':{job}'

	cmd  = f"sbatch --partition=lsst,htc --job-name=plotting -t 0-05:00 -n 2 --mem {slurm_mem}G -D {current_dir} -L sps -o {slogs_path}{task}.out <<?\n"
	cmd += "#!/usr/bin/bash\n"
	cmd += f"source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh\n"
	cmd += f"conda activate /sps/lsst/users/rsolomon/conda_envs/desc_dev/\n"
	
	cmd += f"python {current_dir}/{task}.py {config}"

	res = subprocess.run(cmd, shell=True, capture_output=True)
	job_id = str(res.stdout).split("batch job ")[1].split("\\")[0]
	print(job_id)
	
	return job_id


cfg = sys.argv[1]
#job_id = slurm_submit('mass_richness', cfg)
job_id = slurm_submit('scaling', cfg)
slurm_submit('redshift', cfg, job_list=[job_id])
slurm_submit('matching_metrics', cfg)
