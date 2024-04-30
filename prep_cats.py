
import os, sys, yaml
from src.cluster_challenge.utils import slurm_submit


#current_dir = os.path.dirname(__file__)

config_file = sys.argv[1]
with open(config_file) as fstream :
	cfg = yaml.safe_load(fstream)


job_id = slurm_submit('prepare_catalogs', config_file)

#def slurm_submit(task, config, job_list=[]) :
#	slurm_mem = 8
#	time.sleep(5)
#	
#	slogs_path = f"{current_dir}/slurm_outputs/"
#	if not os.path.exists(slogs_path) :
#		os.makedirs(slogs_path)
#	
#	dep = ""
#	if job_list != [] :
#		dep = " --dependency=afterany"
#		for job in job_list:
#			dep += f':{job}'
#	
#	cmd  = f"sbatch "
#	
#	## sbatch options
#	## \\\\\\\\\\\\\\\\\\\\V////////////////////
#	cmd += f"--partition=lsst,htc "
#	cmd += f"--job-name={task} "
#	cmd += f"-t 0-05:00 "
#	cmd += f"-n 2 "
#	cmd += f"--mem {slurm_mem}G "
#	cmd += f"-D {current_dir} "
#	cmd += f"-L sps "
#	cmd += f"-o {slogs_path}{task}.out "
#	## ////////////////////A\\\\\\\\\\\\\\\\\\\\
#	
#	cmd += f"<<?\n"
#	cmd += "#!/usr/bin/bash\n"
#	cmd += f"source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh\n"
#	cmd += f"conda activate /sps/lsst/users/rsolomon/conda_envs/desc_dev/\n"
#	
#	cmd += f"python {current_dir}/src/cluster_challenge/{task}/{task}.py {current_dir}/{config_file}"
#	
#	res = subprocess.run(cmd, shell=True, capture_output=True)
#	job_id = str(res.stdout).split("batch job ")[1].split("\\")[0]
#	print(job_id)
#	
#	return job_id
#
#
#job_ids = slurm_submit(task, cfg)
#
