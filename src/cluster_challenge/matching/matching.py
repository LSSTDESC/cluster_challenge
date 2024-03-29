import time
import subprocess
import os, sys, yaml

current_dir = os.path.dirname(__file__)


def slurm_submit(task, config) :#cat1, cat2, mt_method, mt_pref, mt_params) :
	slurm_mem = 32
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

	cmd  = f"sbatch "
	
	## sbatch options
	## \\\\\\\\\\\\\\\\\\V//////////////////
	cmd += f"--job-name=matching "
	cmd += f"-t  0-05:00 "
	cmd += f"-n 2 "
	cmd += f"--mem {slurm_mem}G "
	cmd += f"-D {current_dir} "
	cmd += f"-L sps "
	cmd += f"-o {slogs_path}{'_'.join([str(p) for p in mt_params])}.out"
	## //////////////////A\\\\\\\\\\\\\\\\\\

	cmd += f"<<?\n"
	cmd += "#!/usr/bin/bash\n"
	cmd += f"source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh\n"
	cmd += f"conda activate /sps/lsst/users/rsolomon/conda_envs/desc_dev/\n"

	cmd += f"python {current_dir}/{task}.py {config}"

	res = subprocess.run(cmd, shell=True, capture_output=True)
	job_id = str(res.stdout).split('batch job ')[1].split('\\')[0]
	print(f"Matching {job_id}")

	return 0

config_file = sys.argv[1]

slurm_submit('clevar_matching', config_file)

### MATCH COSMODC2 AND WAZP.COSMODC2.TRUEZ
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'wazp_cosmoDC2.truez_v0', 'proximity', 'more_massive', [0.05, 1])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'wazp_cosmoDC2.truez_v0', 'proximity', 'angular_proximity', [0.05, 1])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'wazp_cosmoDC2.truez_v0', 'member', 'more_massive', [0.0])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'wazp_cosmoDC2.truez_v0', 'member', 'more_massive', [0.1])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'wazp_cosmoDC2.truez_v0', 'member', 'shared_member_fraction', [0.0])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'wazp_cosmoDC2.truez_v0', 'member', 'shared_member_fraction', [0.1])


### MATCH COSMODC2 AND WAZP.COSMODC2.FZB
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'wazp_cosmoDC2.fzb_v0', 'proximity', 'more_massive', [0.05, 1])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'wazp_cosmoDC2.fzb_v0', 'proximity', 'angular_proximity', [0.05, 1])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'wazp_cosmoDC2.fzb_v0', 'member', 'more_massive', [0.0])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'wazp_cosmoDC2.fzb_v0', 'member', 'more_massive', [0.1])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'wazp_cosmoDC2.fzb_v0', 'member', 'shared_member_fraction', [0.0])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'wazp_cosmoDC2.fzb_v0', 'member', 'shared_member_fraction', [0.1])

### MATCH COSMODC2 AND PYWAZP.COSMODC2.TRUEZ
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'pywazp_cosmoDC2.truez_v0', 'proximity', 'more_massive', [0.05, 1])

### MATCH COSMODC2 AND WAZP.DC2.FZB
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'wazp_DC2.fzb_v0', 'proximity', 'more_massive', [0.05, 1])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'wazp_DC2.fzb_v0', 'proximity', 'angular_proximity', [0.05, 1])
#
#
### MATCH COSMODC2 AND REDMAPPER.COSMODC2
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'redmapper_cosmoDC2_v0', 'proximity', 'more_massive', [0.05, 1])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'redmapper_cosmoDC2_v0', 'proximity', 'angular_proximity', [0.05, 1])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'redmapper_cosmoDC2_v0', 'member', 'more_massive', [0.0])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'redmapper_cosmoDC2_v0', 'member', 'more_massive', [0.1])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'redmapper_cosmoDC2_v0', 'member', 'shared_member_fraction', [0.0])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'redmapper_cosmoDC2_v0', 'member', 'shared_member_fraction', [0.1])
#
#
### MATCH COSMODC2 AND REDMAPPER.DC2
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'redmapper_DC2_v0', 'proximity', 'more_massive', [0.05, 1])
#slurm_submit('clevar_matching', 'cosmoDC2_v0', 'redmapper_DC2_v0', 'proximity', 'angular_proximity', [0.05, 1])


