import time
import subprocess
import os

current_dir = os.path.dirname(__file__)

print(f'Running from: {current_dir}')

def slurm_submit(task, cat1, cat2, mt_method, mt_pref, mt_params) :
	slurm_mem = 32
	time.sleep(5)

	cats_str = f"{cat1.replace('_', '.')}_{cat2.replace('_', '.')}"
	slurm_output_path = f"{current_dir}/slurm_outputs/{cats_str}_{mt_method}_{mt_pref}_{'_'.join([str(p) for p in mt_params])}.out"
	
	cmd  = f"sbatch --job-name=plotting_{cats_str} -t  0-05:00 -n 2 --mem {slurm_mem}G -D {current_dir} -L sps -o {slurm_output_path} <<?\n"
	cmd += "#!/usr/bin/bash\n"
	cmd += f"source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh\n"
	cmd += f"conda activate /sps/lsst/users/rsolomon/conda_envs/desc_dev/\n"

	cmd += f"python {current_dir}/../{task}.py {cat1} {cat2} {mt_method} {mt_pref} {' '.join([str(p) for p in mt_params])}"

	print(f"{task}.py \t {cat1} \t {cat2} \t {mt_method} \t {mt_pref} \t {' '.join([str(p) for p in mt_params])}")

	res = subprocess.run(cmd, shell=True, capture_output=True)
	job_id = str(res.stdout).split('batch job ')[1].split('\\')[0]
	print(job_id)

	return 0


### COSMODC2 AND WAZP.COSMODC2
#slurm_submit('make_plots', 'cosmoDC2_v0', 'wazp_cosmoDC2.truez_v0', 'proximity', 'more_massive', [0.05, 1])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'wazp_cosmoDC2.truez_v0', 'proximity', 'angular_proximity', [0.05, 1])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'wazp_cosmoDC2.truez_v0', 'member', 'more_massive', [0.0])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'wazp_cosmoDC2.truez_v0', 'member', 'more_massive', [0.1])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'wazp_cosmoDC2.truez_v0', 'member', 'shared_member_fraction', [0.0])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'wazp_cosmoDC2.truez_v0', 'member', 'shared_member_fraction', [0.1])
#
#
### COSMODC2 AND WAZP.COSMODC2.FZB
#slurm_submit('make_plots', 'cosmoDC2_v0', 'wazp_cosmoDC2.fzb_v0', 'proximity', 'more_massive', [0.05, 1])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'wazp_cosmoDC2.fzb_v0', 'proximity', 'angular_proximity', [0.05, 1])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'wazp_cosmoDC2.fzb_v0', 'member', 'more_massive', [0.0])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'wazp_cosmoDC2.fzb_v0', 'member', 'more_massive', [0.1])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'wazp_cosmoDC2.fzb_v0', 'member', 'shared_member_fraction', [0.0])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'wazp_cosmoDC2.fzb_v0', 'member', 'shared_member_fraction', [0.1])


## COSMODC2 AND WAZP.DC2
slurm_submit('make_plots', 'cosmoDC2_v0', 'wazp_DC2.fzb_v0', 'proximity', 'more_massive', [0.05, 1])
slurm_submit('make_plots', 'cosmoDC2_v0', 'wazp_DC2.fzb_v0', 'proximity', 'angular_proximity', [0.05, 1])


### COSMODC2 AND REDMAPPER.COSMODC2
#slurm_submit('make_plots', 'cosmoDC2_v0', 'redmapper_cosmoDC2_v0', 'proximity', 'more_massive', [0.05, 1])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'redmapper_cosmoDC2_v0', 'proximity', 'angular_proximity', [0.05, 1])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'redmapper_cosmoDC2_v0', 'member', 'more_massive', [0.0])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'redmapper_cosmoDC2_v0', 'member', 'more_massive', [0.1])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'redmapper_cosmoDC2_v0', 'member', 'shared_member_fraction', [0.0])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'redmapper_cosmoDC2_v0', 'member', 'shared_member_fraction', [0.1])
#
#
### COSMODC2 AND REDMAPPER.DC2
#slurm_submit('make_plots', 'cosmoDC2_v0', 'redmapper_DC2_v0', 'proximity', 'more_massive', [0.05, 1])
#slurm_submit('make_plots', 'cosmoDC2_v0', 'redmapper_DC2_v0', 'proximity', 'angular_proximity', [0.05, 1])



