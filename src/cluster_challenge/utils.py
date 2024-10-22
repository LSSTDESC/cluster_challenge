import os, sys, yaml, glob, difflib, subprocess, time
import glob, pickle, struct
import numpy as np
import argparse as argparse
from astropy.io import ascii
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt


def create_slurm_script(task, config):
    with open(config) as fstream:
        param_cfg = yaml.safe_load(fstream)
    slurm_cfg = param_cfg['admin']['slurm']

    print ('...task = ', task)

    scr = __file__.replace('utils.py', '')

    logPath = os.path.join(slurm_cfg['logPath'].replace('TMP', scr+'/../../'), param_cfg['name'])
    scriptPath = os.path.join(slurm_cfg['scriptPath'].replace('TMP', scr+'../../'), param_cfg['name'])

    logFile = os.path.join(logPath, slurm_cfg['logFile'][task])
    script = os.path.join(scriptPath, slurm_cfg['scriptFile'][task])

    if not os.path.exists(logPath) :
        os.makedirs(logPath)

    if not os.path.exists(scriptPath) :
        os.makedirs(scriptPath)


    f = open(f"{script}", "w")
    f.write("#!/bin/sh\n")
    f.write(f"#SBATCH --nodes={slurm_cfg['Nnodes']}\n")
    f.write(f"#SBATCH --job-name={task}\n")
    f.write(f"#SBATCH --time={slurm_cfg['time']}\n")
    f.write(f"#SBATCH --partition=lsst,htc\n")
    f.write(f"#SBATCH --ntasks=4\n")
    f.write(f"#SBATCH --output={logFile}\n")
    f.write(f"#SBATCH --mem={slurm_cfg['memory'][task]}GB\n")
    f.write(f"python -u {os.path.join(scr, task, task+'.py')} {config}\n")
    f.close()
    return script


def slurm_submit(task, config, dep=None):

    if (dep is not None) :
        time.sleep(3)

    script = create_slurm_script(task, config)
    if dep is not None:
        cmd = f"sbatch --depend=afterok:{dep} {script}"
    else:
        cmd = f"sbatch {script}"

    res = subprocess.run(cmd, shell=True, capture_output=True)
    job_id = str(res.stdout).split("Submitted batch job ")[1].split("\\")[0]
    return job_id


def run_argparse() :
	CLI = argparse.ArgumentParser()
	CLI.add_argument('cfgFile', type=str)
	CLI.add_argument('-p', '--processes', nargs='*', type=str, default=['matching'])

	return CLI.parse_args()

