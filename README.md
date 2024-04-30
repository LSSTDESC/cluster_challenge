# cluster_challenge
Comparison of cluster algorithm performances on cosmoDC2 and DC2 datasets   
Developed for the DC2 project: https://portal.lsstdesc.org/DESCPub/app/PB/show_project?pid=248   
Common framework for the analysis of cosmoDC2, redMaPPer, WaZP and AMICO catalogs   
Using CLEVAR package   
For CLEVAR installation: pip install -e clevar  (after python setup.py install --user)

## DESC environment
### Create a DESC conda environment (to be done only once)   
source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh    
conda create --clone desc -p desc_july_2022_v0   
conda activate desc_july_2022_v0   
conda install mysql-connector-python -c conda-forge   
### Setup (to be done every time you log in)  
source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh   
conda activate desc_july_2022_v0  

## Prepare Catalogs
It may be necessary or useful to pull only necessary quantities from the original catalog files and to change their headers.
This can be done using some configuration files to help trace back to the original files and quantities:
```bash
python prep_cats.py prepConfigs/example.cfg
```
The prep_cats.py file will submit a slurm job that prepares the catalogs for you.
A bash script and slurm output file will be created for each run with a unique configuration file.
The configuration file contains information on the catalog location, name, and column names that may be changed.
It also contains slurm job information.

## Run Matching
Once you know the locations and header names of the catalogs you want to compare, matching can be done with a separate configuration file.
```bash
python clch_main.py runConfigs/example.cfg
```
The clch_main.py file can be used to simultaneously submit slurm jobs to run both catalog matching and performance scripts.
The default setting is to only run matching since the performance scripts are not finalized.
To run both you can use the -p option:
```bash
python clch_main.py runConfigs/example.cfg -p matching performance
```
The configuration file here similarly has catalog locations and output locations along with details for the matching procedure.
It also contains the catalog headers (cltags and mbtags).
These tags only need to include the parameters you want CLEVAR to be able to access.
The other headers from the input catalog will be available in the matched catalogs but may not be accessible to certain CLEVAR functions.
