#!/bin/bash
#SBATCH --job-name=search.svr.sh
#SBATCH --partition=256GB
#SBATCH --nodes=1
#SBATCH --time=30-00:00:00
#SBATCH --output=./search.svr.sh.log
#SBATCH --error=./search.svr.sh.error

source activate s418336
python ../parameter.search/search.svr.py
