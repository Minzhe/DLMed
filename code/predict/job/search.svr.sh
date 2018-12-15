#!/bin/bash
#SBATCH --job-name=search.svr.cell.sh
#SBATCH --partition=256GB
#SBATCH --nodes=1
#SBATCH --time=30-00:00:00
#SBATCH --output=./search.svr.cell.sh.log
#SBATCH --error=./search.svr.cell.sh.error

source activate s418336
python ../parameter.search/search.svr.cell.py
