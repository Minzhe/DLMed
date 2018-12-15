#!/bin/bash
#SBATCH --job-name=search.rfr.cell.sh
#SBATCH --partition=256GB
#SBATCH --nodes=1
#SBATCH --time=30-00:00:00
#SBATCH --output=./search.rfr.cell.sh.log
#SBATCH --error=./search.rfr.cell.sh.error

source activate s418336
python ../parameter.search/search.rfr.cell.py
