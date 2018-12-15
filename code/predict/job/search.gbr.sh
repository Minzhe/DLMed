#!/bin/bash
#SBATCH --job-name=search.gbr.cell.sh
#SBATCH --partition=256GB
#SBATCH --nodes=1
#SBATCH --time=30-00:00:00
#SBATCH --output=./search.gbr.cell.sh.log
#SBATCH --error=./search.gbr.cell.sh.error

source activate s418336
python ../parameter.search/search.gbr.cell.py
