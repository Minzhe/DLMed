#!/bin/bash
#SBATCH --job-name=search.gbr.sh
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --time=30-00:00:00
#SBATCH --output=./search.gbr.sh.log
#SBATCH --error=./search.gbr.sh.error

source activate s418336
python ../parameter.search/search.gbr.py
