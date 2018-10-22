#!/bin/bash
#SBATCH --job-name=gbr.pred.sh
#SBATCH --partition=256GB
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --output=./gbr.pred.sh.log
#SBATCH --error=./gbr.pred.sh.error

python ../run.predict.gbr.py
