#!/bin/bash
#SBATCH --job-name=rfr.pred.sh
#SBATCH --partition=256GB
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --output=./rfr.pred.sh.log
#SBATCH --error=./rfr.pred.sh.error

python ../run.predict.rfr.py
