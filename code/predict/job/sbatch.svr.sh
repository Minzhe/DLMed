#!/bin/bash
#SBATCH --job-name=svr.pred.sh
#SBATCH --partition=256GB
#SBATCH --nodes=1
#SBATCH --time=30-00:00:00
#SBATCH --output=./svr.pred.sh.log
#SBATCH --error=./svr.pred.sh.error

python ../run.predict.svr.py
