#!/bin/bash
#SBATCH --job-name=mlpr.pred.sh
#SBATCH --partition=256GB
#SBATCH --nodes=1
#SBATCH --time=30-00:00:00
#SBATCH --output=./mlpr.pred.sh.log
#SBATCH --error=./mlpr.pred.sh.error

python ../run.predict.mlpr.py
