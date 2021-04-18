#!/usr/bin/env bash

## Allocate resources
#SBATCH --job-name="3dpolys_le_stats"
##SBATCH --mail-user=todor.gitchev@izb.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --time=0-05:00:00
#SBATCH --mem-per-cpu=16G
##SBATCH --ntasks=1
## SBATCH --nodes=1
## SBATCH --ntasks-per-node=20
## SBATCH --cpus-per-task=20

echo "CALL: 3dpolys_le_stats.py $@"
3dpolys_le_stats.py "$@"
