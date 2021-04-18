#! /bin/bash

## Allocate resources
#SBATCH --job-name="3dpolys_le_runner"
##SBATCH --mail-user=todor.gitchev@izb.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
##  SBATCH --ntasks-per-node=20
##  SBATCH --cpus-per-task=20

echo "CALL: 3dpolys_le_runner.py $@"
3dpolys_le_runner.py "$@"