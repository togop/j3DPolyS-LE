#! /bin/bash

## Allocate resources
#SBATCH --job-name="3dpolys-le_analysis"
#SBATCH --mail-user=todor.gitchev@izb.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
##  SBATCH --ntasks-per-node=20
##  SBATCH --cpus-per-task=20

echo "CALL: mpirun_3dpolys-le.py $@"
mpirun_3dpolys-le.sh "$@"