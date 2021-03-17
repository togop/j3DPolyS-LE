#! /bin/bash

## Allocate resources
#SBATCH --job-name="3dpolys-le"
#SBATCH --mail-user=todor.gitchev@izb.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=6G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=50
#SBATCH --cpus-per-task=1
##  SBATCH --ntasks-per-node=20
##  SBATCH --cpus-per-task=20

mpirun_3dpolys-le.sh "$@"