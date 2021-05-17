# 3DPolyS-LE

3D Polymer Simulation of chromosome folding by modeled loop extrusion, boundary elements and binding sites.

# Instalation

Required packages and libraries:

- **git** client version 2.17.1, only if you use git command to download the repository;
- **gcc** compiler version 7.5.0 or higher;
- **gfortran** compiler version 7.5.0 or higher;
- **MPI** implementation like MPICH and libmpich-dev (Debian/Ubuntu) or openMPI;
- **HDF5** libraries. Debian/Ubuntu: libhdf5-103 libhdf5-cpp-103 libhdf5-dev libhdf5-mpich-dev;
- **GNU make** version 3.81 or higher;
- **CMake** version 3.13.0 or higher;
- **Python** 3.7, all required packages are listed in the requirements.txt file and alternatively in the environment.yml file;
- **Conda** version 4.8.2 or higher.

1. Clone repository

from master branch:
`git clone https://gitlab.com/togop/3DPolyS-LE.git`

from development branch:
`git clone https://gitlab.com/togop/3DPolyS-LE.git -b develop`

2. Install 

run the following commands:

```
cd 3DPolyS-LE 

# optional for Slurm environment
module load Anaconda3
# or manual installation of Miniconda

# optional or as troubleshooting for problems with python environments
make env
conda activate py3dpolys_le
# or
source activate py3dpolys_le
# optional or as troubleshooting for problems with curently system installed cmake 
conda install cmake
# or as troubleshooting in case of problems with the conda version  of cmake
module load CMake

# optional for Slurm environment
# If you use this during installation you need to load them always you want to use the package after new login!  
module load HDF5
module load OpenMPI

make all
```

# Usage

To run a simulation:

Create a copy of an input.cfg file and update the parameters you want.
An example copy of such a configuration file you can find in the package:
https://gitlab.com/togop/3DPolyS-LE/-/blob/develop/py3dpolys_le/data/ce/input.cfg

Be aware to update properly the _[job_runner]_ section according to your system environment.

For Slurm environment you can use such a configuration (also could be found in the example input.cfg):

```
[job_runner]
cmd_run=shell
jobid_re=\d+$
cmd_job_dependency=--dependency=afterany:{jobid}
cmd_prefix=sbatch --job-name=3dpolys_le --time=1-00:00:00 --mem-per-cpu=8G --nodes=1 --ntasks-per-node=1 --cpus-per-task=8 {cmd_job_dependency}

[job_runner_sim]
cmd_prefix=sbatch --job-name=sim_3dpolys_le --time=3-00:00:00 --mem-per-cpu=6G --nodes=1 --ntasks-per-node=50 --cpus-per-task=1 {cmd_job_dependency}

[job_runner_analysis]
cmd_prefix=sbatch --job-name=anl_3dpolys_le --time=1-00:00:00 --mem-per-cpu=16G --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 {cmd_job_dependency}

[job_runner_stats]
cmd_prefix=sbatch --job-name=sts_3dpolys_le --time=1-00:00:00 --mem-per-cpu=16G --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 {cmd_job_dependency}
```

To start a simulation job just run the following command:

 `3dpolys_le_runner run -i my_sim_input.cfg -o ./my_sim_out`

It will start a series of commands including simulation, analysis, and downstream statistical analysis (hic-chi2-min score) and plots (hic, contact-decay).

It is also helpful to save the output of the main 'run' command as it will print out all executed commands and in case of some errors you can rerun the failed one. One way to do that is to save the output in a file:

`3dpolys_le_runner run -i my_sim_input.cfg -o ./my_sim_out &> 3dpolys_le_runner.log`


 To see all supported parameters run:

 `3dpolys_le_runner --help`

Other available commands are:

  `3dpolys_le_stats --help`

  `plot_hic --help`

