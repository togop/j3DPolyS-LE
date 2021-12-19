# 3DPolyS-LE

3D Polymer Simulation of chromosome folding by modeled loop extrusion, boundary elements and loading sites.

# Installation

### Requirements

Packages and libraries:

- **git** client version 2.17.1, only if you use git command to download the repository;
- **gcc** compiler version 7.5.0 or higher;
- **gfortran** compiler version 7.5.0 or higher;
- **MPI** implementation like MPICH and libmpich-dev (Debian/Ubuntu) or openMPI;
- **HDF5** libraries. Debian/Ubuntu: libhdf5-103 libhdf5-cpp-103 libhdf5-dev libhdf5-mpich-dev;
- **GNU make** version 3.81 or higher;
- **CMake** version 3.15.0 or higher;
- **Python** 3.7, all required packages are listed in the requirements.txt file and alternatively in the environment.yml file;
- **Conda** version 4.8.2 or higher.

Make sure you have installed or loaded the required libraries.

For example on a HPC cluster (Slurm) you might need to load the following modules:

###### Conda (https://conda.io)
```
module load Anaconda3
```
Alternative could be installation of Miniconda (https://docs.conda.io/en/latest/miniconda.html)

###### HDF5 (https://www.hdfgroup.org/solutions/hdf5/)
```
module load HDF5
```

###### MPI (Message Passing Interface)
```
# https://www.open-mpi.org/
module load OpenMPI

# OR https://www.mpich.org/
module load mvapich2
```

On Ubuntu linus could be installed like hits:
```
sudo apt-get install mpich
```

###### CMake (https://cmake.org/)
```
module load CMake
```
Alternatively, you can install it using Conda:
```
conda install cmake
```

This module is required only when you build and install the py3DPlyS-LE package. 

Be aware that all required HPC modules have to be loaded before you run simulations.  

### 1. Clone repository

from the master branch:
`git clone https://gitlab.com/togop/3DPolyS-LE.git`

or from the development branch:
`git clone https://gitlab.com/togop/3DPolyS-LE.git -b develop`

### 2. Build and install 

To build and install as Python package, run the following commands:

```
# go to the cloned repository project folder 
cd 3DPolyS-LE 
make all
```

### Troubleshooting
Depending on your installation environment, you might want to create a dedicated Python environment.

Go to the cloned repository project's folder:
```
cd 3DPolyS-LE 
```

Build the default *3DPolyS-LE*'s Python environment *py3dpolys_le*:
```
make env
```
If your default Python version is a bit old you might need to specify a newer version.
In this case, you can install the *py3dpolys_le* like that:
```
conda env create -f environment.yml python=3.9
# or 
conda env create -f requirements_dev.txt python=3.8
```

Active your *py3dpolys_le* environment:
```
conda activate py3dpolys_le
# or
source activate py3dpolys_le
```
Finally, build and install:
```
make all
```


# Usage

To run a simulation:

Create a copy of an input.cfg file and update the parameters you want.
An example copy of such a configuration file you can find in the package:
https://gitlab.com/togop/3DPolyS-LE/-/blob/develop/py3dpolys_le/data/ce/input.cfg
All simulation's parameters are under section *[3dpolys_le]*, here is an example:
```
[3dpolys_le]
# default 3dpolys_le parameters' values
# polymer characteristics
Nchain = 8860
L = 16
# not used yet: kint = 1.17
Ea = 0.
init_mode = z

# measurements
Niter = 250
Nmeas = 3
Ninter = 840000
burnin = 0
burnout = 0
burnoutM = 0

# Loop-Extrusion factors
kb = 2.8e-6
ku = 2e-6
km = 2.7e-3
Nlef = 200
# optional lef_loading_sites.csv: name,position,length,factor
lef_loading_sites = py3dpolys_le/data/ce/dcc_rex-sites_Crane2015_bindings.csv
basal_loading_factor = 0.
# optional boundaries.csv: name,midpoint,impermeability,score,b-position,strand
boundary = py3dpolys_le/data/ce/dcc_mex-sites_boundaries.csv
boundary_direction = 0
z_loop = true
unidirectional = false

# analysis: experiments in silico:
# 1.42 = 100nm
radius_contact = 2.84

# hic-chi2-min:
cmp_chrs=chrX,X,6
exp_cool=./test/data/wt_N2_Moushumi2020_HIC1_5000.cool
tads_boundary=./test/sip_loopanchor_boundaries.csv
```

Be aware to update properly the *[job_runner]* section according to your system environment.

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

In case you want to generate a shell script and execute the single steps one by one, you can generate the run shell script by:

`3dpolys_le_runner run -i my_sim_input.cfg -o ./my_sim_out --cmd_run_file run_my_sim.sh`

To see all supported parameters, run the following command:

 `3dpolys_le_runner --help`

Other available commands are:

  `3dpolys_le_stats --help`

  `plot_hic --help`

  `plot_sijm_stats --help`


