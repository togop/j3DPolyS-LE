# 3DPolyS-LE

3D Polymer Simulation of chromosome folding by modeled loop extrusion, boundary elements and loading sites.

## Citing

Todor Gitchev, Gabriel Zala, Peter Meister, Daniel Jost, 3DPolyS-LE: an accessible simulation framework to model the interplay between chromatin and loop extrusion, Bioinformatics, Volume 38, Issue 24, 15 December 2022, Pages 5454–5456, https://doi.org/10.1093/bioinformatics/btac705

![Figure_1.png](Figure_1.png)
Figure 1. Features provided by 3DPolyS-LE. (A) Input parameters for the simulation framework. For the polymer, its length, as well as the location and permeability of individual loading sites and boundary elements for the loop extrusion factors (LEFs) can be defined. The properties of the LEFs include their mode of extrusion (symmetrical/asymmetrical), the extrusion speed, the number of LEFs per polymer and the capacity of LEFs to cross each other (Z-loop formation). (B) Typical outputs of the simulations: virtual Hi-C data (top) and ChIP-Seq profile (bottom) of loop extruders


# Installation

### Requirements

Packages and libraries:

- **git** client version 2.17.1, only if you use git command to download the repository;
- **gcc** compiler version 7.5.0 or higher;
- **gfortran** compiler version 7.5.0 or higher;
- **MPI** implementation like MPICH and libmpich-dev (Debian/Ubuntu) or openMPI;
- **HDF5** libraries. Debian/Ubuntu: libhdf5-serial-dev or libhdf5-103 libhdf5-cpp-103 libhdf5-dev libhdf5-mpich-dev;
- **GNU make** version 3.81 or higher;
- **CMake** version 3.15.0 or higher;
- **Python** 3.10, all required packages are listed in the requirements.txt file and alternatively in the environment.yml file;
- **UV** version 0.6.16 or higher.

for previous version v2022.9:
 - **Python** 3.7, 3.8. 3.9;
 - **Conda** version 4.8.2 or higher.

Make sure you have installed or loaded the required libraries.

Typically, gfortran is part of gcc.
If missed, on an HPC cluster you can check if available and load the latest version:
```
module avail gcc
module load gcc/8.2.0
```
On an Ubuntu/Debian Linux it can be installed like this:
```
sudo apt-get install gfortran
# or
sudo apt-get install gcc
```

On MacOs, you can use:
```
brew install gcc
```

For example on an HPC cluster (Slurm) you might need to load the following modules:

###### UV (https://github.com/astral-sh/uv)
```
curl -LsSf https://astral.sh/uv/install.sh | sh
```
Or, from PyPI:
```
pip install uv
```

###### HDF5 (https://www.hdfgroup.org/solutions/hdf5/)
```
module load HDF5
```
Alternatively, on an Ubuntu/Debian Linux could be installed like hits:
```
uv pip install h5py
```

On MacOs, you can use:
```
brew install hdf5
```

###### MPI (Message Passing Interface)
```
# https://www.open-mpi.org/
module load OpenMPI

# OR https://www.mpich.org/
module load mvapich2
```

On an Ubuntu/Debian Linux could be installed like hits:
```
sudo apt-get install mpich
```

On MacOs, you can use:
```
brew install mpich
```

###### CMake (https://cmake.org/)
On an HPC cluster (Slurm) you might need to load like this:
```
module load CMake
```
Alternatively, you can install it using UV:
```
uv pip install cmake
```

On MacOs, you can use:
```
brew install cmake
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
make venv
make all
```

#### Alternative installation from the Python Package Index (PyPi) repository.
With prebuild binaries for some Linux and macOS system environments.

If you already use UV as a package manager:
```
curl -LsSf https://astral.sh/uv/install.sh | sh
```
For more datails see https://docs.astral.sh/uv/getting-started/installation/

Creating a dedicated Python environment is also recommended:
```
uv venv
source .venv/bin/activate
```

Installing dependencies (tested on _uv 0.6.16_):
```
uv pip install cooler pyranges numpy pandas matplotlib scipy filelock h5py dask seaborn openmpi setuptools
```

Installing 3DPolyS-LE using pip (experimental):
```
export UV_LINK_MODE=copy
uv pip install py3dpolys-le
```

For a development release, specifying the exact version, use for example:
```
export UV_LINK_MODE=copy
uv pip install -i https://test.pypi.org/simple/ py3dpolys-le==2025.2.dev3
```
On macOS, it is recommended to build from source (see "2. Build and install").

### 3. Test installation

Check commands help with the following commands:

```
3dpolys_le -h
3dpolys_le_runner -h
3dpolys_le_stats -h
plot_hic -h
plot_sim_stats -h
```

Download test configuration files from https://gitlab.com/togop/3DPolyS-LE/-/tree/master/test?ref_type=heads , 
and unzip them in our working folder. You should have a subfolder `test` with all demo configuration files. 

Run a demo simulation:
```
mpirun 3dpolys_le -o:./out/demo_run test/demo_run_shell.cfg
3dpolys_le -o:./out/demo_run -a:./out/demo_run/r2.84 test/demo_run_shell.cfg
3dpolys_le_stats -o ./out/demo_run -a ./out/demo_run/r2.84 -i test/demo_run_shell.cfg -f sim_stats.csv
3dpolys_le_runner multi_decay_plot -o ./out/demo_run -a ./out/demo_run/r2.84 -i test/demo_run_shell.cfg
```

If everything was correctly installed, the complete help should be printed out. Otherwise, check the '6.Troubleshooting' section below.

### 4. Additional outputs

![Figure_2.png](Figure_2.png)
Figure 2. (A) Grid-simulations run with a `3dpolys_le_runner grid_nlef_km --nlef_list --km_list` command and plotted with 
a `plot_sim_stats --stats_file` command to find the best set of parameters that fit the target data. For each parameter set, a Chi2-score is estimated. Example of optimization by varying the
LEF density and extruding speed with synthetic human Hi-C data as target (see Paper's Supp. Methods). (B) Best Hi-C map model predictions from the simulations in (A) (lower part) compared to
target Hi-C data (upper part).

### 5. Example scenarios

![Figure_3.png](Figure_3.png)
Figure 3. Examples of simulated Hi-C maps for several loop extrusion scenarios.

#### 5.1. Continue option

Any simulation's configuration can be continued with different configuration scenario for polymer with the same length (parameter *Nchain*) 
by using the *init_mode* parameter (*im:s=<sim_output_folder>*).
Where the <sim_output_folder> is the output folder of the simulation which will be used as initializing polymer configuration.

Example in a configuration file (**.cfg**):

```
[3dpolys_le]
...
init_mode = s=./out/test_init
...
```

Passing a parameter directly to the *3dpolys_le*, overwriting the configuration file's value:
```
mpirun 3dpolys_le -im:s=./out/test_init -o:./out/test_continue test_continue.cfg 
```

Passing a parameter over *3dpolys_le_runner* to the *3dpolys_le*, overwriting the configuration file's value:
```
3dpolys_le_runner -im s=./out/test_init -o ./out/test_continue -i test_continue.cfg 
```

### 6. Troubleshooting
Depending on your installation environment, you might want to create a dedicated Python environment.

Go to the cloned repository project's folder:
```
cd 3DPolyS-LE 
```

Build the default *3DPolyS-LE*'s Python environment *py3dpolys_le*:
```
make venv
```
If your default Python version is a bit old you might need to specify a newer version.
In this case, you can install the *py3dpolys_le* like that:
```
uv venv --python 3.12
# or
source .venv/bin/activate
uv pip install cooler pyranges numpy pandas matplotlib scipy dask h5py filelock seaborn openmpi setuptools build cmake 
```

Finally, build and install:
```
make all
```

Check your installation as shown in the 'Test installation' section above. 

- In case of such an error:
```
    from py3dpolys_le.3dpolys_le import main
                     ^
SyntaxError: invalid syntax
```
A possible reason could be incompatibility of Python and Pip versions.
We recommend *Pip* 21.2.4 as prove to be compatible version with the current .
You can install desired *Pip* version like that:
```
python -m pip install pip==21.2.4
```
If this doesn't help, consider upgrading to Python 3.10.


- In case of such an error:

```
AttributeError: module 'numpy' has no attribute 'object'. Did you mean: 'object_'?
```

Try downgrading the *numpy* package to a version bellow 1.24 like _numpy=1.23.5_

- In case of such an error, and you are on a linux terminal:

```
qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
```

Try the following command, that also could be included in your ~/.bashrc file:

```
export QT_QPA_PLATFORM=offscreen
```

- Cleaning your build `cmake-build`or any other temporary folders can also help:


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
# optional interaction_sites.csv/tsv: name,position,length,state(1)
#interaction_sites=
# optional: interaction enegry (<0.) used by interaction_sites
# Ei = -1.
# optional lef_loading_sites.csv/tsv: name,position,length,factor
lef_loading_sites = py3dpolys_le/data/ce/dcc_rex-sites_Crane2015_bindings.csv
basal_loading_factor = 0.
# optional boundaries.csv/tsv: name,midpoint,impermeability,score,b-position,strand
boundary = py3dpolys_le/data/ce/dcc_mex-sites_boundaries.csv
boundary_direction = 0
z_loop = true
unidirectional = false

# analysis: experiments in silico:
# 1.42 = 100nm
radius_contact = 2.84
chrom = chrX
# optional: if present will trigger hic3d output with the given factor for a resolution reduction 
# hic3d_factor = 5

# hic-chi2-min:
cmp_chrs=chrX,X,6
exp_cool=./test/data/wt_N2_Moushumi2020_HIC1_5000.cool
tads_boundary=py3dpolys_le/data/ce/tad_boundaries/N2.chrX.allValidPairs.hic.5-10kbLoops.bed
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

Configuration file sections:

*3dpolys_le* comprise all parameters for running simulations and data analysis (see above).

*job_runner* comprise general parameters for scheduling simulation pipeline steps:

- _cmd_run_ defines how to treat the generated steps commands with valid values:
- _shell_ : execute commands
- _stdout_ : print out commands to the standard output
- _file:<file_path>_ : save commands into a file. It can be overwritten by a passed 3dpolys_le_runner’s --cmd_run_file argument.
- _jobid_re_ defines a regular expression to extract a batch job identifier out of an HPC batch runner output.
- _cmd_job_dependency_ defines a template to add an HPC batch job dependency, where {jobid} is a placeholder for the dependency job extracted by the jobid_re regular expression.
- _cmd_prefix_ defines the HPC batch command prefix to be used for starting a job, where {cmd_job_dependency} is a placeholder for dependency jobs as built by the cmd_job_dependency parameter.

*job_runner_sim* contains general parameters for scheduling the first model simulation step:

- cmd_prefix same as in the [job_runner] section but specific for this kind of jobs.

*job_runner_analysis* contains general parameters for scheduling simulation output data analysis HPC batch jobs for generating predicted Chip/HiC/HiC3D data:

- _cmd_prefix_ same as in the [job_runner] section but specific for such kind of HPC batch jobs.

*job_runner_stats* contains general parameters for scheduling the simulation output data analysis  HPC batch jobs for calculating hic-hic2-min score and generating additional plots:

- _cmd_prefix_ same as in the [job_runner] section but specific for such kind of HPC batch jobs.
- _plot_format_ defines plot output format. Possible values supported by Python’s matplotlib like: png, tif, svg.
- _plot_cmap_ defines plotting color pallet. Possible values supported by Python’s matplotlib like: YlGnBu, cool, hot_r, gist_heat_r, afmhot_r, YlOrRd, Greys, gist_yarg.

To run each step separately on a personal computer without utilizing an HPC batch system, you can use such [job_runner*] configuration:

```
[job_runner]
cmd_run=stdout

[job_runner_sim]

[job_runner_analysis]

[job_runner_stats]

```

With such a configuration file, you can run the steps separately like this:

```
# 1) running a simulation
3dpolys_le -o:/paht/to/sim_output_folder /path/to/input.cfg

# To get advantage of multiprocessing and parallelizing trajectory simulations, it is recommended to run this command with an MPI runner:
mpirun 3dpolys_le -o:/path/to/sim_output_folder /path/to/input.cfg 

# 2) running analyse step: generating predicted Chip, HiC/HiC3D
3dpolys_le -o:/paht/to/sim_output_folder -a:/paht/to/analyse_output_folder /path/to/input.cfg

# 3) running predicted data analyse step: calculating hic-chi2-min score and storing it into a file /path/to/sim_stats.csv
3dpolys_le_stats -o /paht/to/sim_output_folder -a /paht/to/analyse_output_folder -i /path/to/input.cfg -f /path/to/sim_stats.csv

# 4) generating additional plots:  contact-decay  comparing simulation with an experimental HiC data
3dpolys_le_runner multi_decay_plot -o /paht/to/sim_output_folder -a /paht/to/analyse_output_folder -i /path/to/input.cfg

```

Alternatively, you can use a singularity image https://cloud.sylabs.io/library/todor/default/py3dpolys_le to run the above commands.
```
singularity pull library://todor/default/py3dpolys_le:latest
```
Afterwards, you can run the above commands using the following prefix:
```
singularity exec -H $HOME -B $PWD py3dpolys_le_latest.sif <my 3dpolys_le command>
```

You might need to adjust the -B parameter to specify bind paths used in the command or the configuration.

Additionally, you can add a batch command prefix to run it in your HPC like IBM’s LSF for example:

```
bsub -n 12 -R "rusage[mem=8192]" mpirun "<my singularity 3dpolys_le command>"
```

For a SLURM HPC system, the batch command prefix could be like this:

`sbatch --mem-per-cpu=8G --nodes=1 --ntasks-per-node=12 cmd.sh mpirun <my singularity 3dpolys_le command>`

Where the content of the cmd.sh file, needed to overcome some SLURM constrains, is simply:

```
#! /bin/bash
"$@"
```

Demo data and example configuration can be found here https://gitlab.com/togop/3DPolyS-LE/-/blob/master/test/demo_run_shell.cfg and the corresponding commands to run the demo https://gitlab.com/togop/3DPolyS-LE/-/blob/master/test/demo_run_commads.txt (update paths accordingly to your environment) with the needed data files in https://gitlab.com/togop/3DPolyS-LE/-/tree/master/test/data .

Be aware, that only the single steps are supported by the singularity image for now, and NOT all  3dpolys_le_runner’s sub-commands (i.e. '3dpolys_le_runner run|grid_nlef_km|new_stats|contact_radius_analysis'). For the unsupported commands you can use the ‘cmd_run’ option in the configuration file (cmd_run=stdout or cmd_run=file:<file_path>) to generate all step commands and execute them manually afterwards. The singularity image supports and has been tested only for OpenMPI, which needs to be available on the host machine.


With a properly configured input.cfg file for your HPC (so far tested only on Slurm) you can start a simulation job including all the above steps with the following command:

 `3dpolys_le_runner run -i my_sim_input.cfg -o ./my_sim_out`

It will start a series of commands including simulation, analysis, and downstream statistical analysis (hic-chi2-min score) and plots (hic, contact-decay).

It is also helpful to save the output of the main 'run' command as it will print out all executed commands and in case of some errors you can rerun the failed one. One way to do that is to save the output in a file:

`3dpolys_le_runner run -i my_sim_input.cfg -o ./my_sim_out &> 3dpolys_le_runner.log`

In case you want to generate a shell script and execute the single steps one by one, you can generate the run shell script by:

`3dpolys_le_runner run -i my_sim_input.cfg -o ./my_sim_out --cmd_run_file run_my_sim.sh`

To see all supported parameters, run the following command:

 `3dpolys_le_runner --help`

For using the Python wrapper, triggering data analysis (predicted ChIP, HiC/HiC3D, chi2-min score) steps afterwards.

Alternatively, the simulation engine directly is also available via:

  `3dpolys_le -h`

Other available commands are:

  `3dpolys_le_stats --help`

  `plot_hic --help`

  `plot_sim_stats --help`

  `hdf5_to_cooler --help`

  `hic_converters --help`

## Presentation

Todor Gitchev: "Dynamic modelling of chromosome folding in C. elegans suggests in vivo z-loops formation" - INC consortium

https://youtu.be/GtcMC4QTvSo


