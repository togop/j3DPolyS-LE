# 3DPolyS-LE

3D Polymer Simulation of chromosome folding by modeled loop extrusion, boundary elements and binding sites.

# Instalation

Required packages and libraries:

- **git** client version 2.17.1;
- **gcc** compiler version 7.5.0 or higher;
- **gfortran** compiler version 7.5.0 or higher;
- **MPI** implementation like MPICH and libmpich-dev (Debian/Ubuntu) or openMPI;
- **HDF5** libraries. Debian/Ubuntu: libhdf5-103 libhdf5-cpp-103 libhdf5-dev libhdf5-mpich-dev;
- **CMake** version 3.13.0 or higher;
- **Python** 3.7, all required packages are listed in the requirements.txt file and alternatively in the environment.yml file;
- **Conda** version 4.8.2 or higher.


from master branch:
`git clone https://gitlab.com/togop/3DPolyS-LE.git`

from development branch:
`git clone https://gitlab.com/togop/3DPolyS-LE.git -b develop`

To install run the following commands:

```
cd 3DPolyS-LE 

# optional or as troubleshooting for problems with python environments
make env
conda activate py3dpolys_le

make all
```

# Usage

To run a simulation:

Create a copy of an input.cfg file and update the parameters you want.
An example copy of such a configuration file you can find in the package:
https://gitlab.com/togop/3DPolyS-LE/-/blob/develop/py3dpolys_le/data/ce/input.cfg

Be aware to update properly the _[job_runner]_ section according to your system environment.

For Slurm environment you can use such a configuration:

```
[job_runner]
cmd_prefix=sbatch --job-name=3dpolys_le --time=5-00:00:00 --mem-per-cpu=6G --nodes=1 --ntasks-per-node=50 --cpus-per-task=1 {cmd_job_dependency}
jobid_re=\d+$
cmd_job_dependency=--dependency=afterany:{jobid}
```

To start the simulation job just run the following command:

 `3dpolys_le_runner run -i my_sim_input.cfg -o ./my_sim_out`
