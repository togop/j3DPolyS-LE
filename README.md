# 3DPolyS-LE

3D Polymer Simulation with Loop Extrusion. Formerly dcc-extrusion project.

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


```
cd 3DPolyS-LE 

# optional
make env
conda activate py3dpolys_le

make all
```

