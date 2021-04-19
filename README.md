# 3DPolyS-LE

3D Polymer Simulation with Loop Extrusion. Formerly dcc-extrusion project.

# Instalation

Required packages and libraries:

- _ **git** client version 2.17.1;_
- _ **gcc** compiler version 7.5.0 or higher;_
- _ **gfortran** compiler version 7.5.0 or higher;_
- _ **MPI** implementation like MPICH and libmpich-dev (Debian/Ubuntu) or openMPI **;** _
- _ **HDF5** libraries. Debian/Ubuntu: libhdf5-103 libhdf5-cpp-103 libhdf5-dev libhdf5-mpich-dev;_
- _ **CMake** version 3.13.0 or higher;_
- _ **Python** 3.7, all required packages are listed in the requirements.txt file and alternatively in the environment.yml file;_
- _ **Conda** version 4.8.2 or higher._


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

