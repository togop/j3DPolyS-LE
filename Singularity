Bootstrap: docker
From: togop/gcc_gfortran_mpi_hdf5_python:latest

%help

Containing py3DPolyS_LE package

%files
    ./Makefile
    ./MANIFEST.in
    ./CMakeLists.txt
    ${SYLABS_TOKEN} ~/.singularity/sylabs-token

%runscript
# Install py3dpolys_le package
    echo "PATH=$PATH:/root/miniconda3/bin/" >> ~/.bashrc
    source ~/.bashrc
    conda init bash
    source ~/.bashrc
    make env
    conda activate py3dpolys_le
    make all
    echo "conda activate py3dpolys_le" >> ~/.bashrc
    source ~/.bashrc
