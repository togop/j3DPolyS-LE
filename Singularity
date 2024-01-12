Bootstrap: docker
From: togop/gcc_gfortran_mpi_hdf5_python:latest

%labels
    Author todor.gitchev@gmail.com
    Version v2024.1

%help
    3D Polymer Simulation of chromosome folding by modeled loop extrusion
    Available commands:
        3dpolys_le
        3dpolys_le_runner
        3dpolys_le_stats
        plot_hic
        plot_sim_stats
        hdf5_to_cooler
        hic_converters

%environment
     SINGULARITYENV_APPEND_PATH=/opt/miniconda3/bin:
     export SINGULARITYENV_APPEND_PATH
     export PATH=/opt/miniconda3/bin:$PATH

%post
    apt-get update
    # Install py3dpolys_le package
    PATH=/opt/miniconda3/bin:$PATH
    export PATH=$PATH
    echo 'export PATH=/opt/miniconda3/bin:$PATH' >> $SINGULARITY_ENVIRONMENT
    conda update -q conda
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --add channels defaults
    # conda install -y -q --file requirements_dev.txt python=3.8
    conda install -y -q numpy cython pandas matplotlib scipy filelock h5py dask cooler pyranges seaborn python=3.9
    conda install -y -q cmake
    # pip install Cython
    # apt-get update && apt-get install -y git
    # skip update: sometimes cause problems
    apt-get install -y git
    git clone https://gitlab.com/togop/3DPolyS-LE.git -b develop
    cd 3DPolyS-LE
    make build
    pip install -e /3DPolyS-LE

%runscript
    echo "Container py3DPolyS_LE was created $NOW"
    echo "Arguments received: $*"
    exec echo "$@"
