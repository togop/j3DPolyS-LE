Bootstrap: docker
From: togop/gcc_gfortran_mpi_hdf5_python:latest

%labels
    Author todor.gitchev@izb.unibe.ch
    Version v2021.6.24

%help

Containing py3DPolyS_LE package

%files
    ./test/data/wt_N2_Moushumi2020_HIC1_5000.cool /data

%environment
     SINGULARITYENV_APPEND_PATH=/opt/miniconda3/bin:
     export SINGULARITYENV_APPEND_PATH

%post
    # Install py3dpolys_le package
    PATH=/opt/miniconda3/bin:$PATH
    export PATH=$PATH
    echo 'export PATH=$PATH' >> $SINGULARITY_ENVIRONMENT
    #echo "PATH=$PATH:/opt/miniconda3/bin/" >> ~/.bashrc
    #source ~/.bashrc
    #conda init bash
    #source ~/.bashrc
    make env
    conda activate py3dpolys_le
    make all
    echo "conda activate py3dpolys_le" >> $SINGULARITY_ENVIRONMENT

%startscript
    /opt/miniconda3/bin/conda activate py3dpolys_le

%runscript
    echo "Container py3DPolyS_LE was created $NOW"
    echo "Arguments received: $*"
    exec echo "$@"
