Bootstrap: docker
From: togop/gcc_gfortran_mpi_hdf5_python:latest

%labels
    Author todor.gitchev@izb.unibe.ch
    Version v2021.6.24

%help

Containing py3DPolyS_LE package

%files
    # maybe better (but it needs to be public) : git clone
    ./CMakeLists.txt
    ./environment.yml
    ./Makefile
    ./MANIFEST.in
    ./mpi.cmake
    ./README.md
    ./requirements_dev.txt
    ./setup.py
    ./py3dpolys_le/3dpolys_le_runner.py
    ./py3dpolys_le/3dpolys_le_stats.py
    ./py3dpolys_le/__init__.py
    ./py3dpolys_le/_logging.py
    ./py3dpolys_le/_version.py
    ./py3dpolys_le/hic_analysis.py
    ./py3dpolys_le/hic_converters.py
    ./py3dpolys_le/job_runner.py
    ./py3dpolys_le/plot_hic.py
    ./py3dpolys_le/src/3dpolys_le.f03
    ./py3dpolys_le/src/analyse.f03
    ./py3dpolys_le/src/kinds.f03
    ./py3dpolys_le/src/lattice_data.f03
    ./py3dpolys_le/src/lib_conf.f90
    ./py3dpolys_le/src/logging.f03
    ./py3dpolys_le/src/polymer_model.f03
    ./py3dpolys_le/src/randomnumber.f03
    ./py3dpolys_le/src/timers.f03
    ./py3dpolys_le/bin/cmd.sh
    ./py3dpolys_le/data/sample_3dpolys_le.cfg
    ./py3dpolys_le/data/ce/dcc_mex-sites_boundaries.csv
    ./py3dpolys_le/data/ce/dcc_rex-sites_Crane2015_bindings.csv
    ./py3dpolys_le/data/ce/input.cfg
    ./test/data/wt_N2_Moushumi2020_HIC1_5000.cool

%environment
     SINGULARITYENV_APPEND_PATH=/opt/miniconda3/bin:
     export SINGULARITYENV_APPEND_PATH

%post
    # Install py3dpolys_le package
    PATH=/opt/miniconda3/bin:$PATH
    export PATH=$PATH
    echo 'export PATH=$PATH' >> $SINGULARITY_ENVIRONMENT
    make env
    . /opt/miniconda3/bin/activate py3dpolys_le
    conda install cmake
    pip install Cython
    make build
    pip install -e /
    echo ". /opt/miniconda3/bin/activate py3dpolys_le" >> $SINGULARITY_ENVIRONMENT
    #echo "conda activate py3dpolys_le" >> $SINGULARITY_ENVIRONMENT

%startscript
    . /opt/miniconda3/bin/activate py3dpolys_le

%runscript
    echo "Container py3DPolyS_LE was created $NOW"
    echo "Arguments received: $*"
    exec echo "$@"
