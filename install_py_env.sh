#!/bin/bash

DCC_EXTRUSION_SCRIPT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo "Install DCC-EXTRUSION-env in " ${DCC_EXTRUSION_SCRIPT}

# dealing with miniconda vs anaconda installations
CONDA_PACKAGE=`which conda`
if [ ! ${CONDA_PACKAGE} ]; 
then 
	echo "Check if anaconda or miniconda is installed in the home directory";
    exit
fi

#conda init bash
#source ${HOME}/.bashrc
conda env remove -n DCC-EXTRUSION-env

PYTHON_VERSION="3.8"

# alternative:
# conda env create -f environment.yml
conda create -n DCC-EXTRUSION-env python=${PYTHON_VERSION}
source ${CONDA_PACKAGE%conda}activate DCC-EXTRUSION-env

${CONDA_PACKAGE%conda} config --set always_yes yes
${CONDA_PACKAGE%conda} config --set quiet true
${CONDA_PACKAGE%conda} update -q conda
${CONDA_PACKAGE%conda} config --add channels conda-forge
${CONDA_PACKAGE%conda} config --add channels bioconda
${CONDA_PACKAGE%conda} config --add channels r

echo "Install $( cat requirements.txt )"
echo "run: ${CONDA_PACKAGE} install --file requirements.txt"
${CONDA_PACKAGE} install --file requirements.txt

# only install conda R if there is no local R otherwise you get conflicts
#R_VERSION="3.6.1"
#Rinstalled=`which R`
#if [ -z ${Rinstalled} ]
#then
#  #'wget https://cran.r-project.org/src/base/R-3/R-${R_VERSION}.tar.gz
#  #tar -xzvf  R-${R_VERSION}.tar.gz
#  conda install -c r r-base r-essentials
#  conda install -c bioconda bioconductor-biocinstaller
#else
#  echo ""
#  echo "Detected existing R installation "${Rinstalled}
#fi

sed -i.bak '/CONDA_ACTIVATE=/d' ~/.bashrc
BIN_DIR="$(pwd)/bin"
sed -i.bak '/export PATH=".*dcc-extrusion\/bin"/d' ~/.bashrc
echo "export CONDA_ACTIVATE=${CONDA_PACKAGE%conda}activate" >> ~/.bashrc
echo "export PATH=\"\$PATH:${BIN_DIR}\"" >> ~/.bashrc

# useful to extend our .bashrc
#export PATH="$PATH:${HOME}/dev/dcc-extrusion/bin"
#alias ade='conda activate DCC-EXTRUSION-env'
#alias dde='conda deactivate'