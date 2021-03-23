#!/bin/bash

echo "Install 3dpolys-le environment"

# dealing with miniconda vs anaconda installations
CONDA_PACKAGE=$(which conda)
if [ ! "${CONDA_PACKAGE}" ];
then 
	echo "Check if anaconda or miniconda is installed in the home directory";
    exit
fi

#conda init bash
#source ${HOME}/.bashrc
conda env remove -n 3dpolys-le

PYTHON_VERSION="3.8"

# alternative:
# conda env create -f environment.yml
conda create -n 3dpolys-le python=${PYTHON_VERSION}
source ${CONDA_PACKAGE%conda}activate 3dpolys-le

${CONDA_PACKAGE%conda} config --set always_yes yes
${CONDA_PACKAGE%conda} config --set quiet true
${CONDA_PACKAGE%conda} update -q conda
${CONDA_PACKAGE%conda} config --add channels conda-forge
${CONDA_PACKAGE%conda} config --add channels bioconda
${CONDA_PACKAGE%conda} config --add channels r

echo "Install $( cat requirements.txt )"
echo "run: ${CONDA_PACKAGE} install --file requirements.txt"
${CONDA_PACKAGE} install --file requirements.txt

sed -i.bak '/CONDA_ACTIVATE=/d' ~/.bashrc
BIN_DIR="$(pwd)/bin"
sed -i.bak '/export PATH=".*3dpolys-le\/bin"/d' ~/.bashrc
echo "export CONDA_ACTIVATE=${CONDA_PACKAGE%conda}activate" >> ~/.bashrc
echo "export PATH=\"\$PATH:${BIN_DIR}\"" >> ~/.bashrc

# useful to extend our .bashrc
#export PATH="$PATH:${HOME}/dev/3dpolys-le/bin"
#alias ade='conda activate 3dpolys-le'
#alias dde='conda deactivate'