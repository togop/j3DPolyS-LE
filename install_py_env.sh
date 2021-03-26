#!/bin/bash

echo "Install 3dpolys-le environment"

#conda env remove -n 3dpolys-le

conda env create -f environment.yml

BIN_DIR="$(pwd)/bin"
sed -i.bak '/export PATH=".*3DPolyS-LE\/bin"/d' ~/.bashrc
echo "export PATH=\"\$PATH:${BIN_DIR}\"" >> ~/.bashrc

echo "The following line was added to your ~/.bashrc:"
echo "export PATH=\"\$PATH:${BIN_DIR}\""

# useful to extend our .bashrc
#export PATH="$PATH:${BIN_DIR}"
#alias a3le='conda activate 3dpolys-le'
#alias d3le='conda deactivate'