FROM debian:stable-slim
MAINTAINER Todor Gitchev <todor.gitchev@izb.unibe.ch>

ENV DEBIAN_FRONTEND noninteractive

# install all
RUN apt-get update --yes
RUN apt-get install --yes gcc
RUN apt-get install --yes gfortran
RUN apt-get install --yes cmake
RUN apt-get install --yes mpich libmpich-dev
RUN apt-get install --yes libhdf5-103 libhdf5-cpp-103 libhdf5-dev libhdf5-mpich-dev
RUN apt-get install --yes curl
#RUN curl -LsSf https://astral.sh/uv/install.sh | env UV_INSTALL_DIR="/opt/uv" sh
RUN curl -LsSf https://astral.sh/uv/install.sh | sh

# check all is there
CMD dcc --version
CMD gfortran --version
CMD cmake --version
CMD mpirun --version
CMD h5cc -showconfig

# Install py3dpolys_le package
CMD echo "PATH=$PATH:/root/uv/" >> ~/.bashrc
CMD source ~/.bashrc
CMD source ~/.bashrc
CMD make venv
CMD source .venv/bin/activate
CMD make all
CMD echo "conda activate py3dpolys_le" >> ~/.bashrc
  #CMD cat ~/.bashrc
CMD source ~/.bashrc
