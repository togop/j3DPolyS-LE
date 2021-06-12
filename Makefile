.PHONY: all build debug env install test doc

build:
	cmake -G "CodeBlocks - Unix Makefiles" -Bcmake-build -S .
	cmake --build cmake-build --target 3dpolys_le -- -j 6

debug:
	cmake -G "CodeBlocks - Unix Makefiles" -Bcmake-build-debug -S . -DCMAKE_BUILD_TYPE=Debug
	cmake --build cmake-build-debug --target 3dpolys_le -- -j 6

env:
	# conda deactivate
	# conda env remove -n py3dpolys_le
	conda env create -f environment.yml
	# conda create -n py3dpolys_le python=3.8
	# conda activate py3dpolys_le

install:
	# conda config --add channels bioconda
	# conda config --add channels conda-forge
	# conda config --add channels defaults
	pip install -e .

doc:
	sphinx-build -b html source build

test:
	# under construction
	pip install pytest
	pip install pytest-cov
	# pip install coverage
	pytest --cov=py3dpolys_le test/

all: build install
