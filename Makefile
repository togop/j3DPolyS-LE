.PHONY: all build debug env install test doc sif

build:
	cmake -G "CodeBlocks - Unix Makefiles" -Bcmake-build -S .
	cmake --build cmake-build --target 3dpolys_le -- -j 6

debug:
	cmake -G "CodeBlocks - Unix Makefiles" -Bcmake-build-debug -S . -DCMAKE_BUILD_TYPE=Debug
	cmake --build cmake-build-debug --target 3dpolys_le -- -j 6

env:
	conda env create -f environment.yml

sif:
	singularity build --force py3DPolyS-LE.sif Singularity

install:
	pip install -e .

doc:
	sphinx-build -b html doc build_doc

test:
	# under construction
	conda install pytest
	conda install pytest-cov
	pytest --cov=py3dpolys_le test/

uninstall:
	conda env remove --name py3dpolys_le

clean:
	cmake --build cmake-build --target clean
	cmake --build cmake-build-debug --target clean

all: build install
