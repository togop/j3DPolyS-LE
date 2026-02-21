.PHONY: all build build-gcc build-gcc-noacc debug env install test doc sif

build:
	rm -rf cmake-build
	cmake -G "CodeBlocks - Unix Makefiles" -Bcmake-build -S . -DENABLE_OPENACC=ON
	cmake --build cmake-build --target 3dpolys_le -- -j 6

# Build with gfortran/gcc only (no nvfortran). Use when mpifort expects nvfortran but it is not in PATH.
# Load gcc, openmpi, and hdf5: e.g. "module purge && module load gcc openmpi hdf5 && make build-gcc"
build-gcc:
	rm -rf cmake-build
	cmake -G "CodeBlocks - Unix Makefiles" -Bcmake-build -S . -DENABLE_OPENACC=ON \
		-DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc
	cmake --build cmake-build --target 3dpolys_le -- -j 6

# Same as build-gcc but without OpenACC. Use if build-gcc fails with "mkoffload: -fopenacc must be set" at link.
# Produces a CPU-only binary (no GPU offload).
build-gcc-noacc:
	rm -rf cmake-build
	cmake -G "CodeBlocks - Unix Makefiles" -Bcmake-build -S . -DENABLE_OPENACC=OFF \
		-DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc
	cmake --build cmake-build --target 3dpolys_le -- -j 6

debug:
	rm -rf cmake-build-debug
	cmake -G "CodeBlocks - Unix Makefiles" -Bcmake-build-debug -S . -DENABLE_OPENACC=ON -DCMAKE_BUILD_TYPE=Debug
	cmake --build cmake-build-debug --target 3dpolys_le -- -j 6

venv:
	uv venv
	uv add -r requirements.txt
	# source .venv/bin/activate
	# uv pip install -r requirements.txt

sif:
	singularity build --force py3DPolyS-LE.sif Singularity

install:
	export UV_LINK_MODE=copy
	uv pip install -e .

package:
	python -m build

doc:
	sphinx-build -b html doc build_doc

test:
	# under construction
	export UV_LINK_MODE=copy
	uv pip install pytest
	uv pip install pytest-cov
	pytest --cov=py3dpolys_le test/

uninstall:
	rm -r .venv

clean:
	cmake --build cmake-build --target clean
	cmake --build cmake-build-debug --target clean

all: build install
