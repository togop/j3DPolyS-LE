.PHONY: all build debug env install test doc sif

build:
	cmake -G "CodeBlocks - Unix Makefiles" -Bcmake-build -S .
	cmake --build cmake-build --target 3dpolys_le -- -j 6

debug:
	cmake -G "CodeBlocks - Unix Makefiles" -Bcmake-build-debug -S . -DCMAKE_BUILD_TYPE=Debug
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
