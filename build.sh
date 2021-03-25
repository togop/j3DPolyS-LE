#!/bin/bash

CMAKE_VERSION=$(cmake --version | grep version | awk '{print $3}')
echo "cmake version: ${CMAKE_VERSION}"

#mkdir -p bin

cmake -G "CodeBlocks - Unix Makefiles" -Bcmake-build -S . # -DCMAKE_BUILD_TYPE=Debug
cmake --build cmake-build --target 3dpolys-le -- -j 6
