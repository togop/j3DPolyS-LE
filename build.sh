#!/usr/bin/env bash

CMAKE_VERSION=$(cmake --version | grep version | awk '{print $3}')
echo "cmake version: ${CMAKE_VERSION}"

mkdir -p bin

if [[ "${CMAKE_VERSION}" > "3.13.0" ]]; then
    echo "cmake version > 3.13.0"
    cmake -G "CodeBlocks - Unix Makefiles" -Bcmake-build -S . # -DCMAKE_BUILD_TYPE=Debug
    cmake --build cmake-build --target 3dpolys-le -- -j 6
else
    echo "cmake version < 3.13.0"
    cmake -G "CodeBlocks - Unix Makefiles" # -DCMAKE_BUILD_TYPE=Debug
    cmake --build . --target 3dpolys-le -- -j 6
fi
