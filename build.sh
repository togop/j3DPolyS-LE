#!/bin/bash

CMAKE_VERSION=$(cmake --version | grep version | awk '{print $3}')
echo "cmake version: ${CMAKE_VERSION}"

#mkdir -p bin

if [[ "$1" = "debug" ]]; then
  cmake -G "CodeBlocks - Unix Makefiles" -Bcmake-build-debug -S . -DCMAKE_BUILD_TYPE=Debug
  cmake --build cmake-build-debug --target 3dpolys-le -- -j 6
else
  cmake -G "CodeBlocks - Unix Makefiles" -Bcmake-build -S . # -DCMAKE_BUILD_TYPE=Debug
  cmake --build cmake-build --target 3dpolys-le -- -j 6
fi
