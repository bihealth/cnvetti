#!/bin/bash
set -exuo pipefail

mkdir build && cd build

dep_prefix=$(readlink -f ../inst-deps)
prefix=$(readlink -f ../inst)

cmake ../core \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$prefix \
    -DCMAKE_CXX_FLAGS="-isystem $dep_prefix/include" \
    -DCMAKE_LIBRARY_PATH=$dep_prefix/lib \
    -DCMAKE_INCLUDE_PATH=$dep_prefix/include
make VERBOSE=1
make install
