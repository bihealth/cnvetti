#!/bin/bash
set -exuo pipefail

mkdir build && cd build

dep_prefix=$(readlink -f ../inst-deps)
prefix=$(readlink -f ../inst)

cmake ../core \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$prefix \
    -DCMAKE_INCLUDE_DIRECTORIES=$dep_prefix/include \
    -DCMAKE_LINK_DIRECTORIES=$dep_prefix/lib
make VERBOSE=1
make install
