#!/bin/bash

mkdir build && pushd build

cmake \
    ../../core/CMakeListst.txt \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$PREFIX

make VERBOSE=1

make install
