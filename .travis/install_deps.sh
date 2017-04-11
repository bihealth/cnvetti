#!/bin/bash
set -exuo pipefail

HTSLIB_VERSION=1.4
CLI11_VERSION=0.8

# Download HTSlib sources and build ---------------------------------------------------------------

if [[ ! -e inst-deps/lib/libhts.so.$HTSLIB_VERSION ]]; then
    mkdir -p tmp && pushd tmp

    wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
    tar xf htslib-${HTSLIB_VERSION}.tar.bz2

    cd htslib-${HTSLIB_VERSION}
    ./configure --prefix=$(readlink -f ../../inst-deps)
    make
    make install

    popd
fi

# Download CLI11 sources and build ----------------------------------------------------------------

if [[ ! -e inst-deps/include/cli11 ]]; then
    mkdir -p tmp && pushd tmp

    wget -O CLI11-0.8.tar.gz \
        https://github.com/CLIUtils/CLI11/archive/v0.8.tar.gz
    tar xf CLI11-${CLI11_VERSION}.tar.gz

    cd CLI11-${CLI11_VERSION}
    install -d ../../inst-deps/include/CLI
    install -t $(readlink -f ../../inst-deps/include/CLI) include/CLI/*

    popd
fi