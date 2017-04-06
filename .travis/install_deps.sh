#!/bin/bash
set -exuo pipefail

HTSLIB_VERSION=1.4

# Download HTSlib sources and build ---------------------------------------------------------------

if [[ ! -e dep-inst/bin/inst/lib/libhts.so ]]; then
    mkdir -p tmp && pushd tmp

    wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
    tar xf htslib-${HTSLIB_VERSION}.tar.bz2

    cd htslib-${HTSLIB_VERSION}
    ./configure --prefix=$(readlink -f ../../inst-deps)
    make
    make install

    popd
fi
