#!/bin/bash
set -exuo pipefail

HTSLIB_VERSION=1.4
TCLAP_VERSION=1.2.1

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

# Download TCLAP sources and build ----------------------------------------------------------------

if [[ ! -e inst-deps/include/tclap ]]; then
    mkdir -p tmp && pushd tmp

    wget https://sourceforge.net/projects/tclap/files/tclap-${TCLAP_VERSION}.tar.gz/download \
        -O tclap-${TCLAP_VERSION}.tar.gz
    tar xf tclap-${TCLAP_VERSION}.tar.gz

    cd tclap-${TCLAP_VERSION}
    ./configure --prefix=$(readlink -f ../../inst-deps)
    make
    make install

    popd
fi
