#!/bin/bash
set -exuo pipefail

export LD_LIBRARY_PATH=$(readlink -f inst-deps)/lib

./inst/bin/cnvetti --help
./inst/bin/cnvetti --version
