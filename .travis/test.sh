#!/bin/bash
set -exuo pipefail

export LD_LIBRARY_PATH=$(readlink -f inst-deps)/lib

./inst/bin/cnvetti --version
./inst/bin/cnvetti --help
./inst/bin/cnvetti coverage --help
./inst/bin/cnvetti normalize --help
./inst/bin/cnvetti background --help
