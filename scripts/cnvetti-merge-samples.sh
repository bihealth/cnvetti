#!/bin/bash

# Helper script for merging coverage BCF file from multiple samples.  The main
# requirement is that all have the same window ID.

usage() { >&2 echo "Usage: $0 -o OUT.bam ONE.bcf TWO.bcf [MORE.bcf]"; }

OUT=

while getopts ":ho" flag; do
    case "$flag" in
        o) OUT=$OPTARG;;
        "h") usage; exit 0;
        "?") >&2 echo "Unknown option $OPTARG";;
        ":") >&2 echo "No value for option $OPTARG";;
        *) >&2 echo "Unknown error while processing"; usage; exit 1;;
    esac

    shift $((OPTIND-1))
done

if [[ -z "$OUT" ]]; then
    >&2 echo "Option '-o' is missing"
    exit 1
fi

O=
if [[ "$OUT" == *.vcf.gz ]]; then
    O="-O z"
elif [[ "$OUT" == *.bcf ]]; then
    O="-O b"
fi

set -x
bcftools merge $O -m id -o "$out" $@
