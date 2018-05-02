#!/bin/bash

# Helper script for converting normalized coverage BCF tools, possibly with
# more than one sample into a ".igv" file for convertion to TGF and ultimately
# display in IGV.

usage() { >&2 echo "Usage: $0 [-m RCOV|COV|NCOV|NCOV2|NCOV_GC|NCOV2_GC|SCOV|SCOV2] IN.bcf OUT.igv"; }

METRIC=

while getopts ":hm:" flag; do
    case "$flag" in
        h) usage; exit 0;;
        m) METRIC=$OPTARG;;
        "?") >&2 echo "Unknown option $OPTARG";;
        ":") >&2 echo "No value for option $OPTARG";;
        *) >&2 echo "Unknown error while processing"; usage; exit 1;;
    esac

    shift $((OPTIND-1))
done

case "$METRIC" in
    RCOV|COV|NCOV|NCOV2|NCOV_GC|NCOV2_GC|SCOV|SCOV2);;
    *) >&2 echo "ERROR: Invalid metric '$METRIC'"; usage; exit 1;;
esac

if [[ $# -ne 2 ]]; then
    usage
    exit 1
fi

IN=$1
OUT=$2

header="chr start end feature"
samples=$(
    bcftools view --header-only "$IN" \
    | grep '^#CHROM' \
    | cut -f 10-)

# clear
>$OUT

case $METRIC in
    NCOV2*|SCOV2)
        echo "#track color=000,000,255 altColor=255,000,000 graphType=points viewLimits=-1.5:1.5 windowingFunction=none" \
        >> $OUT
    ;;
    NCOV_*|SCOV)
        echo "#track color=000,000,000 graphType=points viewLimits=0:1.8 windowingFunction=none" \
        >> $OUT
    ;;
    *)
        echo "#track color=000,000,255 graphType=points windowingFunction=none" \
        >> $OUT
    ;;
esac

#echo "type=COPY_NUMBER" \
#>> $OUT

echo "$header $(for sample in $samples; do echo -n "$sample-$METRIC "; done)" \
| tr ' ' '\t' \
>> $OUT

bcftools query \
    -f "%CHROM\t%POS0\t%END\t$METRIC[\t%$METRIC]\n" \
    $IN \
>> $OUT
