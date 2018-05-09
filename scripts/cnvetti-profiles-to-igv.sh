#!/bin/bash

# Helper script for converting BCF files with coverage windows/profiles to
# the IGV format.

usage() { >&2 echo "Usage: $0 [-f FIELDS] IN.bcf OUT"; }

FIELDS="COV NCOV SCOV NCOV2 SCOV2"

while getopts ":hf:" flag; do
    case "$flag" in
        h) usage; exit 0;;
        f) FIELDS=$OPTARG;;
        "?") >&2 echo "Unknown option $OPTARG";;
        ":") >&2 echo "No value for option $OPTARG";;
        *) >&2 echo "Unknown error while processing"; usage; exit 1;;
    esac

    shift $((OPTIND-1))
done

if [[ $# -ne 2 ]]; then
    usage
    exit 1
fi

IN=$1
OUT=$2

OUT_LIN=
OUT_LOG=

LIN_FIELDS=
LOG_FIELDS=

for field in $FIELDS; do
    case $field in
        COV|NCOV|SCOV)
            OUT_LIN=$OUT.linear.igv
            LIN_FIELDS+=" $field"
            ;;
        NCOV2|SCOV|SCOV2)
            OUT_LOG=$OUT.log2.igv
            LOG_FIELDS+=" $field"
            ;;
        *)
            >&2 echo "ERROR: Invalid metric '$METRIC'";
            usage;
            exit 1
            ;;
    esac
done

header="chr start end feature"

if [[ ! -z "$OUT_LIN" ]]; then
    >&2 echo "Creating $OUT_LIN..."
    echo "#track color=000,000,000 graphType=points viewLimits=0:2.0 windowingFunction=mean" \
    > $OUT_LIN
    echo $header $LIN_FIELDS \
    | tr ' ' '\t' \
    >> $OUT_LIN

    METRICS=%$(echo $LIN_FIELDS | sed -e 's/ /\\t%/g')
    bcftools query \
        -f "%CHROM\t%POS0\t%END\tlinear_coverage[\t$METRICS]\n" \
        $IN \
    >> $OUT_LIN

    >&2 echo "Creating ${OUT_LIN%.igv}..."
    igvtools totdf $OUT_LIN ${OUT_LIN%.igv}.tdf b37
fi

if [[ ! -z "$OUT_LOG" ]]; then
    >&2 echo "Creating $OUT_LOG..."
    echo "#track color=000,000,255 altColor=255,000,000 midColor=000,000,000 graphType=points viewLimits=-1.5:1.5 windowingFunction=mean" \
    >$OUT_LOG
    echo $header $LOG_FIELDS \
    | tr ' ' '\t' \
    >> $OUT_LOG

    METRICS=%$(echo $LOG_FIELDS | sed -e 's/ /\\t%/g')
    bcftools query \
        -f "%CHROM\t%POS0\t%END\tlog_coverage[\t$METRICS]\n" \
        $IN \
    >> $OUT_LOG

    >&2 echo "Creating ${OUT_LOG%.igv}..."
    igvtools totdf $OUT_LOG ${OUT_LOG%.igv}.tdf b37
fi
