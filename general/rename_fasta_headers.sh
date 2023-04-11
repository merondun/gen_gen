#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <lookup_table> <fasta>"
    exit 1
fi

TABLE="$1"
FASTA="$2"

#Remove sed if you want to retain the old chr ID
awk 'FNR==NR { a[">"$1]=$2; next } $1 in a { sub(/>/,">"a[$1]"|",$1)}1' "${TABLE}" "${FASTA}" | sed 's/|.*//g'
