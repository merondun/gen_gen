#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <lookup_table> <fasta>"
    exit 1
fi

TABLE="$1" #simply a table with OLD_ID\tNEW_ID
FASTA="$2"

#Remove sed if you want to retain the old chr ID, otherwise it will append 
awk 'FNR==NR { a[">"$1]=$2; next } $1 in a { sub(/>/,">"a[$1]"|",$1)}1' "${TABLE}" "${FASTA}" | sed 's/|.*//g'
