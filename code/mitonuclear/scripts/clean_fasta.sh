#!/bin/bash

# Check if a file is provided as an argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <fasta_file>"
  exit 1
fi

input_file=$1

# Process the FASTA file
awk '/^>/ {print gensub(/_.*/, "", "g", $0); next} {print}' $input_file
