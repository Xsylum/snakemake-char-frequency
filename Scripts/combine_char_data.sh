#!/bin/bash

output="$1"
shift # shift arguments by one so that we don't use the output file below
declare -A char_occurrence_array

for input in "$@"; do
    while IFS=" " read -r char occs; do
        char_occurrence_array[$char]=$((${char_occurrence_array[$char]} + $occs))
    done < "$input"
done

> "$output" # empty the output file
for key in ${!char_occurrence_array[@]}; do
    echo "$key ${char_occurrence_array[$key]}" >> "$output"
done