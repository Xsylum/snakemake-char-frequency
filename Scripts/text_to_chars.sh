#!/bin/bash

function allCaseCounter () {
    input="$1"
    output="$2"
    declare -A char_occurrence_array

    while IFS= read -r line; do
        same_case_line=$(echo -n "$line" | tr '[:upper:]' '[:lower:]')
        while IFS= read -r -n1 char; do
            [[ -z "$char" ]] && continue # if char is empty, then && checks continue, which skips current loop
            # echo "$char"
            ((char_occurrence_array["$char"]++))
        done <<< "$same_case_line"$'\n'
    done < "$input"

    for key in ${!char_occurrence_array[@]}; do
        echo "$key ${char_occurrence_array[$key]}" >> "$output"
    done
}

allCaseCounter $1 $2