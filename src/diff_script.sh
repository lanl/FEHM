##!/bin/bash
# Usage: diff_script.sh file1 file

# this produces error
## Check for the correct number of arguments
# if [[ "$##" -eq 2 ]]; then
#     echo "Usage: $0 file1 file2"
#     exit 1
# fi

file1=$1
file2=$2

## Check if files exist
if [[ ! -f "$file1" ]]; then
    echo "File 1 does not exist: $file1"
    exit 1
fi

if [[ ! -f "$file2" ]]; then
    echo "File 2 does not exist: $file2"
    exit 1
fi

## Store the differential output in an array, using -w to ignore whitespace
diff_output=$(diff -uw "$file1" "$file2")

## Count the lines that differ
num_diff_lines=$(echo "$diff_output" | grep -c '^[-+]' )

## Output only if there are differing lines
if [[ "$num_diff_lines" -gt 0 ]]; then
    echo "$file1 Number of differing lines: $num_diff_lines"
#    echo "Differing lines (up to 5):"
    echo "$diff_output" | grep '^[-+]' | head -n 10 
    echo " "
fi

