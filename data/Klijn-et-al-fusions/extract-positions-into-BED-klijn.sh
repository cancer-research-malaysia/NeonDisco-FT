#!/bin/bash

# Script to convert genomic fusion/breakpoint data to two BED files using awk
# Usage: ./script.sh input_file.txt
# Output: Creates input_file_5prime.bed and input_file_3prime.bed

input_file="$1"

# Check if input file is provided
if [ -z "$input_file" ]; then
    echo "Usage: $0 <input_file>" >&2
    echo "Example: $0 data.txt" >&2
    echo "Output: Creates data_5prime.bed and data_3prime.bed" >&2
    exit 1
fi

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found!" >&2
    exit 1
fi

# Create output filenames
base_name=$(basename "$input_file" | sed 's/\.[^.]*$//')
output_5prime="${base_name}_5prime.bed"
output_3prime="${base_name}_3prime.bed"

echo "Processing $input_file..."
echo "Creating 5' positions BED file: $output_5prime"
echo "Creating 3' positions BED file: $output_3prime"

# Process 5' positions (columns 1, 2, 5)
awk 'BEGIN {FS="\t"; OFS="\t"} 
     NR > 1 {
         end5 = $2 + 3
         print "chr"$1, $2, end5, $5
     }' "$input_file" > "$output_5prime"

# Process 3' positions (columns 3, 4, 6)  
awk 'BEGIN {FS="\t"; OFS="\t"} 
     NR > 1 {
         end3 = $4 + 3
         print "chr"$3, $4, end3, $6
     }' "$input_file" > "$output_3prime"

echo "Done! Created:"
echo "  - $output_5prime (5' positions with gene names)"
echo "  - $output_3prime (3' positions with gene names)"

