#!/bin/bash

# Check if correct number of arguments provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 <5prime.bed> <3prime.bed> <output.txt>"
    echo "Example: $0 5prime_breakpoints.bed 3prime_breakpoints.bed fusion_output.txt"
    exit 1
fi

# Assign arguments to variables
FIVEPRIME_BED="$1"
THREEPRIME_BED="$2"
OUTPUT_FILE="$3"

# Check if input files exist
if [ ! -f "$FIVEPRIME_BED" ]; then
    echo "Error: 5' BED file '$FIVEPRIME_BED' not found!"
    exit 1
fi

if [ ! -f "$THREEPRIME_BED" ]; then
    echo "Error: 3' BED file '$THREEPRIME_BED' not found!"
    exit 1
fi

# Use awk to combine the files
awk '
BEGIN {
    # Print header
    print "chrUp\tbreakpointUp\tchrDw\tbreakpointDw\tfusionGenePair\tbreakpointID"
    
    # Read all lines from the 3prime file into an array
    while ((getline line < ARGV[2]) > 0) {
        threeprime_lines[++threeprime_count] = line
    }
    close(ARGV[2])
    
    # Reset line counter
    line_num = 0
}
{
    # Increment line counter
    line_num++
    
    # Check if we have a corresponding 3prime line
    if (line_num > threeprime_count) {
        print "Error: 5prime file has more lines than 3prime file" > "/dev/stderr"
        exit 1
    }
    
    # Split the corresponding 3prime line
    split(threeprime_lines[line_num], threeprime, "\t")
    
    # Extract fields from 5prime (current line) and 3prime
    chr5 = substr($1, 4)      # Remove "chr" prefix
    pos5 = $2                 # 5prime position
    gene5 = $4                # 5prime gene
    
    chr3 = substr(threeprime[1], 4)  # Remove "chr" prefix from 3prime
    pos3 = threeprime[2]             # 3prime position  
    gene3 = threeprime[4]            # 3prime gene
    
    # Create fusionGenePair and breakpointID
    fusionGenePair = gene5 "::" gene3
    breakpointID = chr5 ":" pos5 "-" chr3 ":" pos3
    
    # Print the combined row
    print chr5 "\t" pos5 "\t" chr3 "\t" pos3 "\t" fusionGenePair "\t" breakpointID
}
END {
    # Check if files had same number of lines
    if (line_num < threeprime_count) {
        print "Warning: 3prime file has more lines than 5prime file" > "/dev/stderr"
    }
}
' "$FIVEPRIME_BED" "$THREEPRIME_BED" > "$OUTPUT_FILE"

echo "Fusion data combined successfully!"
echo "Input files:"
echo "  5' BED: $FIVEPRIME_BED ($(wc -l < "$FIVEPRIME_BED") lines)"
echo "  3' BED: $THREEPRIME_BED ($(wc -l < "$THREEPRIME_BED") lines)"
echo "Output: $OUTPUT_FILE ($(wc -l < "$OUTPUT_FILE") lines including header)"
echo ""
echo "First few lines of output:"
head -n 4 "$OUTPUT_FILE"

