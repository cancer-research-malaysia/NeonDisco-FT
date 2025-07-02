#!/usr/bin/env python3
"""
Script to update fusion breakpoint coordinates using hg38 lifted coordinates
from reference BED files.
"""

import pandas as pd
import argparse
import sys
from pathlib import Path

def load_bed_coordinates(bed_file):
    """
    Load BED file and return list of coordinates (second column).
    """
    try:
        # Read BED file (assuming tab-separated: chr, start, end, gene)
        bed_df = pd.read_csv(bed_file, sep='\t', header=None, 
                           names=['chromosome', 'start', 'end', 'gene'])
        
        # Return list of start positions (second column)
        return bed_df['start'].tolist()
    
    except Exception as e:
        print(f"Error loading BED file {bed_file}: {e}")
        return []

def update_coordinates(fusion_file, bed_5prime, bed_3prime, output_file):
    """
    Update fusion coordinates using paired BED files (row-by-row replacement).
    """
    try:
        # Load fusion data
        fusion_df = pd.read_csv(fusion_file, sep='\t')
        
        # Load coordinates from BED files
        print("Loading 5' BED coordinates...")
        coords_5prime = load_bed_coordinates(bed_5prime)
        print(f"Loaded {len(coords_5prime)} coordinates from 5' BED file")
        
        print("Loading 3' BED coordinates...")
        coords_3prime = load_bed_coordinates(bed_3prime)
        print(f"Loaded {len(coords_3prime)} coordinates from 3' BED file")
        
        # Check if row counts match
        fusion_rows = len(fusion_df)
        if len(coords_5prime) != fusion_rows:
            print(f"WARNING: 5' BED file has {len(coords_5prime)} rows, but fusion file has {fusion_rows} rows")
        if len(coords_3prime) != fusion_rows:
            print(f"WARNING: 3' BED file has {len(coords_3prime)} rows, but fusion file has {fusion_rows} rows")
        
        # Determine how many rows we can safely update
        max_rows = min(fusion_rows, len(coords_5prime), len(coords_3prime))
        if max_rows < fusion_rows:
            print(f"Will only update first {max_rows} rows due to BED file size mismatch")
        
        # Create a copy of the dataframe to modify
        updated_df = fusion_df.copy()
        
        # Direct row-by-row replacement
        print("Performing row-by-row coordinate replacement...")
        for i in range(max_rows):
            updated_df.at[i, "Position 5'"] = coords_5prime[i]
            updated_df.at[i, "Position 3'"] = coords_3prime[i]
        
        # Add new columns: breakpointID and fusionGenePair
        print("Adding breakpointID and fusionGenePair columns...")
        updated_df['breakpointID'] = (
            updated_df["Chromosome 5'"].astype(str) + ':' + 
            updated_df["Position 5'"].astype(str) + '-' +
            updated_df["Chromosome 3'"].astype(str) + ':' + 
            updated_df["Position 3'"].astype(str)
        )
        
        updated_df['fusionGenePair'] = (
            updated_df["Gene 5'"].astype(str) + '::' + 
            updated_df["Gene 3'"].astype(str)
        )
        
        # Save updated data
        updated_df.to_csv(output_file, sep='\t', index=False)
        
        # Print summary
        print(f"\nUpdate Summary:")
        print(f"Total fusion entries: {fusion_rows}")
        print(f"Rows updated: {max_rows}")
        print(f"5' coordinates replaced: {max_rows}")
        print(f"3' coordinates replaced: {max_rows}")
        print(f"Added breakpointID and fusionGenePair columns")
        
        print(f"\nUpdated data saved to: {output_file}")
        
        # Show a sample of the changes
        if max_rows > 0:
            print(f"\nSample of updates (first 3 rows):")
            print("Row | Gene 5' | New Pos 5' | Gene 3' | New Pos 3' | breakpointID | fusionGenePair")
            print("-" * 90)
            for i in range(min(3, max_rows)):
                new_5 = updated_df.at[i, "Position 5'"]
                new_3 = updated_df.at[i, "Position 3'"]
                gene_5 = updated_df.at[i, "Gene 5'"]
                gene_3 = updated_df.at[i, "Gene 3'"]
                breakpoint_id = updated_df.at[i, "breakpointID"]
                fusion_pair = updated_df.at[i, "fusionGenePair"]
                print(f"{i+1:3} | {gene_5:7} | {new_5:10} | {gene_3:7} | {new_3:10} | {breakpoint_id:20} | {fusion_pair}")
        
        # Show column info
        print(f"\nFinal output columns: {list(updated_df.columns)}")
        
    except Exception as e:
        print(f"Error processing files: {e}")
        return False
    
    return True

def main():
    parser = argparse.ArgumentParser(
        description="Update fusion breakpoint coordinates using paired BED files (row-by-row replacement)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python liftover_coordinates.py fusion_data.txt 5prime_hg38.bed 3prime_hg38.bed -o updated_fusion_data.txt
  
Input file format (tab-separated):
  Chromosome 5'	Position 5'	Chromosome 3'	Position 3'	Gene 5'	Gene 3'	Cell line
  
BED file format (tab-separated):
  chr19	49993179	49993182	RPL13A
        """
    )
    
    parser.add_argument('fusion_file', 
                       help='Input fusion data file (tab-separated)')
    parser.add_argument('bed_5prime', 
                       help="5' reference BED file with hg38 coordinates")
    parser.add_argument('bed_3prime', 
                       help="3' reference BED file with hg38 coordinates")
    parser.add_argument('-o', '--output', 
                       default='updated_coordinates_Klijn_fusion_data.tsv',
                       help='Output file name (default: updated_coordinates_Klijn_fusion_data.tsv)')
    
    args = parser.parse_args()
    
    # Check if input files exist
    for file_path in [args.fusion_file, args.bed_5prime, args.bed_3prime]:
        if not Path(file_path).exists():
            print(f"Error: File not found: {file_path}")
            sys.exit(1)
    
    # Run the coordinate update
    success = update_coordinates(args.fusion_file, args.bed_5prime, 
                               args.bed_3prime, args.output)
    
    if success:
        print("Coordinate update completed successfully!")
    else:
        print("Coordinate update failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()
    