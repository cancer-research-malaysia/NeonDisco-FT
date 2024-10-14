#!/usr/bin/env python3

import os
import re
import sys
import polars as pl
from pathlib import Path

def extract_sample_id(filename, prefix):
    pattern = rf'{re.escape(prefix)}_(\d+)\.tsv$'
    match = re.search(pattern, os.path.basename(filename))
    return match.group(1) if match else None

def list_files(path, prefix):
    return {extract_sample_id(str(file), prefix): str(file) 
            for file in Path(path).rglob('*.tsv') 
            if file.is_file() and extract_sample_id(str(file), prefix)}

def wrangle_df(file_path, sample_id, tool_name):
    lazy_df = pl.scan_csv(file_path, separator="\t")
    lazy_df = lazy_df.select([
            (pl.col('#gene1') + "::" + pl.col('gene2') + '__' + pl.col('breakpoint1').str.replace("chr", "") + "-" + pl.col('breakpoint2').str.replace("chr", "")).alias("fusionTranscriptID"),
            (pl.col('#gene1') + "::" + pl.col('gene2')).alias("fusionGeneID"),
            (pl.col('breakpoint1').str.replace("chr", "") + "-" + pl.col('breakpoint2').str.replace("chr", "")).alias("breakpointPair"),
            (pl.col('strand1(gene/fusion)').str.split("/").list.get(1)).alias("strand1"),
            (pl.col('strand2(gene/fusion)').str.split("/").list.get(1)).alias("strand2"),
            'site1', 
            'site2', 
            'type', 
            'confidence',
            pl.lit(sample_id).alias("sampleID"),
            pl.lit(tool_name).alias("toolID")
        ])
    return lazy_df

def main():
    input_path = os.path.abspath(sys.argv[1])
    prefix = sys.argv[2]

	# Check if the path exists
    if not os.path.exists(input_path):
        print(f"Error: The path '{input_path}' does not exist.")
        sys.exit(1)

    all_files = list_files(input_path, prefix)
    print(f'{len(all_files)} files in total.')

    print('Setting tool name...')
    
    tool_name = {'arr': 'Arriba', 'fc': 'FusionCatcher'}.get(prefix)
    if not tool_name:
        print(f'Error: Unknown prefix for tool name {prefix}. Aborting...')
        sys.exit(1)

    print(f'Reading {tool_name} TSV files by creating a list of lazy Frames...')

    # for i, (sample_id, file_path) in enumerate(all_files.items()):
    #     if i < 5:
    #         df = wrangle_df(file_path, sample_id, tool_name)
    #         print(df.collect())

    # Create a list of lazy DataFrames
    lazy_dfs = [wrangle_df(file_path, sample_id, tool_name) 
            for sample_id, file_path in all_files.items()]
    print(f"List of lazy Frames for all {len(all_files)} files has been created. Concatenating...")
    # Concatenate all lazy DataFrames
    combined_lazy_df = pl.concat(lazy_dfs)
    print("Concatenation completed. Printing...")
    # Sort by Sample ID then collect 
    results = combined_lazy_df.sort("sampleID").collect()
    print(results)
    # save as parquet and tsv
    print(f"Saving as parquet and tsv files...")
    results.write_parquet(f"data/{tool_name}-fusiontranscript-raw-list.parquet")
    results.write_csv(f"data/{tool_name}-fusiontranscript-raw-list.tsv", separator="\t")

    print("Done.")

if __name__ == "__main__":
    main()