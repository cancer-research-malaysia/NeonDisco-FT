#!/usr/bin/env python3

import os
import re
import sys
import polars as pl
from pathlib import Path

### uncomment this function for TCGA raw FT files
def extract_sample_id(filename, suffix):
    pattern = rf'^(TCGA-\w{{2}}-\w{{4}})_{re.escape(suffix)}\.tsv$'
    match = re.search(pattern, os.path.basename(filename))
    return match.group(1) if match else None

def list_files(path, suffix):
    return {extract_sample_id(str(file), suffix): str(file) 
            for file in Path(path).glob('*.tsv') 
            if file.is_file() and extract_sample_id(str(file), suffix)}

def extract_fuscat_breakpoint(s):
    return s.str.split(':').list.slice(0, 2).list.join(':')

def natural_sort_key(s):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', s)]

def wrangle_df(file_path, sample_id, tool_name):
    lazy_df = pl.scan_csv(file_path, separator="\t")
    match tool_name:
        case 'Arriba':
            return lazy_df.select([
            (pl.col('#gene1') + "::" + pl.col('gene2') + '__' + pl.col('breakpoint1').str.replace("chr", "") + "-" + pl.col('breakpoint2').str.replace("chr", "")).alias("fusionTranscriptID"),
            (pl.col('#gene1') + "::" + pl.col('gene2')).alias("fusionGeneID"),
            (pl.col('breakpoint1').str.replace("chr", "") + "-" + pl.col('breakpoint2').str.replace("chr", "")).alias("breakpointPair"),
            (pl.col('strand1(gene/fusion)').str.split("/").list.get(1)).alias("strand1"),
            (pl.col('strand2(gene/fusion)').str.split("/").list.get(1)).alias("strand2"),
            'site1', 
            'site2', 
            'type', 
            'confidence',
            pl.lit(sample_id).cast(pl.Utf8).alias("sampleID"),
            pl.lit(tool_name).alias("toolID")
            ])
        case 'FusionCatcher':
            base_columns = [
            (pl.col('Gene_1_symbol(5end_fusion_partner)') + "::" + pl.col('Gene_2_symbol(3end_fusion_partner)') + '__' + extract_fuscat_breakpoint(pl.col('Fusion_point_for_gene_1(5end_fusion_partner)')) + "-" + extract_fuscat_breakpoint(pl.col('Fusion_point_for_gene_2(3end_fusion_partner)'))).alias("fusionTranscriptID"),
            (pl.col('Gene_1_symbol(5end_fusion_partner)') + "::" + pl.col('Gene_2_symbol(3end_fusion_partner)')).alias("fusionGeneID"),
            (extract_fuscat_breakpoint(pl.col('Fusion_point_for_gene_1(5end_fusion_partner)')) + "-" + extract_fuscat_breakpoint(pl.col('Fusion_point_for_gene_2(3end_fusion_partner)'))).alias("breakpointPair"),
            (pl.col('Fusion_point_for_gene_1(5end_fusion_partner)').str.split(":").list.get(2)).alias("strand1"),
            (pl.col('Fusion_point_for_gene_2(3end_fusion_partner)').str.split(":").list.get(2)).alias("strand2")
            ]
            # Check if 'Predicted_effect' column exists
            if 'Predicted_effect' in lazy_df.collect_schema().names():
                # print(f'Sample ID: {sample_id} has Predicted_effect column')
                predicted_effect_columns = [
                    pl.col('Predicted_effect').str.extract(r'^([^/]+)(?:/|$)').alias('site1'),
                    pl.when(pl.col('Predicted_effect').str.contains('/'))
                        .then(pl.col('Predicted_effect').str.extract(r'/(.+)$'))
                        .otherwise(pl.lit('.'))
                        .alias('site2')
                    ]
            else:
                # print(f'Sample ID: {sample_id} does not have Predicted_effect column')
                predicted_effect_columns = [
                    pl.lit('.').alias("site1"),
                    pl.lit('.').alias("site2")
                    ]
            
            return lazy_df.select(base_columns + predicted_effect_columns + [pl.lit('.').alias("type"),
            pl.lit('.').alias("confidence"),
            pl.lit(sample_id).cast(pl.Utf8).alias("sampleID"),
            pl.lit(tool_name).alias("toolID")]) 
        case _:
            raise ValueError(f"Unsupported tool name: {tool_name}")

def main():
    input_path = os.path.abspath(sys.argv[1])
    file_suffix = sys.argv[2]

	# Check if the path exists
    if not os.path.exists(input_path):
        print(f"Error: The path '{input_path}' does not exist.")
        sys.exit(1)

    all_files = list_files(input_path, file_suffix)
    print(f'{len(all_files)} files in total.')

    print('Setting tool name...')
    
    tool_name = {'arr': 'Arriba', 'fc': 'FusionCatcher'}.get(file_suffix)
    if not tool_name:
        print(f'Error: Unknown prefix for tool name {file_suffix}. Aborting...')
        sys.exit(1)

    print(f'Reading {tool_name} TSV files by creating a list of lazy Frames...')

    ###### TESTING BLOCK ##########
    # print(all_files)
    # for i, (sample_id, file_path) in enumerate(all_files.items()):
    #     if i < 10:
    #         df = wrangle_df(file_path, sample_id, tool_name)
    #         print(df.collect())
    
    # ###### BATCH MODE BLOCK ##########
    # Create a list of lazy DataFrames
    lazy_dfs = [wrangle_df(file_path, sample_id, tool_name) for sample_id, file_path in all_files.items()]
    print(f"List of lazy Frames for all {len(all_files)} files has been created. Concatenating...")
    # Concatenate all lazy DataFrames
    combined_lazy_df = pl.concat(lazy_dfs, rechunk=True)
    print("Concatenation completed. Collecting...")

    # Sort by Sample ID, drop that column, then collect
    results = combined_lazy_df.sort("sampleID").with_columns(
        [pl.col(col).cast(pl.Categorical) for col in ['fusionTranscriptID', 'fusionGeneID', 'breakpointPair', 'strand1', 'strand2', 'site1', 'site2', 'type', 'confidence', 'toolID']]).with_columns(
            pl.col("sampleID")).collect()
    print(results)

    # save as parquet and tsv
    print(f"Saving as parquet and tsv files...")
    results.write_parquet(f"output/TCGANormals/{tool_name}-TCGANormals-fusiontranscript-raw-list.parquet")
    results.write_csv(f"output/TCGANormals/{tool_name}-TCGANormals-fusiontranscript-raw-list.tsv", separator="\t")

    print("Done.")

if __name__ == "__main__":
    main()