#!/usr/bin/env python3

import os
import sys
import polars as pl

def check_file_empty(file_path, file_description):
    """
    Check if a parquet file is empty or contains no data rows.
    Returns True if empty, False if it contains data.
    """
    try:
        df = pl.scan_parquet(file_path).collect()
        if df.height == 0:
            print(f"WARNING: {file_description} is empty (no data rows)")
            return True
        return False
    except Exception as e:
        print(f"ERROR: Could not read {file_description}: {e}")
        return True

def consolidate_duplicate_rows(df, groupby_cols):
    """
    Consolidate rows that have the same values in groupby_cols by filling
    NA values with non-NA values from other rows in the same group.
    
    Args:
        df: Polars DataFrame
        groupby_cols: List of column names to group by
    
    Returns:
        Polars DataFrame with consolidated rows
    """
    
    # Handle empty DataFrame case
    if df.height == 0:
        print("Input DataFrame is empty (no data rows). Returning empty DataFrame with same schema.")
        return df
    
    # Check if groupby columns exist
    missing_cols = [col for col in groupby_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Groupby columns not found in DataFrame: {missing_cols}")
    
    # Define aggregation strategy for non-groupby columns only
    # (groupby columns are automatically included in the result)
    agg_exprs = []
    
    for col in df.columns:
        if col not in groupby_cols:
            # For non-groupby columns, take the first non-null value
            # Use drop_nulls().first() which is more reliable than filter + first
            agg_exprs.append(
                pl.col(col).drop_nulls().first().alias(col)
            )
    
    # Group by the specified columns and aggregate
    result = df.group_by(groupby_cols).agg(agg_exprs)
    
    return result

def create_empty_output_files(output_filename):
    """
    Create empty output files with proper headers when no fusions are found.
    """
    empty_df = pl.DataFrame({
        'fusionTranscriptID': [],
        'fusionGenePair': [],
        'breakpointID': [],
        '5pStrand': [],
        '3pStrand': [],
        'originalTool': [],
        'sampleID': [],
        'sampleNum': [],
        'sampleNum_Padded': [],
        '5pSite_ARR': [],
        '3pSite_ARR': [],
        'mutationType_ARR': [],
        'confidenceLabel_ARR': [],
        'predictedEffect_FC': [],
        'fusionPairAnnotation_FC': [],
        'fusionPairAnnotation_SF': [],
        'junctionReadCount_SF': [],
        'spanningFragCount_SF': [],
        'largeAnchorSupport_SF': [],
        'detectedBy': [],
        'toolOverlapCount': [],
        'foundInCCLE&InternalCLs': [],
        'fusionGenePair_FusIns': []
    })
    
    # Write empty TSV
    empty_df.write_csv(f"{output_filename}.tsv", separator='\t')
    print(f"Empty results file saved to {output_filename}.tsv")
    
    # Write empty txt file for FusionInspector
    with open(f"{output_filename}-unique-genePairs-for-FusIns.txt", 'w') as f:
        pass  # Create empty file
    print(f"Empty fusion gene pairs file saved to {output_filename}-unique-genePairs-for-FusIns.txt")

def main():
    """
    Process and filter fusion transcript data based on:
    1. Consolidating duplicate entries (filling NA values across rows with same fusion IDs)
    2. Removing duplicate entries
    3. Adding tool overlap information
    4. Checking presence in CCLE and internal cell lines
    5. Filtering out fusions found in panel of normals
    6. Formatting for FusionInspector compatibility
    The script takes command-line arguments for input files and outputs the results in TSV format and a txt format for Fusion Inspector.
    """
    # Parse command-line arguments
    sample_name = sys.argv[1]
    parquet_input_file = os.path.abspath(sys.argv[2])
    panel_of_normals_file = os.path.abspath(sys.argv[3])
    ccle_internal_cell_line_file = os.path.abspath(sys.argv[4])
    output_filename = os.path.abspath(sys.argv[5])
    
    # Print parameters for debugging
    print(f"Sample name: {sample_name}")
    print(f"Parquet file of collated FTs: {parquet_input_file}")
    print(f"Parquet file of panel of Normals (TCGA Normals) FTs: {panel_of_normals_file}")
    print(f"Parquet file of CCLE + internal cell line FTs: {ccle_internal_cell_line_file}")
    print(f"Output filename: {output_filename}")

    # Check if main input file is empty - if so, create empty outputs and exit
    if check_file_empty(parquet_input_file, "Collated fusion transcript data"):
        print(f"No fusion transcripts detected for sample {sample_name}. Creating empty output files...")
        create_empty_output_files(output_filename)
        print("Processing complete - no fusions detected!")
        return
    
    # Load the collated fusion transcript data
    print("Loading collated fusion transcript data...")
    collated_df = pl.scan_parquet(parquet_input_file).collect()
    
    # Step 1: Consolidate duplicate rows by filling NA values across rows with the same fusion identifiers
    print("Consolidating duplicate rows (filling NA values)...")
    groupby_cols = ['fusionTranscriptID', 'fusionGenePair', 'breakpointID', '5pStrand', '3pStrand', 'originalTool']
    
    # Check if required columns exist for consolidation
    missing_cols = [col for col in groupby_cols if col not in collated_df.columns]
    if missing_cols:
        print(f"Warning: Expected consolidation columns not found: {missing_cols}")
        print(f"Available columns: {collated_df.columns}")
        print("Skipping consolidation step...")
        consolidated_df = collated_df
    else:
        consolidated_df = consolidate_duplicate_rows(collated_df, groupby_cols)
        print(f"Consolidation complete: {len(collated_df)} -> {len(consolidated_df)} rows")
    
    # Step 2: Filter for unique rows based on fusionTranscriptID and originalTool
    print("Filtering for unique rows based on fusionTranscriptID and originalTool...")
    unique_collated_df = consolidated_df.unique(subset=["fusionTranscriptID", "originalTool"])
    
    # Step 3: Create a new dataframe with unique fusionTranscriptIDs and list of tools that detected them
    print("Creating tool overlap information...")
    tool_group_df = (
        unique_collated_df
        .group_by('fusionTranscriptID')
        .agg(
            pl.col('originalTool').unique().alias('detectedBy'),
            pl.col('originalTool').unique().count().alias('toolOverlapCount')
        )
    )
    
    # Step 4: Create tool overlap information and combine data from multiple tools
    print("Creating tool overlap information and combining multi-tool data...")
    
    # First, create tool overlap summary
    tool_group_df = (
        unique_collated_df
        .group_by('fusionTranscriptID')
        .agg(
            pl.col('originalTool').unique().alias('detectedBy'),
            pl.col('originalTool').unique().count().alias('toolOverlapCount')
        )
    )
    
    # Then, consolidate data across tools for the same fusion
    # For each fusion, combine non-null values from all tools
    non_tool_cols = [col for col in unique_collated_df.columns 
                     if col not in ['fusionTranscriptID', 'originalTool']]
    
    fusion_agg_exprs = []
    for col in non_tool_cols:
        # Take the first non-null value across all tools for this fusion
        fusion_agg_exprs.append(
            pl.col(col).drop_nulls().first().alias(col)
        )
    
    # Combine data across tools for each fusion
    combined_fusion_df = (
        unique_collated_df
        .group_by('fusionTranscriptID')
        .agg(fusion_agg_exprs)
    )
    
    # Join with tool overlap information
    unique_fusions_df = combined_fusion_df.join(tool_group_df, on='fusionTranscriptID')
    
    # Step 5: Load CCLE & Internal Cell Line FT data
    print("Loading CCLE & Internal Cell Line FT data...")
    if check_file_empty(ccle_internal_cell_line_file, "CCLE & Internal Cell Line data"):
        ccle_set = set()  # Empty set if file is empty
    else:
        ccle_df = pl.scan_parquet(ccle_internal_cell_line_file).collect()
        ccle_set = set(ccle_df['breakpointID'].to_list())

    # Step 6: Add 'foundInCCLE&InternalCLs' column to unique fusions
    print("Adding 'foundInCCLE&InternalCLs' column to unique fusions...")
    ccle_added_df = unique_fusions_df.with_columns(
        pl.when(pl.col('breakpointID').is_in(ccle_set)).then(True).otherwise(False).alias('foundInCCLE&InternalCLs')
    )
    
    # Step 7: Load Panel of Normals (TCGA Normals) data
    print("Loading Panel of Normals data...")
    if check_file_empty(panel_of_normals_file, "Panel of Normals data"):
        pon_set = set()  # Empty set if file is empty
    else:
        pon_df = pl.scan_parquet(panel_of_normals_file).collect()
        pon_set = set(pon_df['breakpointID'].to_list())

    # Step 8: Filter out breakpoints that appear in the Panel of Normals
    print("Filtering out breakpoints that appear in the Panel of Normals...")
    normfilt_df = ccle_added_df.filter(~pl.col('breakpointID').is_in(pon_set))

    # Check if filtering removed all fusions
    if normfilt_df.height == 0:
        print("WARNING: All fusions were filtered out by Panel of Normals filtering.")
        print("Creating empty output files...")
        create_empty_output_files(output_filename)
        print("Processing complete - all fusions filtered out!")
        return

    # Step 9: Create FusionInspector format column
    print("Creating FusionInspector format column...")
    final_result_df = normfilt_df.with_columns(
        pl.col('fusionGenePair').cast(pl.Utf8).str.replace('::', '--').alias('fusionGenePair_FusIns')
    )
    
    # Step 10: Format detectedBy column and apply consensus filtering
    print("Formatting output and applying consensus filtering...")
    export_df = final_result_df.with_columns([
        pl.col('detectedBy').list.eval(pl.element().cast(pl.Utf8)).list.join(" | ").alias('detectedBy')
    ])

    # Filter for rows based on 'toolOverlapCount'
    # > 0 is default to keep the union of all fusions detected by at least one tool
    # change this value for different filtering
    export_consensus_df = export_df.filter(pl.col('toolOverlapCount') > 0)

    # Check if consensus filtering removed all fusions
    if export_consensus_df.height == 0:
        print("WARNING: No fusions passed the consensus filter.")
        print("Creating empty output files...")
        create_empty_output_files(output_filename)
        print("Processing complete - no consensus fusions found!")
        return

    # Step 11: Save results
    print(f"Saving filtered results to {output_filename}.tsv...")
    export_consensus_df.write_csv(f"{output_filename}.tsv", separator='\t')
    print(f"Results saved to {output_filename}.tsv")
    
    # Write unique fusion gene pairs for FusionInspector
    export_consensus_df.select('fusionGenePair_FusIns').unique().write_csv(f"{output_filename}-unique-genePairs-for-FusIns.txt", include_header=False)
    print(f"Unique fusion gene pairs for FusionInspector saved to {output_filename}-unique-genePairs-for-FusIns.txt")
    
    print("Processing complete!")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: wrangle-and-filter-FTs--nf.py <sample id> <parquet file of combined FTs> <panel of normals FTs parquet file> <ccle+internal cell line FTs parquet file> <output filename>")
        sys.exit(1)
    main()
