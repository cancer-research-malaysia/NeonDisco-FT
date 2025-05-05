import marimo

__generated_with = "0.13.0"
app = marimo.App()


@app.cell
def _(mo):
    mo.md(r"### **Processing of Raw List of Gene Fusion Transcripts from Arriba and FusionCatcher on CCLE+internal Cell Lines**")
    return

@app.cell
def _(mo):
    mo.md(
        r"""
        This notebook details the processes (semi-automated) done to further process the raw output files from Arriba and FusionCatcher fusion transcript callers. 
    
        1. Define a function (adapted from the external Python executable called `combine-ft-nf.py`) inside the marimo notebook to process the raw output files from Arriba and FusionCatcher. This function will take the following parameters:
        	- `sample_id`: the sample ID (e.g. `373T`)
        	- `arr_file`: the Arriba raw output file
        	- `arr_tool_suffix`: the tool suffix matching the file name (e.g. `arr`)
        	- `fc_file`: the FusionCatcher raw output file
        	- `fc_tool_suffix`: the tool suffix matching the file name (e.g. `fc`)
        	- `output_dir`: the output directory to save the processed files
    
        2. Then, load up the two datasets on Jupyter Notebook and concatenate the dataframes together so that Arriba+FusionCatcher unfiltered FT data are combined into one data table and saved in one `.parquet` and `.tsv` file. Do the same for the `TCGANormals` panel of normals.
        """
    )
    return

@app.cell
def _():
    import os
    import re
    import polars as pl
    from pathlib import Path
    
    def extract_sample_id(filename, suffix):
        pattern = rf'^(\d+)T_{re.escape(suffix)}\.tsv$'
        match = re.search(pattern, os.path.basename(filename))
        return match.group(1) if match else None
    
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
                (pl.col('#gene1') + "::" + pl.col('gene2')).alias("fusionGenePair"),
                (pl.col('breakpoint1').str.replace("chr", "") + "-" + pl.col('breakpoint2').str.replace("chr", "")).alias("breakpointID"),
                (pl.col('strand1(gene/fusion)').str.split("/").list.get(1)).alias("strand1"),
                (pl.col('strand2(gene/fusion)').str.split("/").list.get(1)).alias("strand2"),
                'site1', 
                'site2', 
                'type', 
                'confidence',
                pl.lit(sample_id).alias("sampleID"),
                pl.lit(sample_id).cast(pl.Utf8).str.zfill(3).alias("sampleID_padded"),
                pl.lit(tool_name).alias("toolID")
                ])
            case 'FusionCatcher':
                # Handle NaN values in gene symbol columns by replacing with gene IDs
                gene1_expr = (
                    pl.when(pl.col('Gene_1_symbol(5end_fusion_partner)').is_null() | (pl.col('Gene_1_symbol(5end_fusion_partner)') == ""))
                    .then(pl.col('Gene_1_id(5end_fusion_partner)'))
                    .otherwise(pl.col('Gene_1_symbol(5end_fusion_partner)'))
                )
                
                gene2_expr = (
                    pl.when(pl.col('Gene_2_symbol(3end_fusion_partner)').is_null() | (pl.col('Gene_2_symbol(3end_fusion_partner)') == ""))
                    .then(pl.col('Gene_2_id(3end_fusion_partner)'))
                    .otherwise(pl.col('Gene_2_symbol(3end_fusion_partner)'))
                )
                
                base_columns = [
                (gene1_expr + "::" + gene2_expr + '__' + extract_fuscat_breakpoint(pl.col('Fusion_point_for_gene_1(5end_fusion_partner)')) + "-" + extract_fuscat_breakpoint(pl.col('Fusion_point_for_gene_2(3end_fusion_partner)'))).alias("fusionTranscriptID"),
                (gene1_expr + "::" + gene2_expr).alias("fusionGenePair"),
                (extract_fuscat_breakpoint(pl.col('Fusion_point_for_gene_1(5end_fusion_partner)')) + "-" + extract_fuscat_breakpoint(pl.col('Fusion_point_for_gene_2(3end_fusion_partner)'))).alias("breakpointID"),
                (pl.col('Fusion_point_for_gene_1(5end_fusion_partner)').str.split(":").list.get(2)).alias("strand1"),
                (pl.col('Fusion_point_for_gene_2(3end_fusion_partner)').str.split(":").list.get(2)).alias("strand2")
                ]
                # Check if 'Predicted_effect' column exists
                if 'Predicted_effect' in lazy_df.collect_schema().names():
                    predicted_effect_columns = [
                        pl.col('Predicted_effect').str.extract(r'^([^/]+)(?:/|$)').alias('site1'),
                        pl.when(pl.col('Predicted_effect').str.contains('/'))
                            .then(pl.col('Predicted_effect').str.extract(r'/(.+)$'))
                            .otherwise(pl.lit('.'))
                            .alias('site2')
                        ]
                else:
                    predicted_effect_columns = [
                        pl.lit('.').alias("site1"),
                        pl.lit('.').alias("site2")
                        ]
                
                return lazy_df.select(base_columns + predicted_effect_columns + [pl.lit('.').alias("type"),
                pl.lit('.').alias("confidence"),
                pl.lit(sample_id).alias("sampleID"),
                pl.lit(sample_id).cast(pl.Utf8).str.zfill(3).alias("sampleID_padded"),
                pl.lit(tool_name).alias("toolID")]) 
            case _:
                raise ValueError(f"Unsupported tool name: {tool_name}")
    
    def combine_fusion_data(sample_name, arr_file_path, fc_file_path, output_dir=None):
        """
        Combine fusion transcript data from Arriba and FusionCatcher into unified format.
        
        Parameters:
        - sample_name: Sample identifier
        - arr_file_path: Path to Arriba TSV file
        - fc_file_path: Path to FusionCatcher TSV file
        - output_dir: Directory to save output files (default: current directory)
        
        Returns:
        - Parquet and TSV dataframe of combined results
        """
        # Create tuples of filepath and tool suffix
        input_tuples = [(os.path.abspath(arr_file_path), 'arr'), (os.path.abspath(fc_file_path), 'fc')]
        
        # Print parameters for tracking
        print(f"Sample name: {sample_name}")
        print(f"Input files: {[path for path, _ in input_tuples]}")
        
        # Initialize empty list to store the lazy DataFrames
        lazy_dfs = []
        
        for input_path, suffix in input_tuples:
            if not Path(input_path).exists():
                raise FileNotFoundError(f"File {input_path} does not exist.")
            
            print('Setting tool name...')
            tool_name = {'arr': 'Arriba', 'fc': 'FusionCatcher'}.get(suffix)
            
            # Extract sample id
            sample_id = extract_sample_id(input_path, suffix) or sample_name
            print(f'Reading {tool_name} of {sample_name} TSV file...(sample ID: {sample_id})')
            
            # Create a lazy dataframe for each file
            lazy_df = wrangle_df(input_path, sample_id, tool_name)
            
            # Append the lazy DataFrame to the list
            lazy_dfs.append(lazy_df)
        
        # Concatenate all lazy DataFrames
        print("Concatenating lazy DataFrames from Arriba and FusionCatcher...")
        combined_lazy_df = pl.concat(lazy_dfs, rechunk=True)
        
        print("Concatenation completed. Collecting...")
        
        # Sort by Sample ID (padded), drop that column, then collect
        results = combined_lazy_df.sort("sampleID_padded").drop("sampleID_padded").with_columns(
            [pl.col(col).cast(pl.Categorical) for col in ['fusionTranscriptID', 'fusionGenePair', 'breakpointID', 'strand1', 'strand2', 'site1', 'site2', 'type', 'confidence', 'toolID']]).with_columns(
                pl.col("sampleID").cast(pl.Int64)).collect()
        
        # Save files if output directory is provided
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            parquet_path = output_dir / f"{sample_name}-combined-tool-FT-UNFILTERED.parquet"
            tsv_path = output_dir / f"{sample_name}-combined-tool-FT-UNFILTERED.tsv"
            
            print(f"Saving as parquet and tsv files to {output_dir}...")
            results.write_parquet(parquet_path)
            results.write_csv(tsv_path, separator="\t")
            
            print(f"Files saved to:\n- {parquet_path}\n- {tsv_path}")
        
        print("Done.")
        return results
    
    # Export the function to the global namespace
    return combine_fusion_data


@app.cell
def _(combine_fusion_data):
    # Example usage of the function with UI elements
    import marimo as mo
    
    # Create UI elements for the input parameters
    sample_name_input = mo.ui.text(value="373T", label="Sample Name")
    arr_file_input = mo.ui.text(
        value="../../data/CCLE+internal/CCLE/249-CCLE-ARR/A-253_arriba-output/arriba-fusions.tsv", 
        label="Arriba CCLE Test File Path"
    )
    fc_file_input = mo.ui.text(
        value="../../data/CCLE+internal/CCLE/249-CCLE-FC/A-253_fusioncatcher-output/final-list_candidate-fusion-genes.txt", 
        label="FusionCatcher CCLE Test File Path"
    )
    output_dir_input = mo.ui.text(
        value="../../output/CCLE+internal", 
        label="Output Directory"
    )
    
    # Display the inputs in a UI
    mo.vstack([
        mo.md("### Combine Fusion Transcript Data"),
        mo.md("Enter the parameters to combine Arriba and FusionCatcher data:"),
        sample_name_input,
        arr_file_input,
        fc_file_input,
        output_dir_input
    ])
    
    return arr_file_input, fc_file_input, output_dir_input, sample_name_input

@app.cell
def _(arr_file_input, combine_fusion_data, fc_file_input, mo, output_dir_input, sample_name_input):
    # Create a button to run the function
    run_button = mo.ui.button(label="Run Combination")
    
    # Display the button
    mo.hstack([run_button])
    
    # Return the button to be used in subsequent cells
    return run_button

@app.cell
def _(arr_file_input, combine_fusion_data, fc_file_input, mo, output_dir_input, run_button, sample_name_input):
    # Handle the button click
    if run_button.value:
        try:
            # Use the function with the input values
            results = combine_fusion_data(
                sample_name_input.value,
                arr_file_input.value,
                fc_file_input.value,
                output_dir_input.value
            )
            
            # Show a preview of the results
            output_display = mo.vstack([
                mo.md("### Results Preview"),
                mo.md(f"Combined {len(results)} fusion transcript records."),
                results.head(10)
            ])
        except Exception as e:
            mo.md(f"**Error:** {str(e)}")
    else:
        mo.md("Click the button to run the combination.")


# Direct function call without UI
@app.cell
def _(combine_fusion_data):
    # Example of direct call without UI
    # Uncomment to use this approach instead
    '''
    sample_name = "373T"
    arr_file = "../data/minimal-test/373T_arr.tsv"
    fc_file = "../data/minimal-test/373T_fc.tsv"
    output_dir = "../output/combined"
    
    results = combine_fusion_data(sample_name, arr_file, fc_file, output_dir)
    results.head()
    '''
    pass

if __name__ == "__main__":
    app.run()
