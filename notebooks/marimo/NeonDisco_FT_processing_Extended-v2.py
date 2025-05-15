

import marimo

__generated_with = "0.13.0"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## **Processing of Raw List of Gene Fusion Transcripts from Arriba and FusionCatcher on CCLE+internal Cell Lines**""")
    return


@app.cell(hide_code=True)
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

        2. Then, load up the two datasets on Marimo Notebook and concatenate the dataframes together so that Arriba+FusionCatcher unfiltered FT data are combined into one dataframe and saved in one `.parquet` and `.tsv` file.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""Below is the functionified version of the external Py executable. (*toggle hidden*)""")
    return


@app.cell(hide_code=True)
def _():
    import os
    import re
    import polars as pl
    from pathlib import Path


    def extract_fuscat_breakpoint(s):
        return s.str.split(":").list.slice(0, 2).list.join(":")


    def natural_sort_key(s):
        return [
            int(c) if c.isdigit() else c.lower() for c in re.split(r"(\d+)", s)
        ]


    def wrangle_df(file_path, sample_id, tool_name):
        lazy_df = pl.scan_csv(file_path, separator="\t")
        match tool_name:
            case "Arriba":
                return lazy_df.select(
                    [
                        (
                            pl.col("#gene1")
                            + "::"
                            + pl.col("gene2")
                            + "__"
                            + pl.col("breakpoint1").str.replace("chr", "")
                            + "-"
                            + pl.col("breakpoint2").str.replace("chr", "")
                        ).alias("fusionTranscriptID"),
                        (pl.col("#gene1") + "::" + pl.col("gene2")).alias(
                            "fusionGenePair"
                        ),
                        (
                            pl.col("breakpoint1").str.replace("chr", "")
                            + "-"
                            + pl.col("breakpoint2").str.replace("chr", "")
                        ).alias("breakpointID"),
                        (
                            pl.col("strand1(gene/fusion)")
                            .str.split("/")
                            .list.get(1)
                        ).alias("strand1"),
                        (
                            pl.col("strand2(gene/fusion)")
                            .str.split("/")
                            .list.get(1)
                        ).alias("strand2"),
                        "site1",
                        "site2",
                        "type",
                        "confidence",
                        pl.lit(sample_id).alias("sampleID"),
                        pl.lit(tool_name).alias("toolID"),
                    ]
                )
            case "FusionCatcher":
                # Handle NaN values in gene symbol columns by replacing with gene IDs
                gene1_expr = (
                    pl.when(
                        pl.col("Gene_1_symbol(5end_fusion_partner)").is_null()
                        | (pl.col("Gene_1_symbol(5end_fusion_partner)") == "")
                    )
                    .then(pl.col("Gene_1_id(5end_fusion_partner)"))
                    .otherwise(pl.col("Gene_1_symbol(5end_fusion_partner)"))
                )

                gene2_expr = (
                    pl.when(
                        pl.col("Gene_2_symbol(3end_fusion_partner)").is_null()
                        | (pl.col("Gene_2_symbol(3end_fusion_partner)") == "")
                    )
                    .then(pl.col("Gene_2_id(3end_fusion_partner)"))
                    .otherwise(pl.col("Gene_2_symbol(3end_fusion_partner)"))
                )

                base_columns = [
                    (
                        gene1_expr
                        + "::"
                        + gene2_expr
                        + "__"
                        + extract_fuscat_breakpoint(
                            pl.col("Fusion_point_for_gene_1(5end_fusion_partner)")
                        )
                        + "-"
                        + extract_fuscat_breakpoint(
                            pl.col("Fusion_point_for_gene_2(3end_fusion_partner)")
                        )
                    ).alias("fusionTranscriptID"),
                    (gene1_expr + "::" + gene2_expr).alias("fusionGenePair"),
                    (
                        extract_fuscat_breakpoint(
                            pl.col("Fusion_point_for_gene_1(5end_fusion_partner)")
                        )
                        + "-"
                        + extract_fuscat_breakpoint(
                            pl.col("Fusion_point_for_gene_2(3end_fusion_partner)")
                        )
                    ).alias("breakpointID"),
                    (
                        pl.col("Fusion_point_for_gene_1(5end_fusion_partner)")
                        .str.split(":")
                        .list.get(2)
                    ).alias("strand1"),
                    (
                        pl.col("Fusion_point_for_gene_2(3end_fusion_partner)")
                        .str.split(":")
                        .list.get(2)
                    ).alias("strand2"),
                ]
                # Check if 'Predicted_effect' column exists
                if "Predicted_effect" in lazy_df.collect_schema().names():
                    predicted_effect_columns = [
                        pl.col("Predicted_effect")
                        .str.extract(r"^([^/]+)(?:/|$)")
                        .alias("site1"),
                        pl.when(pl.col("Predicted_effect").str.contains("/"))
                        .then(pl.col("Predicted_effect").str.extract(r"/(.+)$"))
                        .otherwise(pl.lit("."))
                        .alias("site2"),
                    ]
                else:
                    predicted_effect_columns = [
                        pl.lit(".").alias("site1"),
                        pl.lit(".").alias("site2"),
                    ]

                return lazy_df.select(
                    base_columns
                    + predicted_effect_columns
                    + [
                        pl.lit(".").alias("type"),
                        pl.lit(".").alias("confidence"),
                        pl.lit(sample_id).alias("sampleID"),
                        pl.lit(tool_name).alias("toolID"),
                    ]
                )
            case _:
                raise ValueError(f"Unsupported tool name: {tool_name}")


    def combine_fusion_data(
        sample_name, arr_file_path, fc_file_path, output_dir=None
    ):
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
        input_tuples = [
            (os.path.abspath(arr_file_path), "arr"),
            (os.path.abspath(fc_file_path), "fc"),
        ]

        # Print parameters for tracking
        print(f"Sample name: {sample_name}")
        print(f"Input files: {[path for path, _ in input_tuples]}")

        # Initialize empty list to store the lazy DataFrames
        lazy_dfs = []

        for input_path, suffix in input_tuples:
            if not Path(input_path).exists():
                raise FileNotFoundError(f"File {input_path} does not exist.")

            print("Setting tool name...")
            tool_name = {"arr": "Arriba", "fc": "FusionCatcher"}.get(suffix)

            # Extract sample id
            sample_id = sample_name
            print(
                f"Reading {tool_name} of {sample_name} TSV file...(sample ID: {sample_id})"
            )

            # Create a lazy dataframe for each file
            lazy_df = wrangle_df(input_path, sample_id, tool_name)

            # Append the lazy DataFrame to the list
            lazy_dfs.append(lazy_df)

        # Concatenate all lazy DataFrames
        print("Concatenating lazy DataFrames from Arriba and FusionCatcher...")
        combined_lazy_df = pl.concat(lazy_dfs, rechunk=True)

        print("Concatenation completed. Collecting...")

        # collect
        results = (
            combined_lazy_df.with_columns(
                [
                    pl.col(col).cast(pl.Categorical)
                    for col in [
                        "fusionTranscriptID",
                        "fusionGenePair",
                        "breakpointID",
                        "strand1",
                        "strand2",
                        "site1",
                        "site2",
                        "type",
                        "confidence",
                        "toolID",
                    ]
                ]
            )
            .with_columns(pl.col("sampleID"))
            .collect()
        )

        # Save files if output directory is provided
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            parquet_path = (
                output_dir / f"{sample_name}-combined-tool-FT-UNFILTERED.parquet"
            )
            tsv_path = (
                output_dir / f"{sample_name}-combined-tool-FT-UNFILTERED.tsv"
            )

            print(f"Saving as parquet and tsv files to {output_dir}...")
            results.write_parquet(parquet_path)
            results.write_csv(tsv_path, separator="\t")

            print(f"Files saved to:\n- {parquet_path}\n- {tsv_path}")

        print("Done.")
        return results


    # Export the function to the global namespace
    return combine_fusion_data, pl


@app.cell(hide_code=True)
def _():
    # Define sample files to test using UI
    import marimo as mo

    sample_name_input = mo.ui.text(value='A-253', label="Sample Name")

    output_dir_input = mo.ui.text(value='/home/ec2-user/repos/FT-NeonDisco/output/CCLE+internal', label="Output Directory")

    # Define file paths as functions of the sample name
    def get_arriba_path(sample_name):
        return f"/home/ec2-user/repos/FT-NeonDisco/data/CCLE+internal/CCLE/249-CCLE-ARR/{sample_name}_arriba-output/arriba-fusions.tsv"

    def get_fusioncatcher_path(sample_name):
        return f"/home/ec2-user/repos/FT-NeonDisco/data/CCLE+internal/CCLE/249-CCLE-FC/{sample_name}_fusioncatcher-output/final-list_candidate-fusion-genes.txt"

    return (
        get_arriba_path,
        get_fusioncatcher_path,
        mo,
        output_dir_input,
        sample_name_input,
    )


@app.cell(hide_code=True)
def _(
    get_arriba_path,
    get_fusioncatcher_path,
    mo,
    output_dir_input,
    sample_name_input,
):

    # Create reactive UI elements that update when sample_name changes
    arr_file_input = mo.ui.text(
        value=get_arriba_path(sample_name_input.value),
        label="Arriba CCLE Test File Path"
    )

    fc_file_input = mo.ui.text(
        value=get_fusioncatcher_path(sample_name_input.value),
        label="FusionCatcher CCLE Test File Path"
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

    return arr_file_input, fc_file_input


@app.cell
def _(
    arr_file_input,
    combine_fusion_data,
    fc_file_input,
    output_dir_input,
    sample_name_input,
):
    # Direct call without UI

    results = combine_fusion_data(sample_name_input.value, arr_file_input.value, fc_file_input.value, output_dir_input.value)

    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Collating Combined FT List into Cohort-wide Dataset

        We can then use a simple loop based on a manifest file of all the 249 chosen CCLE datasets + 16 internal cell line (cancer) datasets to **(1)** generate individual combined FT file, and then **(2)** concatenate all of them into one cohort-wide dataframe.
        ## STEP 1
        """
    )
    return


@app.cell
def _(combine_fusion_data):
    def _():
        # define manifest file of CCLE+internal
        import os

        combined_manifest_f = os.path.abspath('/home/ec2-user/repos/FT-NeonDisco/data/CCLE+internal/CCLE+internalCL-processed.IDlist.txt')

        internal_f = os.path.abspath('/home/ec2-user/repos/FT-NeonDisco/data/CCLE+internal/16-inhouse-cell-lines-processed.IDlist.txt')
        with open(internal_f) as f:
            internal_set = set(line.strip() for line in f)
        print(internal_set)


        # read the combined manifest file
        with open(combined_manifest_f) as f:
            for line in f:
                # get the sample ID
                sample_idx = line.strip()
                if sample_idx not in internal_set:
                    # get the Arriba and FusionCatcher file paths
                    arr_fpath = f"/home/ec2-user/repos/FT-NeonDisco/data/CCLE+internal/CCLE/249-CCLE-ARR/{sample_idx}_arriba-output/arriba-fusions.tsv"
                    fc_fpath = f"/home/ec2-user/repos/FT-NeonDisco/data/CCLE+internal/CCLE/249-CCLE-FC/{sample_idx}_fusioncatcher-output/final-list_candidate-fusion-genes.txt"
                else:
                    arr_fpath = f"/home/ec2-user/repos/FT-NeonDisco/data/CCLE+internal/IN-HOUSE/16-inhouse-ARR/{sample_idx}_arriba-output/arriba-fusions.tsv"
                    fc_fpath = f"/home/ec2-user/repos/FT-NeonDisco/data/CCLE+internal/IN-HOUSE/16-inhouse-FC/{sample_idx}_fusioncatcher-output/final-list_candidate-fusion-genes.txt"

                # call the function
                combine_fusion_data(sample_idx, arr_fpath, fc_fpath, output_dir='/home/ec2-user/repos/FT-NeonDisco/output/CCLE+internal')
    _()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## STEP 2""")
    return


@app.cell
def _():
    # read the output directory path, and for each file, read the parquet file into a Polars dataframe, then concatenate all of them into one dataframe
    def _():
        import os
        import polars as pl
        from pathlib import Path

        # Define the output directory
        output_dir = Path('/home/ec2-user/repos/FT-NeonDisco/output/CCLE+internal')
        # Get a list of all Parquet files in the output directory
        parquet_files = list(output_dir.glob("*.parquet"))
        # Initialize an empty list to store DataFrames
        dfs = []
        # Loop through each Parquet file and read it into a Polars DataFrame
        for file in parquet_files:
            # Read the Parquet file into a Polars DataFrame
            df = pl.read_parquet(file)
            # Append the DataFrame to the list
            dfs.append(df)
        # Concatenate all DataFrames in the list into a single DataFrame
        combined_df = pl.concat(dfs, rechunk=True)
        # Return the combined DataFrame
        return combined_df
    # Call the function to get the combined DataFrame
    combined_df = _()
    # Display the first few rows of the combined DataFrame

    # save the combined DataFrame to a Parquet file
    combined_df.write_parquet('/home/ec2-user/repos/FT-NeonDisco/output/CCLE+internal/01-CCLE+internal-ALL-FT-UNFILTERED.parquet')
    return


@app.cell
def _():
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Filtering The Combined Cohort-wide FT dataset

        We can first check the combined cohort-wide FT dataframe from CCLE+internal cell lines.
        """
    )
    return


@app.cell
def _(pl):
    # load up the concatenated dataframe

    concat_df = pl.read_parquet('/home/ec2-user/repos/FT-NeonDisco/output/CCLE+internal/01-CCLE+internal-ALL-FT-UNFILTERED.parquet')

    return (concat_df,)


@app.cell
def _(concat_df, pl):
    df = (
            concat_df
                .group_by("breakpointID")
                .agg(pl.len(), 
                    pl.col("sampleID").unique().alias("sampleIDs"),
                    pl.col("toolID")
                    )
                .sort("len", descending=True)
         )


    return


@app.cell
def _(mo):
    mo.md(r"""We need to get a unique `breakpointID` column BUT also maintain the `sampleID` information. Maybe try aggregating by `sampleID` then do unique on `breakpointID` then return just these two columns to avoid column collision.""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Let's convert the Polars dataframe into a Pandas one for easy processing.""")
    return


@app.cell
def _():
    # now use TCGA-Normals to remove non-tumor-specific breakpoints
    import pandas as pd

    panel_of_norms_df = pd.read_csv('/home/ec2-user/repos/FT-NeonDisco/output/TCGANormals/Arr-and-FC_TCGANormals-UNIQUE-breakpointID-list.tsv', sep='\t')

    return (panel_of_norms_df,)


@app.cell
def _(panel_of_norms_df):
    # then get the unique breakpoints as set

    panel_of_normset = set(panel_of_norms_df['breakpointID'])

    return (panel_of_normset,)


@app.cell
def _(panel_of_normset):
    # test
    bp = '2:62263691-22:23745219'

    print('YAS!' if bp in panel_of_normset else 'NAUR...')
    return


if __name__ == "__main__":
    app.run()
