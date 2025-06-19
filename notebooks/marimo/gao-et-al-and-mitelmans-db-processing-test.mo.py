import marimo

__generated_with = "0.13.15"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    return mo, pl


@app.cell
def _(pl):
    # load tsv into polars
    gao_df = pl.read_csv("/home/ec2-user/repos/FT-NeonDisco/data/Gao-et-al-TCGA-fusion-set/Gao-et-al-TCGA-fusion-set.tsv.txt", separator='\t')
    gao_df
    return (gao_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Rename The Columns

    First let's rename all the columns in the raw Gao et al. extracted fusion set.
    """
    )
    return


@app.cell
def _(gao_df):
    gao_pdf = gao_df.rename({"Breakpoint1":"breakpoint1", "Breakpoint2":"breakpoint2", "Fusion":"fusionGenePair", "Junction":"junctionReads", "Spanning":"spanningReads", "Cancer":"TCGACancerCode", "Sample":"TCGASampleID" })
    gao_pdf
    return (gao_pdf,)


@app.cell
def _(mo):
    mo.md(
        r"""
    ### Create `fusionTranscriptID`

    Then create a new column called `fusionTranscriptID` by combining breakpoint information with gene name pair. Then combine breakpoint columns into one, and separate strandedness into columns.
    """
    )
    return


@app.cell
def _(gao_pdf, pl):
    gao_pdf_recol = gao_pdf.with_columns([
        # Extract chromosome (remove 'chr' prefix)
        pl.col("breakpoint1").str.extract(r"chr(\w+):", 1).alias("chr1"),
        pl.col("breakpoint2").str.extract(r"chr(\w+):", 1).alias("chr2"),

        # Extract position
        pl.col("breakpoint1").str.extract(r":(\d+):", 1).alias("pos1"),
        pl.col("breakpoint2").str.extract(r":(\d+):", 1).alias("pos2"),

        # Extract strand
        pl.col("breakpoint1").str.extract(r"([+-])$", 1).alias("5pStrand"),
        pl.col("breakpoint2").str.extract(r"([+-])$", 1).alias("3pStrand"),
    ]).select([
        (pl.col("fusionGenePair").str.replace("--", "::") + "__" + 
         pl.concat_str([pl.col("chr1"), pl.col("pos1")], separator=":") + "-" + 
         pl.concat_str([pl.col("chr2"), pl.col("pos2")], separator=":")).alias("fusionTranscriptID"),
        pl.col("fusionGenePair").str.replace("--", "::"),
        (pl.concat_str([pl.col("chr1"), pl.col("pos1")], separator=":") + "-" + 
         pl.concat_str([pl.col("chr2"), pl.col("pos2")], separator=":")).alias("breakpointID"),
        pl.col("5pStrand"),
        pl.col("3pStrand"),
        pl.col("junctionReads"),
        pl.col("spanningReads"),
        pl.col("TCGACancerCode"),
        pl.col("TCGASampleID")
    ])

    gao_pdf_recol
    return


@app.cell
def _():
    # gao_pdf_recol.write_csv("Gao-et-al-TCGA-fusion-set.recolumned.tsv", separator="\t")

    # gao_pdf_recol.write_parquet("Gao-et-al-TCGA-fusion-set.recolumned.parquet")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Processing Mitelman's DB of Gene Fusions for NeonDisco

    We are going to take schema files from Mitelman's DB and wrangle it into a simple TSV containing unique list of detected gene fusions from cancers. Note that one of the files, `/home/ec2-user/repos/FT-NeonDisco/data/mitelmans-db/mitelmans-db-fusions-v2025.txt`
    """
    )
    return


@app.cell
def _(pl):
    # load tsvs into polars
    mitelmans_fusions_series = pl.read_csv("/home/ec2-user/repos/FT-NeonDisco/data/mitelmans-db/mitelmans-db-fusions-v2025.txt").to_dict()
    return (mitelmans_fusions_series,)


@app.cell
def _(mitelmans_fusions_series):
    mitelmans_fusions_series
    return


@app.cell
def _(mitelmans_fusions_series):
    # turn to set
    mitelmans_fus_set = set(mitelmans_fusions_series['uniqueFusions'])
    print(mitelmans_fus_set)
    return (mitelmans_fus_set,)


@app.cell
def _(mitelmans_fus_set):
    if 'PARG::BMS1' in mitelmans_fus_set:
        print('YASSSS')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Load Mitelman's DB Morphology Codes TSV & Morphology Annotation TSV""")
    return


@app.cell
def _(pl):
    # load tsv into polars
    mitelmans_morpho_df = pl.read_csv("/home/ec2-user/repos/FT-NeonDisco/data/mitelmans-db/MBCA.TXT.DATA.mitelmans-genesym-morph-db.tsv", separator='\t').cast({"Morph": pl.Utf8})
    mitelmans_morpho_df
    return (mitelmans_morpho_df,)


@app.cell
def _(pl):
    # load tsv
    mitelmans_morpho_annot_df = pl.read_csv("/home/ec2-user/repos/FT-NeonDisco/data/mitelmans-db/KODER.TXT.DATA.mitelmans-morph-codes.tsv", separator='\t').cast({"Kod": pl.Utf8})
    mitelmans_morpho_annot_df
    return (mitelmans_morpho_annot_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Cross-Referencing The Files to Create a Clean Annotated Fusion DataFrame

    Now we can combine the data from the three files into one clean TSV.
    """
    )
    return


@app.cell
def _(mitelmans_fus_set, mitelmans_morpho_annot_df, mitelmans_morpho_df, pl):
    # Convert the fusion set to a dataframe
    fusions_df = pl.DataFrame({"fusionGenePair": list(mitelmans_fus_set)})

    # Since GeneShort is unreliable, we need to find another way to link fusions to morphology codes
    # Let's look at what linking fields we have available in the morphology dataframe

    # First, let's see if we can extract valid gene fusions from GeneShort where they exist
    morpho_clean = mitelmans_morpho_df.filter(
        pl.col("GeneShort").is_not_null() & 
        pl.col("GeneShort").str.contains("::")
    ).select(["GeneShort", "Morph"]).rename({"GeneShort": "fusionGenePair"})

    # Now create the combined dataframe by joining fusions with morphology data
    mitelman_combined_df = fusions_df.join(
        morpho_clean,
        on="fusionGenePair",
        how="left"
    ).join(
        mitelmans_morpho_annot_df.select(["Kod", "Benamning"]),
        left_on="Morph",
        right_on="Kod", 
        how="left"
    ).select([
        "fusionGenePair",
        pl.col("Morph").alias("tumorCode"),
        pl.col("Benamning").alias("tumorAnnotation")
    ]).unique()

    mitelman_combined_df
    return fusions_df, mitelman_combined_df


@app.cell
def _(fusions_df, mitelman_combined_df, pl):
    mitelmans_cleanedup_df = mitelman_combined_df.group_by("fusionGenePair").agg([
            pl.col("tumorAnnotation").drop_nulls().str.concat(" | ").alias("tumorAnnotations")
                ]).with_columns([
                    pl.when(pl.col("tumorAnnotations").str.len_chars() == 0)
                        .then(pl.lit("NA"))
                        .otherwise(pl.col("tumorAnnotations"))
                        .alias("tumorAnnotations")
                ]).select([
                    "fusionGenePair",
                    "tumorAnnotations"
                ])

    # Show some stats
    print(f"Total unique fusions: {len(fusions_df)}")
    print(f"Final condensed dataframe rows: {mitelmans_cleanedup_df.height}")

    mitelmans_cleanedup_df
    return (mitelmans_cleanedup_df,)


@app.cell
def _(mitelmans_cleanedup_df):
    # Save to TSV file
    mitelmans_cleanedup_df.sort("fusionGenePair").write_csv("/home/ec2-user/repos/FT-NeonDisco/output/mitelman-et-al-Cancer-Fusion-set/Mitelman_Fusions_with_tumAnnots.tsv", separator="\t")
    mitelmans_cleanedup_df.sort("fusionGenePair").write_parquet("/home/ec2-user/repos/FT-NeonDisco/output/mitelman-et-al-Cancer-Fusion-set/Mitelman_Fusions_with_tumAnnots.parquet")
    print("Saved to: output/mitelman-et-al-Cancer-Fusion-set/Mitelman_Fusions_with_tumAnnots.tsv")
    return


if __name__ == "__main__":
    app.run()
