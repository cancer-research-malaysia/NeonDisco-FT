import marimo

__generated_with = "0.12.10"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    return mo, pl


@app.cell
def _(pl):
    # load tsv into polars
    gao_df = pl.read_csv("/Users/sufyazi/Library/CloudStorage/OneDrive-CancerResearchMalaysia/NeonDisco/onedrive-cloud-tmp/Gao-et-al-TCGA-fusion-set.tsv.txt", separator='\t')
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
    return (gao_pdf_recol,)


@app.cell
def _(gao_pdf_recol):
    gao_pdf_recol.write_csv("Gao-et-al-TCGA-fusion-set.recolumned.tsv", separator="\t")

    gao_pdf_recol.write_parquet("Gao-et-al-TCGA-fusion-set.recolumned.parquet")
    return


if __name__ == "__main__":
    app.run()
