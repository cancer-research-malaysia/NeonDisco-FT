### **Fusion Transcript (FT) Raw Data Exploration**

#### **Introduction**

This notebook details the processes (semi-automated) done to further process the raw output files from Arriba and FusionCatcher fusion transcript callers. 

1. Run the `pypolars-process-ft-tsv.py` script to generate fusion transcript list from Arriba and FusionCatcher output files. The script takes a mandatory input of path to the directory where sample-specific fusion call output files from Arriba or FusionCatcher are stored as the first argument, and the specific string that is used to identify tool name (`arr` for Arriba fusion transcript call output file prefix, for instance). 

	For example:
	> ``` pypolars-process-ft-tsv.py data/FTmyBRCAs_raw/Arriba arr ```

	Do the same for the FusionCatcher raw output files, as well as the same Arriba and FusionCatcher output files generated from the processing 113 TCGA-Normals (to use as a panel of normals for FT filtering).

2. Then, load up the two datasets on Jupyter Notebook and concatenate the dataframes together so that Arriba+FusionCatcher unfiltered FT data are combined into one data table and saved in one `.parquet` and `.tsv` file. Do the same for the `TCGANormals` panel of normals.


```python
# first, import packages
import polars as pl
import pandas as pd

from itables import init_notebook_mode, show
#init_notebook_mode(all_interactive=True)
import itables.options as opt
opt.maxBytes = "100KB"

pd.set_option('display.html.table_schema', False)
pd.set_option('display.html.use_mathjax', False)
```


```python
# load up MyBrCa datasets
arr_mdf = pl.scan_parquet('output/MyBrCa/Arriba-FT-all-unfilt-list-v2.parquet')
fc_mdf = pl.scan_parquet('output/MyBrCa/FusionCatcher-FT-all-unfilt-list-v2.parquet')

# now load TCGANormals
arr_norms_mdf = pl.scan_parquet('output/TCGANormals/Arriba-Normal-FT-all-unfilt-list-v2.parquet')
fc_norms_mdf = pl.scan_parquet('output/TCGANormals/FusionCatcher-Normal-FT-all-unfilt-list-v2.parquet')
```

##### **Dataset 1A** (MyBrCa): Arriba unfiltered


```python
from IPython.display import display, HTML, Markdown

display(HTML("Arriba MyBrCa datatable dimension: " + f"<b>{arr_mdf.collect().shape}</b>"))

pl.DataFrame.to_pandas(arr_mdf.collect(), use_pyarrow_extension_array=True).head().to_markdown()

```


Arriba MyBrCa datatable dimension: <b>(49465, 11)</b>





    "|    | fusionTranscriptID                       | fusionGeneID     | breakpointID            | strand1   | strand2   | site1           | site2             | type                  | confidence   |   sampleID | toolID   |\n|---:|:-----------------------------------------|:-----------------|:------------------------|:----------|:----------|:----------------|:------------------|:----------------------|:-------------|-----------:|:---------|\n|  0 | TRMT11::SMG6__6:125986622-17:2244719     | TRMT11::SMG6     | 6:125986622-17:2244719  | +         | -         | CDS/splice-site | CDS/splice-site   | translocation         | high         |          1 | Arriba   |\n|  1 | STAG3::MEF2C-AS1__7:100189570-5:88919251 | STAG3::MEF2C-AS1 | 7:100189570-5:88919251  | +         | -         | CDS             | intron            | translocation/5'-5'   | low          |          1 | Arriba   |\n|  2 | MAPK13::C1QL1__6:36132629-17:44965446    | MAPK13::C1QL1    | 6:36132629-17:44965446  | +         | +         | CDS             | intron            | translocation/5'-5'   | low          |          1 | Arriba   |\n|  3 | STX16::NPEPL1__20:58673711-20:58691724   | STX16::NPEPL1    | 20:58673711-20:58691724 | +         | +         | CDS/splice-site | 5'UTR/splice-site | deletion/read-through | low          |          1 | Arriba   |\n|  4 | MAPK13::NMT1__6:36132629-17:44965446     | MAPK13::NMT1     | 6:36132629-17:44965446  | +         | +         | CDS             | intron            | translocation         | low          |          1 | Arriba   |"



##### **Dataset 1B** (MyBrCa): FusionCatcher unfiltered


```python
from IPython.display import display, HTML

display(HTML("FusionCatcher MyBrCa datatable dimension: " + f"<b>{fc_mdf.collect().shape}</b>"))

pl.DataFrame.to_pandas(fc_mdf.collect(), use_pyarrow_extension_array=True).head()
```


FusionCatcher MyBrCa datatable dimension: <b>(31364, 11)</b>





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe tex2jax_ignore">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>fusionTranscriptID</th>
      <th>fusionGeneID</th>
      <th>breakpointID</th>
      <th>strand1</th>
      <th>strand2</th>
      <th>site1</th>
      <th>site2</th>
      <th>type</th>
      <th>confidence</th>
      <th>sampleID</th>
      <th>toolID</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>SIDT2::TAGLN__11:117195915-11:117203002</td>
      <td>SIDT2::TAGLN</td>
      <td>11:117195915-11:117203002</td>
      <td>+</td>
      <td>+</td>
      <td>CDS(truncated)</td>
      <td>UTR</td>
      <td>.</td>
      <td>.</td>
      <td>2</td>
      <td>FusionCatcher</td>
    </tr>
    <tr>
      <th>1</th>
      <td>AZGP1::GJC3__7:99971746-7:99923603</td>
      <td>AZGP1::GJC3</td>
      <td>7:99971746-7:99923603</td>
      <td>-</td>
      <td>-</td>
      <td>in-frame</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>2</td>
      <td>FusionCatcher</td>
    </tr>
    <tr>
      <th>2</th>
      <td>NPEPPS::TBC1D3__17:47592545-17:38191030</td>
      <td>NPEPPS::TBC1D3</td>
      <td>17:47592545-17:38191030</td>
      <td>+</td>
      <td>-</td>
      <td>CDS(complete)</td>
      <td>UTR</td>
      <td>.</td>
      <td>.</td>
      <td>2</td>
      <td>FusionCatcher</td>
    </tr>
    <tr>
      <th>3</th>
      <td>CYP4F11::CYP4F23P__19:15914762-19:15583580</td>
      <td>CYP4F11::CYP4F23P</td>
      <td>19:15914762-19:15583580</td>
      <td>-</td>
      <td>+</td>
      <td>CDS(truncated)</td>
      <td>exonic(no-known-CDS)</td>
      <td>.</td>
      <td>.</td>
      <td>2</td>
      <td>FusionCatcher</td>
    </tr>
    <tr>
      <th>4</th>
      <td>SLC49A3::ATP5ME__4:686532-4:673401</td>
      <td>SLC49A3::ATP5ME</td>
      <td>4:686532-4:673401</td>
      <td>-</td>
      <td>-</td>
      <td>out-of-frame</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>3</td>
      <td>FusionCatcher</td>
    </tr>
  </tbody>
</table>
</div>



##### **Dataset 2A** (TCGA Normals): Arriba unfiltered


```python
from IPython.display import display, HTML

display(HTML("Arriba TCGA-Normals datatable dimension: " + f"<b>{arr_norms_mdf.collect().shape}</b>"))

pl.DataFrame.to_pandas(arr_norms_mdf.collect(), use_pyarrow_extension_array=True).head()
```

##### **Dataset 2B** (TCGA Normals): FusionCatcher unfiltered


```python
from IPython.display import display, HTML

display(HTML("FusionCatcher TCGA-Normals datatable dimension: " + f"<b>{fc_norms_mdf.collect().shape}</b>"))

pl.DataFrame.to_pandas(fc_norms_mdf.collect(), use_pyarrow_extension_array=True).head()
```

#### **Concatenate Arriba and FusionCatcher Datasets**

Now, we can merge the two dataframes into one masterFrame for each cohort data (MyBrCa & TCGA panel of normals) using Polars' `concat` (vertical concatenation is the default, where two dataframes sharing the exact same columns would be joined together, adding all rows of dataframe 1 and 2 vertically).


```python
%%capture --no-stdout --no-display

joined_df = pl.concat(
    [
        arr_mdf.collect(),
        fc_mdf.collect()
    ]
)
from IPython.display import display, HTML

display(HTML("Concatenated MyBrCa Arriba+FusionCatcher datatable dimension: " + f"<b>{joined_df.shape}</b>"))

pl.DataFrame.to_pandas(joined_df, use_pyarrow_extension_array=True).head()
```

Do the same with the TCGA panel of normal FTs.


```python
%%capture --no-stdout --no-display

joined_norms_df = pl.concat(
    [
        arr_norms_mdf.collect(),
        fc_norms_mdf.collect()
    ]
)
from IPython.display import display, HTML

display(HTML("Concatenated TCGA Normals Arriba+FusionCatcher datatable dimension: " + f"<b>{joined_norms_df.shape}</b>"))

pl.DataFrame.to_pandas(joined_norms_df, use_pyarrow_extension_array=True).head()
```

#### **Filter MyBrCa Merged Datatable using Panel of Normals**
Now we can filter the unfiltered FT datatable by discarding those that are present in TCGA Normal datatable.



```python
# first load the parquet files
# load up MyBrCa datasets
mybrca_ccdf = pl.scan_parquet('output/MyBrCa/Arr_FC-concat-FT-all-unfilt-list-v2.parquet')

# now load TCGANormals
norms_ccdf = pl.scan_parquet('output/TCGANormals/Arr_FC-Normals-concat-FT-all-unfilt-list-v2.parquet')
```


```python
mybrca_ccdf_pan = pl.DataFrame.to_pandas(mybrca_ccdf.collect(), use_pyarrow_extension_array=True)

mybrca_ccdf_pan
```


```python
norms_ccdf_pan = pl.DataFrame.to_pandas(norms_ccdf.collect(), use_pyarrow_extension_array=True)

norms_ccdf_pan
```

Use Polars' `filter` expression with `is_in` and the negation `~` to keep only unique rows for column `breakpointID` in MyBrCa dataframe that are NOT in the `breakpointID` in TCGANormals dataframe. 


```python
normfilt_mybrca_ccdf = mybrca_ccdf.collect().filter(~pl.col('breakpointID').is_in(norms_ccdf.collect()['breakpointID']))

print(normfilt_mybrca_ccdf.shape)
normfilt_mybrca_ccdf.head()
```

Use `group_by` and `n_unique()` in Polars to create a count table for unique `breakpointID` and how "shared" it is across our MyBrCa cohort. 

Here, we subset the dataframe into just `breakpointID` and `sampleID` and then use `group_by` on the `breakpointID` column, then counting number of unique occurences of each unique `breakpointID` in the `sampleID` column. 

This would return a count of unique samples (patients) one particular unique breakpoint appears in. I call this the `sharednessDegree`.

> **NOTE:** This is the best way to address miscounting breakpoints that appear in multiple rows due to differences in gene naming but they are only seen in one sample. Using other counting strategies such as window function (`.over` method) will count these duplicate rows as separate entities when in reality they are the same breakpoint seen in just one patient.


```python
normfilt_mybrca_sharedness = normfilt_mybrca_ccdf.select(pl.col(["breakpointID", "sampleID"])).group_by("breakpointID").n_unique().rename({"sampleID": "sharednessDegree"})

normfilt_mybrca_sharedness.sort("sharednessDegree", descending=True)
```

#### **Plotting the Sharedness Degree**

We have used Polars to easily group and count the number of patients sharing a particular breakpoint ID for each unique breakpoint ID as above, let's formalize that again by using Pandas instead.

First, subset the filtered dataframe to just the two columns we are interested in using Polars, but this time prepend the "P" string to all values of the `sampleID` column, then convert to Pandas for visualization.


```python
bp_sample_array = normfilt_mybrca_ccdf.select(
    pl.col("breakpointID"),
    pl.concat_str(pl.lit("P"), pl.col("sampleID")).alias("sampleID")
)

bpsample_pdf = bp_sample_array.to_pandas()
bpsample_pdf
```

Due to the annotation redundancy in `fusionGeneID` column in the original df, we now have rows in `breakpointID` and `sampleID` that are repeated (i.e. `6:36132629-17:44965446	P1` as seen above). Let's filter these out, as they represent the same putative FT.


```python
# Drop duplicates based on both columns
bpsample_pdf_unique = bpsample_pdf.drop_duplicates()

# see how many duplicates were removed
print("Original number of rows:", len(bpsample_pdf))
print("Number of rows after removing duplicates:", len(bpsample_pdf_unique))
print("Number of duplicates removed:", len(bpsample_pdf) - len(bpsample_pdf_unique))
```


```python
bpsample_pdf_unique
```

Now group by each unique `breakpointID` and count how many `sampleID` is associated with this breakpoint (*sharedness degree*).


```python
# Group by breakpointID and count unique sampleIDs
breakpoint_counts = bpsample_pdf_unique.groupby('breakpointID')['sampleID'].nunique().reset_index()

# Rename the column for clarity
breakpoint_counts = breakpoint_counts.rename(columns={'sampleID': 'sharednessDegree'})

show(breakpoint_counts, maxBytes=0)

```

Then we can count the number of unique `breakpointID`s for each `sharednessDegree` value.


```python
# Count the number of unique breakpointIDs for each sharednessDegree value
sharedness_counts = (
    breakpoint_counts
    .groupby('sharednessDegree')
    .agg(
        unique_bp_count=('breakpointID', 'nunique')
    )
    .reset_index()
    .sort_values('sharednessDegree')
)

sharedness_counts.reset_index(drop=True)
```


```python
import matplotlib.pyplot as plt
import seaborn as sns

# Create the bar plot
plt.figure(figsize=(10, 8), dpi=300)
sns.barplot(x=sharedness_counts['sharednessDegree'],y=sharedness_counts['unique_bp_count'], color='steelblue')

# Add value labels on top of the bars
for i, v in enumerate(sharedness_counts['unique_bp_count']):
    plt.text(i, v, str(v), color='black', ha='center', fontweight='bold', fontsize=8)
    
# Set labels and title
plt.xlabel('Sharedness Degree (Number of Patients A Unique FT Is Observed)')
plt.ylabel('Count of Unique FTs')
plt.title('Frequency of Unique Tumor-Specific FTs by Sharedness Degree')

# Rotate x-axis labels for better readability
plt.xticks(rotation=90)

plt.show()
```

### Using Graph Theory to Investigate Bipartite Relationship

We can use graph theory to explore the underlying bipartite network between unique `breakpointID` and `sampleID`.

#### Design an Analysis Class
Create a complex class called `NetworkAnalyzer` to do graph network analysis between `breakpointID` and `sampleID`. 


```python
import numpy as np
import networkx as nx
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.spatial.distance import pdist, squareform
from networkx.algorithms.bipartite import density as bipartite_density
from scipy.sparse import csr_matrix

class NetworkAnalyzer:
	def __init__(self, df=None, patient_col=None, breakpoint_col=None, 
					precomputed_matrix=None, patients=None, breakpoints=None):
		"""
		Initialize NetworkAnalyzer with either DataFrame or precomputed matrix.
		
		Args:
			df (pl.DataFrame.Polars, optional): Input DataFrame HAS TO BE IN POLARS
			patient_col (str, optional): Column name for patients
			breakpoint_col (str, optional): Column name for breakpoints
			precomputed_matrix (csr_matrix, optional): Pre-computed sparse adjacency matrix
			patients (list, optional): List of patient IDs (required if using precomputed_matrix)
			breakpoints (list, optional): List of breakpoint IDs (required if using precomputed_matrix)
		"""
		if precomputed_matrix is not None:
			if patients is None or breakpoints is None:
				raise ValueError("Must provide patients and breakpoints lists with precomputed matrix")
			# Keep matrix in sparse format
			self.adj_matrix_sparse = precomputed_matrix
			self.patients = patients
			self.breakpoints = breakpoints
		elif df is not None and patient_col and breakpoint_col:
			self.df = df
			self.patient_col = patient_col
			self.breakpoint_col = breakpoint_col
			
			# Get unique sets
			self.patients = sorted(df[patient_col].unique().to_list())
			self.breakpoints = sorted(df[breakpoint_col].unique().to_list())
			
			# Create sparse adjacency matrix
			self.adj_matrix_sparse, self.patient_idx_dict, self.breakpoint_idx_dict = self._create_adjacency_matrix()
		else:
			raise ValueError("Must provide either DataFrame with column names or precomputed matrix with labels")

		# Don't calculate metrics immediately - do it lazily
		self._metrics_calculated = False
		
	def _ensure_metrics_calculated(self):
		"""Calculate metrics if they haven't been calculated yet."""
		if not self._metrics_calculated:
			self._calculate_metrics()
			self._metrics_calculated = True

	def _create_adjacency_matrix(self) -> csr_matrix:
		"""Create the sparse adjacency matrix from the input DataFrame."""
		matrix = np.zeros((len(self.patients), len(self.breakpoints)))
		connections = self.df.group_by(self.patient_col).agg(
			pl.col(self.breakpoint_col).alias('breakpoints')
		).to_dict(as_series=False)

		patient_idx = {p: i for i, p in enumerate(self.patients)}
		breakpoint_idx = {b: i for i, b in enumerate(self.breakpoints)}

		for i, patient in enumerate(connections[self.patient_col]):
			for bp in connections['breakpoints'][i]:
				matrix[patient_idx[patient]][breakpoint_idx[bp]] = 1

		return csr_matrix(matrix), patient_idx, breakpoint_idx

	def _calculate_metrics(self):
		"""Calculate various network metrics."""
		# Convert to dense only when needed for specific calculations
		dense_matrix = self.adj_matrix_sparse.toarray()
		
		self.patient_degrees = np.asarray(self.adj_matrix_sparse.sum(axis=1)).flatten()
		self.breakpoint_degrees = np.asarray(self.adj_matrix_sparse.sum(axis=0)).flatten()
		
		# Only calculate similarity matrices if needed for visualization
		self.patient_similarity = squareform(pdist(dense_matrix, metric='jaccard'))
		self.breakpoint_similarity = squareform(pdist(dense_matrix.T, metric='jaccard'))

		# Create bipartite graph more efficiently
		G = nx.Graph()
		G.add_nodes_from(range(len(self.patients)), bipartite=0)
		G.add_nodes_from(range(len(self.patients), len(self.patients) + len(self.breakpoints)), bipartite=1)
		
		# Add edges using sparse matrix coordinates
		rows, cols = self.adj_matrix_sparse.nonzero()
		edges = zip(rows, cols + len(self.patients))
		G.add_edges_from(edges)

		self.density = bipartite_density(G, range(len(self.patients), len(self.patients) + len(self.breakpoints)))
		
		# Calculate centrality
		centrality = nx.degree_centrality(G)
		self.breakpoint_centrality = [centrality[i + len(self.patients)] for i in range(len(self.breakpoints))]

	def save_matrix(self, filename):
		"""
		Save the adjacency matrix in CSR format along with patient and breakpoint labels.
		
		Args:
			filename (str): Base filename to save the data (without extension)
		"""
		# Save the sparse matrix
		sparse_matrix = self.adj_matrix_sparse
		np.savez(f"{filename}_adjac_matrix.npz",
					data=sparse_matrix.data,
					indices=sparse_matrix.indices,
					indptr=sparse_matrix.indptr,
					shape=sparse_matrix.shape)
		
		# Save the labels
		np.save(f"{filename}_matrix_label_patients.npy", np.array(self.patients))
		np.save(f"{filename}_matrix_label_breakpoints.npy", np.array(self.breakpoints))

	@classmethod
	def load_from_files(cls, filename):
		"""
		Load a NetworkAnalyzer instance from saved files.
		
		Args:
			filename (str): Base filename (without extension) used when saving
			
		Returns:
			NetworkAnalyzer: New instance with loaded data
		"""
		# Load the sparse matrix
		loader = np.load(f"{filename}_adjac_matrix.npz")
		matrix = csr_matrix((loader['data'], loader['indices'], loader['indptr']), shape=loader['shape'])
		
		# Load the labels
		patients = np.load(f"{filename}_matrix_label_patients.npy").tolist()
		breakpoints = np.load(f"{filename}_matrix_label_breakpoints.npy").tolist()

		return cls(precomputed_matrix=matrix, patients=patients, breakpoints=breakpoints)
	
	def create_adjacency_matrix_plot(self, top_bins: list = None, bottom_bins: list = None) -> go.Figure:
		"""Create a standalone plot of the patient-breakpoint adjacency matrix."""
		# Ensure metrics are calculated before accessing them
		self._ensure_metrics_calculated()

		# Convert to dense only when needed for specific calculations
		dense_matrix = self.adj_matrix_sparse.toarray()

		# If no bins are provided, plot the full adjacency matrix
		if not top_bins and not bottom_bins:
			top_bins = list(range(len(self.breakpoints)))
			top_breakpoints = self.breakpoints
			top_matrix = dense_matrix
		elif top_bins:
			top_breakpoints = [self.breakpoints[i] for i in top_bins]
			top_matrix = dense_matrix[:, top_bins]

		# Create the top breakpoints plot
		fig = go.Figure(
			data=go.Heatmap(
				z=top_matrix,
				x=top_breakpoints,
				y=self.patients,
				colorscale="Blues",
				showscale=True,
				hoverongaps=False,
				hoverinfo='text',
				text=[[f"Patient: {p}<br>Breakpoint: {b}<br>Connected: {'Yes' if top_matrix[i][j] else 'No'}"
						for j, b in enumerate(top_breakpoints)]
						for i, p in enumerate(self.patients)],
				colorbar=dict(title="Connection"),
				name="Top Breakpoints"
			)
		)

		# If bottom bins are provided, add the bottom breakpoints plot
		if bottom_bins:
			bottom_breakpoints = [self.breakpoints[i] for i in bottom_bins]
			bottom_matrix = dense_matrix[:, bottom_bins]

			fig.add_trace(
				go.Heatmap(
					z=bottom_matrix,
					x=bottom_breakpoints,
					y=self.patients,
					colorscale="Blues",
					showscale=True,
					hoverongaps=False,
					hoverinfo='text',
					text=[[f"Patient: {p}<br>Breakpoint: {b}<br>Connected: {'Yes' if bottom_matrix[i][j] else 'No'}"
							for j, b in enumerate(bottom_breakpoints)]
							for i, p in enumerate(self.patients)],
					colorbar=dict(title="Connection"),
					name="Bottom Breakpoints"
				)
			)

		fig.update_layout(
			height=1000,
			width=1400,
			title=dict(
				text="Adjacency Matrix",
				x=0.5,
				y=0.95,
				font=dict(size=18)
			),
			xaxis_title="Breakpoints",
			yaxis_title="Patients",
			template="simple_white"
		)

		return fig

	def create_degree_distribution_plot(self, top_bins: list = None, bottom_bins: list = None) -> go.Figure:
		"""Create a standalone plot of the breakpoint degree distribution."""
		# Ensure metrics are calculated before accessing them
		self._ensure_metrics_calculated()

		# If no bins are provided, plot the full degree distribution
		if not top_bins and not bottom_bins:
			top_bins = list(range(len(self.breakpoints)))
			top_breakpoints = self.breakpoints
			top_degrees = [self.breakpoint_degrees[i] for i in top_bins]

		elif top_bins:
			top_breakpoints = [self.breakpoints[i] for i in top_bins]
			top_degrees = [self.breakpoint_degrees[i] for i in top_bins]

		# Create the top breakpoints plot
		fig = go.Figure(
			data=go.Bar(
				x=top_breakpoints,
				y=top_degrees,
				hovertext=[f"Breakpoint: {bp}<br>Connected to {deg} patients"
							for bp, deg in zip(top_breakpoints, top_degrees)],
				hoverinfo='text',
				marker_color='rgb(158,202,225)',
				marker_line_color='rgb(8,48,107)',
				marker_line_width=1.5,
				name="Degree Distribution"
			)
		)

		# If bottom bins are provided, add the bottom breakpoints plot
		if bottom_bins:
			bottom_breakpoints = [self.breakpoints[i] for i in bottom_bins]
			bottom_degrees = [self.breakpoint_degrees[i] for i in bottom_bins]

			fig.add_trace(
				go.Bar(
					x=bottom_breakpoints,
					y=bottom_degrees,
					hovertext=[f"Breakpoint: {bp}<br>Connected to {deg} patients"
								for bp, deg in zip(bottom_breakpoints, bottom_degrees)],
					hoverinfo='text',
					marker_color='rgb(158,202,225)',
					marker_line_color='rgb(8,48,107)',
					marker_line_width=1.5,
					name="Bottom Breakpoints"
				)
			)

		fig.update_layout(
			height=1000,
			width=1400,
			title=dict(
				text="Breakpoint Degree Distribution",
				x=0.5,
				y=0.95,
				font=dict(size=18)
			),
			xaxis_title="Breakpoints",
			yaxis_title="Number of Patients",
			template="simple_white"
		)

		return fig
	
	def create_patient_similarity_matrix_plot(self) -> go.Figure:
		"""Create a standalone plot of the patient similarity matrix."""
		# Ensure metrics are calculated before accessing them
		self._ensure_metrics_calculated()

		# Extract the lower triangular portion of the matrix (excluding the diagonal)
		patient_similarity_lower = np.tril(1 - self.patient_similarity, -1)

		# Create a mask for the upper triangular portion (excluding the diagonal)
		patient_similarity_mask = np.tri(len(self.patients), len(self.patients), k=1, dtype=bool)

		# Create the heatmap data, setting the upper triangular portion to the maximum value
		patient_similarity_data = np.where(patient_similarity_mask, np.max(patient_similarity_lower), 1 - self.patient_similarity)

		fig = go.Figure(
			data=go.Heatmap(
				z=patient_similarity_data,
				x=self.patients,
				y=self.patients,
				colorscale="Viridis",
				showscale=True,
				colorbar=dict(title="Similarity"),
				name="Patient Similarity"
			)
		)

		fig.update_layout(
			height=1000,
			width=1400,
			title=dict(
				text="Patient Similarity Matrix",
				x=0.5,
				y=0.95,
				font=dict(size=18)
			),
			xaxis_title="Patients",
			yaxis_title="Patients",
			template="simple_white"
		)

		return fig

	def create_breakpoint_cooccurrence_plot(self) -> go.Figure:
		"""Create a standalone plot of the breakpoint co-occurrence matrix."""
		# Ensure metrics are calculated before accessing them
		self._ensure_metrics_calculated()

		# Extract the lower triangular portion of the matrix (excluding the diagonal)
		breakpoint_similarity_lower = np.tril(1 - self.breakpoint_similarity, -1)

		# Create a mask for the upper triangular portion (excluding the diagonal)
		breakpoint_similarity_mask = np.tri(len(self.breakpoints), len(self.breakpoints), k=1, dtype=bool)

		# Create the heatmap data, setting the upper triangular portion to the maximum value
		breakpoint_similarity_data = np.where(breakpoint_similarity_mask, np.max(breakpoint_similarity_lower), 1 - self.breakpoint_similarity)

		fig = go.Figure(
			data=go.Heatmap(
				z=breakpoint_similarity_data,
				x=self.breakpoints,
				y=self.breakpoints,
				colorscale="Viridis",
				showscale=True,
				colorbar=dict(title="Co-occurrence"),
				name="Breakpoint Co-occurrence"
			)
		)

		fig.update_layout(
			height=1000,
			width=1400,
			title=dict(
				text="Breakpoint Co-occurrence Matrix",
				x=0.5,
				y=0.95,
				font=dict(size=18)
			),
			xaxis_title="Breakpoints",
			yaxis_title="Breakpoints",
			template="simple_white"
		)

		return fig

	def create_dashboard(self) -> go.Figure:
		"""Create a comprehensive visualization dashboard."""
		# Ensure metrics are calculated before accessing them
		self._ensure_metrics_calculated()
		# Convert to dense only when needed for specific calculations
		dense_matrix = self.adj_matrix_sparse.toarray()
		fig = make_subplots(
			rows=2, cols=2,
			subplot_titles=("Patient-Breakpoint Adjacency Matrix", 
							"Breakpoint Degree Distribution",
							"Patient Similarity Matrix", 
							"Breakpoint Co-occurrence Matrix"),
			specs=[[{"type": "heatmap"}, {"type": "bar"}],
					[{"type": "heatmap"}, {"type": "heatmap"}]]
		)

		# 1. Adjacency Matrix with custom hover text
		hover_text = [[f"Patient: {p}<br>Breakpoint: {b}<br>Connected: {'Yes' if dense_matrix[i][j] else 'No'}"
						for j, b in enumerate(self.breakpoints)]
						for i, p in enumerate(self.patients)]

		fig.add_trace(
			go.Heatmap(
				z=dense_matrix,
				x=self.breakpoints,
				y=self.patients,
				colorscale="Blues",
				showscale=False,
				hoverongaps=False,
				hoverinfo='text',
				text=hover_text,
				colorbar=dict(title="Connection"),
				name="Connections"
			),
			row=1, col=1
		)

		# 2. Degree Distribution with custom hover
		hover_text = [f"Breakpoint: {bp}<br>Connected to {deg} patients"
						for bp, deg in zip(self.breakpoints, self.breakpoint_degrees)]

		fig.add_trace(
			go.Bar(
				x=self.breakpoints,
				y=self.breakpoint_degrees,
				hovertext=hover_text,
				hoverinfo='text',
				marker_color='rgb(158,202,225)',
				marker_line_color='rgb(8,48,107)',
				marker_line_width=1.5,
				name="Breakpoint Degrees"
			),
			row=1, col=2
		)

		# 3. Patient Similarity Matrix
		fig.add_trace(
			go.Heatmap(
				z=1 - self.patient_similarity,
				x=self.patients,
				y=self.patients,
				colorscale="Viridis",
				showscale=False,
				colorbar=dict(title="Similarity"),
				name="Patient Similarity"
			),
			row=2, col=1
		)

		# 4. Breakpoint Co-occurrence
		fig.add_trace(
			go.Heatmap(
				z=1 - self.breakpoint_similarity,
				x=self.breakpoints,
				y=self.breakpoints,
				colorscale="Viridis",
				showscale=True,
				colorbar=dict(title="Co-occurrence"),
				name="Breakpoint Co-occurrence"
			),
			row=2, col=2
		)

		fig.update_layout(
			height=1000,
			width=1200,
			title=dict(
				text="Network Analysis Dashboard",
				x=0.5,
				y=0.95,
				font=dict(size=24)
			),
			showlegend=False,
			template="simple_white"
		)

		font_size = 14
		fig.update_xaxes(title_text="Breakpoints", title_font=dict(size=font_size), row=1, col=1)
		fig.update_yaxes(title_text="Patients", title_font=dict(size=font_size), row=1, col=1)
		fig.update_xaxes(title_text="Breakpoints", title_font=dict(size=font_size), row=1, col=2)
		fig.update_yaxes(title_text="Number of Patients", title_font=dict(size=font_size), row=1, col=2)
		fig.update_xaxes(title_text="Patients", title_font=dict(size=font_size), row=2, col=1)
		fig.update_yaxes(title_text="Patients", title_font=dict(size=font_size), row=2, col=1)
		fig.update_xaxes(title_text="Breakpoints", title_font=dict(size=font_size), row=2, col=2)
		fig.update_yaxes(title_text="Breakpoints", title_font=dict(size=font_size), row=2, col=2)

		return fig
	
	def get_breakpoint_bins(self, top_percentile: float = 0.001, bottom_percentile: float = 0.001) -> tuple:
		"""
		Calculate the indexes of the breakpoints at the specified percentiles.
		Returns a tuple of two lists:
		- The first list contains the indexes of the top `top_percentile` breakpoints by degree.
		- The second list contains the indexes of the bottom `bottom_percentile` breakpoints by degree.
		"""
		# Ensure metrics are calculated before accessing them
		self._ensure_metrics_calculated()
		sorted_degrees = sorted(self.breakpoint_degrees)
		top_cutoff = int(len(sorted_degrees) * top_percentile)
		bottom_cutoff = int(len(sorted_degrees) * (1 - bottom_percentile))

		top_bins = [i for i in range(top_cutoff)]
		bottom_bins = [i for i in range(bottom_cutoff, len(sorted_degrees))]

		return top_bins, bottom_bins

	def print_summary_stats(self):
		# Ensure metrics are calculated before accessing them
		self._ensure_metrics_calculated()
		print(f"Network Summary Statistics:")
		print(f"---------------------------")
		print(f"Number of Patients: {len(self.patients)}")
		print(f"Number of Breakpoints: {len(self.breakpoints)}")
		print(f"Network Density: {self.density:.3f}")
		print(f"Average Patient Degree: {np.mean(self.patient_degrees):.2f}")
		print(f"Average Breakpoint Degree: {np.mean(self.breakpoint_degrees):.2f}")
		print(f"\nTop Breakpoints by Degree:")
		for bp, degree in sorted(zip(self.breakpoints, self.breakpoint_degrees), 
								key=lambda x: x[1], reverse=True)[:5]:
			print(f"  {bp}: {degree}")
		print(f"\nTop 10 Breakpoints by Degree Centrality:")
    	# Create list of (breakpoint, centrality) tuples and sort by centrality
		centrality_pairs = list(zip(self.breakpoints, self.breakpoint_centrality))
		sorted_by_centrality = sorted(centrality_pairs, key=lambda x: x[1], reverse=True)
		
		# Print top 10
		for bp, centrality in sorted_by_centrality[:10]:
			print(f"  {bp}: {centrality:.3f}")
```

#### Test Class on Toy Data


```python
# create toy data
np.random.seed(420)  # for reproducibility

patients = [f'P{i}' for i in range(1, 21)]  # 20 patients
breakpoints = [f'BP{i}' for i in range(1, 16)]  # 15 breakpoints

# Create random connections (each patient has 2-6 breakpoints)
data = []
for patient in patients:
	num_breakpoints = np.random.randint(2, 7)
	patient_breakpoints = np.random.choice(breakpoints, size=num_breakpoints, replace=False)
	for bp in patient_breakpoints:
		data.append({'patient_id': patient, 'breakpoint_id': bp})

# Create Polars DataFrame
df = pl.DataFrame(data)

# now test the class
# Create and display visualization
analyzer = NetworkAnalyzer(df, patient_col='patient_id', breakpoint_col='breakpoint_id')
fig = analyzer.create_dashboard()
# fig.show()

# Print summary statistics
analyzer.print_summary_stats()
```

##### Instantiate NetworkAnalyzer class on MyBrCa Data

Now instantiate the class we built on our normal-filtered, unique-breakpoint-only subset dataFrame from MyBrCa FT data.


```python
# This is DF with redundant duplicate breakpointID-sampleID pairing
print(bp_sample_array.shape)
bp_sample_array.head()
```


```python
# This is DF with redundant duplicate breakpointID-sampleID pairing FILTERED OUT
bpsample_poldf = pl.from_pandas(bpsample_pdf_unique)
print(bpsample_poldf.shape)
bpsample_poldf.head()
```

##### Saving adjacency matrix to file as precomputed npz file


```python
# # call the method save_matrix directly
# analyzer_my.save_matrix('output/adj_matrix_ALLDEG_mybrca_v2')

# # try reloading into an instance using the decorated class method load_from_file
# analyzer_my_reloaded = NetworkAnalyzer.load_from_files('output/adj_matrix_ALLDEG_mybrca_v2')

# # Print summary statistics
# analyzer_my.print_summary_stats()
```

#### Keeping FT Breakpoints Seen in More than 10 Patients (1% of the MyBrCa Cohort)

The distribution of the Sharedness Degree of each unique breakpoint, is as expected, skewed towards having a lot of unique, patient-specific connections, and very few shared breakpoints across patients. 

We can try to visualize the adjacency matrix, but because of the massive matrix dimension we have (**988 patients x 43927 unique breakpoints**), it is best to first filter out patient-specific breakpoints first. In fact, due to the fact that the putative FT neoantigen distribution is so skewed towards individualized presence, let's create a filtering threshold of keeping only the breakpoint IDs that are seen in **more than 9 patients (approximately 1% of the MyBrCa cohort)**. 



```python
## go back to the sharednessDegree Pandas dataFrame
# # select rows that has sharednessDegree > 10

bp_sharedness_gt9 = pl.from_pandas(breakpoint_counts).filter(pl.col('sharednessDegree') > 9)

show(bp_sharedness_gt9, maxBytes=0)

#### or directly use the Polars dataframe normfilt_mybrca_sharedness
# 
# bp_sharedness_gt9 = normfilt_mybrca_sharedness.filter(pl.col('sharednessDegree') > 9)
```

Now we can use the unique, filtered, thresholded elements in the `breakpointID` column of the filtered dataFrame above, as the filtering list to keep only these same breakpoints in the other dataFrame used to instantiate `NetworkAnalyzer`. 


```python
# bp_sharedness_gt9 is the dataframe with unique breakpointIDs to be used as filter
# bpsample_poldf is the dataframe to be filtered

filt_bpsample_poldf = bpsample_poldf.filter(
    pl.col("breakpointID").is_in(bp_sharedness_gt9["breakpointID"])
)

show(filt_bpsample_poldf, maxBytes=0)
```


```python
analyzer_my_filt = NetworkAnalyzer(filt_bpsample_poldf, patient_col='sampleID', breakpoint_col='breakpointID')
```


```python
# Print summary statistics
analyzer_my_filt.print_summary_stats()
```


```python
# call the method save_matrix directly; UNCOMMENT TO SAVE
# analyzer_my_filt.save_matrix('output/adj_matrix_GT9FILT_mybrca_v2')
```

##### Plot Adjacency Matrix


```python
plot = analyzer_my_filt.create_adjacency_matrix_plot()
plot.show()
```


```python
plot = analyzer_my_filt.create_degree_distribution_plot()
plot.show()
```


```python
plot = analyzer_my_filt.create_patient_similarity_matrix_plot()
plot.show()
```


```python
plot = analyzer_my_filt.create_breakpoint_cooccurrence_plot()
plot.show()
```


```python
filt_matrix = analyzer_my_filt.adj_matrix_sparse.toarray()
```


```python
def find_minimal_covering_subset(adjacency_matrix, coverage_threshold, label_to_index_dict=None):
    """
    Find minimal subset of rows that covers at least coverage_threshold fraction of columns,
    with handling for dictionary-based label mapping.

    Parameters:
    -----------
    adjacency_matrix : np.ndarray
        Binary matrix where rows are members of set A and columns are members of set B
        1 indicates overlap, 0 indicates no overlap
    coverage_threshold : float
        Fraction of set B that needs to be covered (between 0 and 1)
    label_to_index_dict : dict, optional
        Dictionary mapping labels (strings) to indices (int)
        Example: {'label1': 0, 'label2': 1, ...}

    Returns:
    --------
    tuple
        (selected_indices, selected_labels, actual_coverage)
        - selected_indices: List of numerical indices of selected rows
        - selected_labels: List of original labels corresponding to the indices
        - actual_coverage: Achieved coverage fraction
    """
    # Create reverse mapping from index to label
    if label_to_index_dict is not None:
        index_to_label = {v: k for k, v in label_to_index_dict.items()}
    
    # Work with transpose of the matrix
    working_matrix = adjacency_matrix.T
    num_rows, num_cols = working_matrix.shape
    
    # Calculate target coverage
    target_coverage = int(np.ceil(num_cols * coverage_threshold))
    print(f"Target coverage: {target_coverage} columns out of {num_cols}")
    
    # Initialize tracking variables
    selected_rows = []
    covered_cols = np.zeros(num_cols, dtype=bool)
    
    while np.sum(covered_cols) < target_coverage:
        # Calculate coverage gains for remaining rows
        available_rows = [i for i in range(num_rows) if i not in selected_rows]
        
        if not available_rows:
            break
            
        coverage_gains = np.array([
            np.sum(~covered_cols & (working_matrix[i] == 1))
            for i in available_rows
        ])
        
        if np.max(coverage_gains) == 0:
            print("No more improvements possible")
            break
            
        # Select the row that covers the most new columns
        best_row_idx = available_rows[np.argmax(coverage_gains)]
        selected_rows.append(best_row_idx)
        
        # Print progress with label if available
        if label_to_index_dict is not None:
            label = index_to_label[best_row_idx]
            new_coverage = np.sum(~covered_cols & (working_matrix[best_row_idx] == 1))
            print(f"Selected {label} (index {best_row_idx}) covering {new_coverage} new columns")
        
        # Update covered columns
        covered_cols = covered_cols | (working_matrix[best_row_idx] == 1)
    
    # Calculate actual coverage achieved
    actual_coverage = np.sum(covered_cols) / num_cols
    print(f"Achieved {actual_coverage:.2%} coverage")
    
    # Convert indices to labels if dictionary provided
    if label_to_index_dict is not None:
        selected_labels = [index_to_label[idx] for idx in selected_rows]
    else:
        selected_labels = None
    
    return selected_rows, selected_labels, actual_coverage

```


```python
breakpoint_dict = analyzer_my_filt.breakpoint_idx_dict
filt_bp_indices, filt_bp_labels, coverage = find_minimal_covering_subset(filt_matrix, coverage_threshold=1.0, label_to_index_dict=breakpoint_dict)
print(f"Minimal coverage: {coverage}")
print(f"Minimal set of breakpoints: {filt_bp_indices}; Length of set: {len(filt_bp_indices)}")
print(f"Labels of the minimal cover set of breakpoints: {filt_bp_labels}")
```


```python
patient_dict = analyzer_my_filt.patient_idx_dict
print(patient_dict)
```


```python
# subset original matrix
minimal_set_cover_subset_matrix = filt_matrix[:, filt_bp_indices]
minimal_set_cover_subset_matrix.shape
```


```python
# subset the original df used to generate analyzer_my_filt instance

filt_bpsample_minsetcover_poldf = bpsample_poldf.filter(
    pl.col("breakpointID").is_in(filt_bp_labels)
)

show(filt_bpsample_minsetcover_poldf, maxBytes=0)

```


```python
# create a new NetworkAnalyzer instance
analyzer_subset_filt = NetworkAnalyzer(filt_bpsample_minsetcover_poldf, patient_col="sampleID", breakpoint_col="breakpointID")
```


```python
# Print summary statistics
analyzer_subset_filt.print_summary_stats()
```


```python
plot = analyzer_subset_filt.create_adjacency_matrix_plot()
plot.show()
```

#### **Filtering Out non-TNBCs**

We can filter out rows corresponding to `sampleID` more than 172, because these are not TNBC samples.




