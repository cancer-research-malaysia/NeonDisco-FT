### **Fusion Transcript Data Wrangling**

#### Concatenating and Filtering Raw Fusion Transcript Output from Arriba and FusionCatcher

This notebook details the processes (semi-automated) done to further process the raw output files from Arriba and FusionCatcher fusion transcript callers. 

1. Run the `wrangle-ft-tsv.py` script to generate fusion transcript list from Arriba and FusionCatcher output files. The script takes a mandatory input of path to the directory where sample-specific fusion call output files from Arriba or FusionCatcher are stored as the first argument, and the specific string that is used to identify tool name (`arr` for Arriba fusion transcript call output file prefix, for instance). 

For example:
> ``` wrangle-ft-tsv.py data/FTmyBRCAs_raw/Arriba arr ```


```python

```
