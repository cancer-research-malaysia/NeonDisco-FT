### **Fusion Transcript Data Wrangling**

#### Concatenating and Filtering Raw Fusion Transcript Output from Arriba and FusionCatcher

This notebook details the processes (semi-automated) done to further process the raw output files from Arriba and FusionCatcher fusion transcript callers. 

1. Run the `wrangle-ft-tsv.py` script to generate fusion transcript list from Arriba and FusionCatcher output files. The script takes a mandatory input of path to the directory where sample-specific fusion call output files from Arriba or FusionCatcher are stored as the first argument, and the specific string that is used to identify tool name (`arr` for Arriba fusion transcript call output file prefix, for instance). 

For example:
> ``` wrangle-ft-tsv.py data/FTmyBRCAs_raw/Arriba arr ```


```python
import polars as pl
```


```python
# load up Arriba and FusionCatcher merged dataframes lazily
arriba_mdf = pl.scan_parquet('data/Arriba-fusiontranscript-raw-list.parquet')
fc_mdf = pl.scan_parquet('data/FusionCatcher-fusiontranscript-raw-list.parquet')
```


```python
arriba_mdf.collect()
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (49_465, 11)</small><table border="1" class="dataframe"><thead><tr><th>fusionTranscriptID</th><th>fusionGeneID</th><th>breakpointPair</th><th>strand1</th><th>strand2</th><th>site1</th><th>site2</th><th>type</th><th>confidence</th><th>sampleID</th><th>toolID</th></tr><tr><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>i64</td><td>cat</td></tr></thead><tbody><tr><td>&quot;TRMT11::SMG6__6:125986622-17:2…</td><td>&quot;TRMT11::SMG6&quot;</td><td>&quot;6:125986622-17:2244719&quot;</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;translocation&quot;</td><td>&quot;high&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STAG3::MEF2C-AS1__7:100189570-…</td><td>&quot;STAG3::MEF2C-AS1&quot;</td><td>&quot;7:100189570-5:88919251&quot;</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation/5&#x27;-5&#x27;&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAPK13::C1QL1__6:36132629-17:4…</td><td>&quot;MAPK13::C1QL1&quot;</td><td>&quot;6:36132629-17:44965446&quot;</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation/5&#x27;-5&#x27;&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1__20:58673711-20:…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20:58691724&quot;</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;5&#x27;UTR/splice-site&quot;</td><td>&quot;deletion/read-through&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAPK13::NMT1__6:36132629-17:44…</td><td>&quot;MAPK13::NMT1&quot;</td><td>&quot;6:36132629-17:44965446&quot;</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;DENND5B::AC087311.1(22711),SYT…</td><td>&quot;DENND5B::AC087311.1(22711),SYT…</td><td>&quot;12:31479608-12:33016465&quot;</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intergenic&quot;</td><td>&quot;inversion&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;LINC01145::AC245100.2__1:14520…</td><td>&quot;LINC01145::AC245100.2&quot;</td><td>&quot;1:145201150-1:148436753&quot;</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;exon&quot;</td><td>&quot;exon&quot;</td><td>&quot;duplication/5&#x27;-5&#x27;&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;NET1::RNF169__10:5412820-11:74…</td><td>&quot;NET1::RNF169&quot;</td><td>&quot;10:5412820-11:74834676&quot;</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;translocation&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAN2C1::SIN3A__15:75366522-15:…</td><td>&quot;MAN2C1::SIN3A&quot;</td><td>&quot;15:75366522-15:75375872&quot;</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;duplication&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;LINC02224::MRPS30-DT__5:446584…</td><td>&quot;LINC02224::MRPS30-DT&quot;</td><td>&quot;5:44658462-5:44777328&quot;</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;exon/splice-site&quot;</td><td>&quot;exon/splice-site&quot;</td><td>&quot;duplication&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr></tbody></table></div>




```python
print(fc_mdf.collect())
```

    shape: (31_364, 11)
    ┌─────────────┬─────────────┬────────────┬─────────┬───┬──────┬────────────┬──────────┬────────────┐
    │ fusionTrans ┆ fusionGeneI ┆ breakpoint ┆ strand1 ┆ … ┆ type ┆ confidence ┆ sampleID ┆ toolID     │
    │ criptID     ┆ D           ┆ Pair       ┆ ---     ┆   ┆ ---  ┆ ---        ┆ ---      ┆ ---        │
    │ ---         ┆ ---         ┆ ---        ┆ cat     ┆   ┆ cat  ┆ cat        ┆ i64      ┆ cat        │
    │ cat         ┆ cat         ┆ cat        ┆         ┆   ┆      ┆            ┆          ┆            │
    ╞═════════════╪═════════════╪════════════╪═════════╪═══╪══════╪════════════╪══════════╪════════════╡
    │ SIDT2::TAGL ┆ SIDT2::TAGL ┆ 11:1171959 ┆ +       ┆ … ┆ .    ┆ .          ┆ 2        ┆ FusionCatc │
    │ N__11:11719 ┆ N           ┆ 15-11:1172 ┆         ┆   ┆      ┆            ┆          ┆ her        │
    │ 5915-11:…   ┆             ┆ 03002      ┆         ┆   ┆      ┆            ┆          ┆            │
    │ AZGP1::GJC3 ┆ AZGP1::GJC3 ┆ 7:99971746 ┆ -       ┆ … ┆ .    ┆ .          ┆ 2        ┆ FusionCatc │
    │ __7:9997174 ┆             ┆ -7:9992360 ┆         ┆   ┆      ┆            ┆          ┆ her        │
    │ 6-7:9992…   ┆             ┆ 3          ┆         ┆   ┆      ┆            ┆          ┆            │
    │ NPEPPS::TBC ┆ NPEPPS::TBC ┆ 17:4759254 ┆ +       ┆ … ┆ .    ┆ .          ┆ 2        ┆ FusionCatc │
    │ 1D3__17:475 ┆ 1D3         ┆ 5-17:38191 ┆         ┆   ┆      ┆            ┆          ┆ her        │
    │ 92545-17…   ┆             ┆ 030        ┆         ┆   ┆      ┆            ┆          ┆            │
    │ CYP4F11::CY ┆ CYP4F11::CY ┆ 19:1591476 ┆ -       ┆ … ┆ .    ┆ .          ┆ 2        ┆ FusionCatc │
    │ P4F23P__19: ┆ P4F23P      ┆ 2-19:15583 ┆         ┆   ┆      ┆            ┆          ┆ her        │
    │ 15914762…   ┆             ┆ 580        ┆         ┆   ┆      ┆            ┆          ┆            │
    │ SLC49A3::AT ┆ SLC49A3::AT ┆ 4:686532-4 ┆ -       ┆ … ┆ .    ┆ .          ┆ 3        ┆ FusionCatc │
    │ P5ME__4:686 ┆ P5ME        ┆ :673401    ┆         ┆   ┆      ┆            ┆          ┆ her        │
    │ 532-4:67…   ┆             ┆            ┆         ┆   ┆      ┆            ┆          ┆            │
    │ …           ┆ …           ┆ …          ┆ …       ┆ … ┆ …    ┆ …          ┆ …        ┆ …          │
    │ CTBS::GNG5_ ┆ CTBS::GNG5  ┆ 1:84563257 ┆ -       ┆ … ┆ .    ┆ .          ┆ 991      ┆ FusionCatc │
    │ _1:84563257 ┆             ┆ -1:8450197 ┆         ┆   ┆      ┆            ┆          ┆ her        │
    │ -1:84501…   ┆             ┆ 0          ┆         ┆   ┆      ┆            ┆          ┆            │
    │ MRPS30-DT:: ┆ MRPS30-DT:: ┆ 5:44808642 ┆ -       ┆ … ┆ .    ┆ .          ┆ 991      ┆ FusionCatc │
    │ LINC02224__ ┆ LINC02224   ┆ -5:4465855 ┆         ┆   ┆      ┆            ┆          ┆ her        │
    │ 5:448086…   ┆             ┆ 7          ┆         ┆   ┆      ┆            ┆          ┆            │
    │ NBEA::CR382 ┆ NBEA::CR382 ┆ 13:3507085 ┆ +       ┆ … ┆ .    ┆ .          ┆ 991      ┆ FusionCatc │
    │ 287.1__13:3 ┆ 287.1       ┆ 2-21:10127 ┆         ┆   ┆      ┆            ┆          ┆ her        │
    │ 5070852-…   ┆             ┆ 330        ┆         ┆   ┆      ┆            ┆          ┆            │
    │ HACL1::COLQ ┆ HACL1::COLQ ┆ 3:15563358 ┆ -       ┆ … ┆ .    ┆ .          ┆ 991      ┆ FusionCatc │
    │ __3:1556335 ┆             ┆ -3:1548963 ┆         ┆   ┆      ┆            ┆          ┆ her        │
    │ 8-3:1548…   ┆             ┆ 7          ┆         ┆   ┆      ┆            ┆          ┆            │
    │ PRCP::TASOR ┆ PRCP::TASOR ┆ 11:8290023 ┆ -       ┆ … ┆ .    ┆ .          ┆ 992      ┆ FusionCatc │
    │ 2__11:82900 ┆ 2           ┆ 5-10:57396 ┆         ┆   ┆      ┆            ┆          ┆ her        │
    │ 235-10:5…   ┆             ┆ 18         ┆         ┆   ┆      ┆            ┆          ┆            │
    └─────────────┴─────────────┴────────────┴─────────┴───┴──────┴────────────┴──────────┴────────────┘


Now, we can merge the two dataframes into one masterFrame. Use Polars' `concat` (vertical concatenation is the default, where two dataframes sharing the exact same columns would be joined together, adding all rows of dataframe 1 and 2 vertically).


```python
joined_df = pl.concat(
    [
        arriba_mdf.collect(),
        fc_mdf.collect()
    ]
)

joined_df
```

    sys:1: CategoricalRemappingWarning: Local categoricals have different encodings, expensive re-encoding is done to perform this merge operation. Consider using a StringCache or an Enum type if the categories are known in advance





<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (80_829, 11)</small><table border="1" class="dataframe"><thead><tr><th>fusionTranscriptID</th><th>fusionGeneID</th><th>breakpointPair</th><th>strand1</th><th>strand2</th><th>site1</th><th>site2</th><th>type</th><th>confidence</th><th>sampleID</th><th>toolID</th></tr><tr><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>i64</td><td>cat</td></tr></thead><tbody><tr><td>&quot;TRMT11::SMG6__6:125986622-17:2…</td><td>&quot;TRMT11::SMG6&quot;</td><td>&quot;6:125986622-17:2244719&quot;</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;translocation&quot;</td><td>&quot;high&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STAG3::MEF2C-AS1__7:100189570-…</td><td>&quot;STAG3::MEF2C-AS1&quot;</td><td>&quot;7:100189570-5:88919251&quot;</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation/5&#x27;-5&#x27;&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAPK13::C1QL1__6:36132629-17:4…</td><td>&quot;MAPK13::C1QL1&quot;</td><td>&quot;6:36132629-17:44965446&quot;</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation/5&#x27;-5&#x27;&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1__20:58673711-20:…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20:58691724&quot;</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;5&#x27;UTR/splice-site&quot;</td><td>&quot;deletion/read-through&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAPK13::NMT1__6:36132629-17:44…</td><td>&quot;MAPK13::NMT1&quot;</td><td>&quot;6:36132629-17:44965446&quot;</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;CTBS::GNG5__1:84563257-1:84501…</td><td>&quot;CTBS::GNG5&quot;</td><td>&quot;1:84563257-1:84501970&quot;</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;in-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;MRPS30-DT::LINC02224__5:448086…</td><td>&quot;MRPS30-DT::LINC02224&quot;</td><td>&quot;5:44808642-5:44658557&quot;</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;exonic(no-known-CDS)&quot;</td><td>&quot;exonic(no-known-CDS)&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;NBEA::CR382287.1__13:35070852-…</td><td>&quot;NBEA::CR382287.1&quot;</td><td>&quot;13:35070852-21:10127330&quot;</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS(truncated)&quot;</td><td>&quot;exonic(no-known-CDS)&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;HACL1::COLQ__3:15563358-3:1548…</td><td>&quot;HACL1::COLQ&quot;</td><td>&quot;3:15563358-3:15489637&quot;</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;PRCP::TASOR2__11:82900235-10:5…</td><td>&quot;PRCP::TASOR2&quot;</td><td>&quot;11:82900235-10:5739618&quot;</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>992</td><td>&quot;FusionCatcher&quot;</td></tr></tbody></table></div>



Now sort `sampleID` in ascending order.


```python
joined_df.sort("sampleID")
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (80_829, 11)</small><table border="1" class="dataframe"><thead><tr><th>fusionTranscriptID</th><th>fusionGeneID</th><th>breakpointPair</th><th>strand1</th><th>strand2</th><th>site1</th><th>site2</th><th>type</th><th>confidence</th><th>sampleID</th><th>toolID</th></tr><tr><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>i64</td><td>cat</td></tr></thead><tbody><tr><td>&quot;TRMT11::SMG6__6:125986622-17:2…</td><td>&quot;TRMT11::SMG6&quot;</td><td>&quot;6:125986622-17:2244719&quot;</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;translocation&quot;</td><td>&quot;high&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STAG3::MEF2C-AS1__7:100189570-…</td><td>&quot;STAG3::MEF2C-AS1&quot;</td><td>&quot;7:100189570-5:88919251&quot;</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation/5&#x27;-5&#x27;&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAPK13::C1QL1__6:36132629-17:4…</td><td>&quot;MAPK13::C1QL1&quot;</td><td>&quot;6:36132629-17:44965446&quot;</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation/5&#x27;-5&#x27;&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1__20:58673711-20:…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20:58691724&quot;</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;5&#x27;UTR/splice-site&quot;</td><td>&quot;deletion/read-through&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAPK13::NMT1__6:36132629-17:44…</td><td>&quot;MAPK13::NMT1&quot;</td><td>&quot;6:36132629-17:44965446&quot;</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;LINC01145::AC245100.2__1:14520…</td><td>&quot;LINC01145::AC245100.2&quot;</td><td>&quot;1:145201150-1:148436753&quot;</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;exon&quot;</td><td>&quot;exon&quot;</td><td>&quot;duplication/5&#x27;-5&#x27;&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;NET1::RNF169__10:5412820-11:74…</td><td>&quot;NET1::RNF169&quot;</td><td>&quot;10:5412820-11:74834676&quot;</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;translocation&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAN2C1::SIN3A__15:75366522-15:…</td><td>&quot;MAN2C1::SIN3A&quot;</td><td>&quot;15:75366522-15:75375872&quot;</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;CDS/splice-site&quot;</td><td>&quot;duplication&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;LINC02224::MRPS30-DT__5:446584…</td><td>&quot;LINC02224::MRPS30-DT&quot;</td><td>&quot;5:44658462-5:44777328&quot;</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;exon/splice-site&quot;</td><td>&quot;exon/splice-site&quot;</td><td>&quot;duplication&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PRCP::TASOR2__11:82900235-10:5…</td><td>&quot;PRCP::TASOR2&quot;</td><td>&quot;11:82900235-10:5739618&quot;</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>992</td><td>&quot;FusionCatcher&quot;</td></tr></tbody></table></div>




```python
all_ft_counts = joined_df.select(pl.col("fusionTranscriptID").value_counts(sort=True))
all_ft_counts.unnest("fusionTranscriptID")
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (46_254, 2)</small><table border="1" class="dataframe"><thead><tr><th>fusionTranscriptID</th><th>count</th></tr><tr><td>cat</td><td>u32</td></tr></thead><tbody><tr><td>&quot;CTBS::GNG5__1:84563257-1:84501…</td><td>661</td></tr><tr><td>&quot;AZGP1::GJC3__7:99971746-7:9992…</td><td>608</td></tr><tr><td>&quot;NPEPPS::TBC1D3__17:47592545-17…</td><td>608</td></tr><tr><td>&quot;TMED7::TICAM2__5:115616318-5:1…</td><td>428</td></tr><tr><td>&quot;SIDT2::TAGLN__11:117195915-11:…</td><td>412</td></tr><tr><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;AL021546.1::DYNLL1__12:1204571…</td><td>1</td></tr><tr><td>&quot;CNOT1::C16ORF78__16:58629728-1…</td><td>1</td></tr><tr><td>&quot;RPS19::AC067930.9__19:41872769…</td><td>1</td></tr><tr><td>&quot;GTPBP3::AC097717.1__19:1733959…</td><td>1</td></tr><tr><td>&quot;COX4I1::AC026954.2__16:8580127…</td><td>1</td></tr></tbody></table></div>




```python
genelevel_ft_counts = joined_df.select(pl.col("fusionGeneID").value_counts(sort=True))
genelevel_ft_counts.unnest("fusionGeneID")
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (36_460, 2)</small><table border="1" class="dataframe"><thead><tr><th>fusionGeneID</th><th>count</th></tr><tr><td>cat</td><td>u32</td></tr></thead><tbody><tr><td>&quot;TVP23C::CDRT4&quot;</td><td>1448</td></tr><tr><td>&quot;RBM14::RBM4&quot;</td><td>835</td></tr><tr><td>&quot;AZGP1::GJC3&quot;</td><td>830</td></tr><tr><td>&quot;SMG1::NPIPB5&quot;</td><td>765</td></tr><tr><td>&quot;CTBS::GNG5&quot;</td><td>663</td></tr><tr><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;AL021546.1::DYNLL1&quot;</td><td>1</td></tr><tr><td>&quot;CNOT1::C16ORF78&quot;</td><td>1</td></tr><tr><td>&quot;RPS19::AC067930.9&quot;</td><td>1</td></tr><tr><td>&quot;GTPBP3::AC097717.1&quot;</td><td>1</td></tr><tr><td>&quot;COX4I1::AC026954.2&quot;</td><td>1</td></tr></tbody></table></div>




```python

```
