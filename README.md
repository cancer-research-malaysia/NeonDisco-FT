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
arriba_mdf = pl.scan_parquet('output/myBRCA/Arriba-fusiontranscript-raw-list.parquet')
fc_mdf = pl.scan_parquet('output/myBRCA/FusionCatcher-fusiontranscript-raw-list.parquet')
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
<small>shape: (49_465, 11)</small><table border="1" class="dataframe"><thead><tr><th>fusionTranscriptID</th><th>fusionGeneID</th><th>breakpointPair</th><th>strand1</th><th>strand2</th><th>site1</th><th>site2</th><th>type</th><th>confidence</th><th>sampleID</th><th>toolID</th></tr><tr><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>i64</td><td>cat</td></tr></thead><tbody><tr><td>&quot;TRMT11::SMG6__…</td><td>&quot;TRMT11::SMG6&quot;</td><td>&quot;6:125986622-17…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;translocation&quot;</td><td>&quot;high&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STAG3::MEF2C-A…</td><td>&quot;STAG3::MEF2C-A…</td><td>&quot;7:100189570-5:…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation/…</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAPK13::C1QL1_…</td><td>&quot;MAPK13::C1QL1&quot;</td><td>&quot;6:36132629-17:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation/…</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1_…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAPK13::NMT1__…</td><td>&quot;MAPK13::NMT1&quot;</td><td>&quot;6:36132629-17:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;TMEM63A::SRP9(…</td><td>&quot;TMEM63A::SRP9(…</td><td>&quot;1:225853629-1:…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;intergenic&quot;</td><td>&quot;deletion/read-…</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCMT1::AL13609…</td><td>&quot;PCMT1::AL13609…</td><td>&quot;6:149773169-6:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;intergenic&quot;</td><td>&quot;deletion&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;RCOR3::AL59104…</td><td>&quot;RCOR3::AL59104…</td><td>&quot;1:211274262-1:…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;exon&quot;</td><td>&quot;inversion/3&#x27;-3…</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;TRAF5::C4BPB__…</td><td>&quot;TRAF5::C4BPB&quot;</td><td>&quot;1:211326889-1:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;CDS/splice-sit…</td><td>&quot;duplication&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;SIPA1L2::FMN2_…</td><td>&quot;SIPA1L2::FMN2&quot;</td><td>&quot;1:232415494-1:…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;inversion&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;ATP5ME::GTF2IR…</td><td>&quot;ATP5ME::GTF2IR…</td><td>&quot;4:673303-7:748…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;translocation&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PRPF4B::SLC22A…</td><td>&quot;PRPF4B::SLC22A…</td><td>&quot;6:4037570-6:35…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;intergenic&quot;</td><td>&quot;duplication&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;5:141357621-5:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;5:141357621-5:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;5:141357621-5:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;5:141357621-5:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;5:141357621-5:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;RSL24D1::RNU6-…</td><td>&quot;RSL24D1::RNU6-…</td><td>&quot;15:55196840-10…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intergenic&quot;</td><td>&quot;translocation&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;AL353899.1(156…</td><td>&quot;AL353899.1(156…</td><td>&quot;1:157453420-1:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;intergenic&quot;</td><td>&quot;3&#x27;UTR&quot;</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;DENND5B::AC087…</td><td>&quot;DENND5B::AC087…</td><td>&quot;12:31479608-12…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intergenic&quot;</td><td>&quot;inversion&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;LINC01145::AC2…</td><td>&quot;LINC01145::AC2…</td><td>&quot;1:145201150-1:…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;exon&quot;</td><td>&quot;exon&quot;</td><td>&quot;duplication/5&#x27;…</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;NET1::RNF169__…</td><td>&quot;NET1::RNF169&quot;</td><td>&quot;10:5412820-11:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;translocation&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAN2C1::SIN3A_…</td><td>&quot;MAN2C1::SIN3A&quot;</td><td>&quot;15:75366522-15…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;duplication&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;LINC02224::MRP…</td><td>&quot;LINC02224::MRP…</td><td>&quot;5:44658462-5:4…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;exon/splice-si…</td><td>&quot;exon/splice-si…</td><td>&quot;duplication&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr></tbody></table></div>




```python
fc_mdf.collect()
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (31_364, 11)</small><table border="1" class="dataframe"><thead><tr><th>fusionTranscriptID</th><th>fusionGeneID</th><th>breakpointPair</th><th>strand1</th><th>strand2</th><th>site1</th><th>site2</th><th>type</th><th>confidence</th><th>sampleID</th><th>toolID</th></tr><tr><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>i64</td><td>cat</td></tr></thead><tbody><tr><td>&quot;SIDT2::TAGLN__…</td><td>&quot;SIDT2::TAGLN&quot;</td><td>&quot;11:117195915-1…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;UTR&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>2</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;AZGP1::GJC3__7…</td><td>&quot;AZGP1::GJC3&quot;</td><td>&quot;7:99971746-7:9…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;in-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>2</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;NPEPPS::TBC1D3…</td><td>&quot;NPEPPS::TBC1D3…</td><td>&quot;17:47592545-17…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS(complete)&quot;</td><td>&quot;UTR&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>2</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;CYP4F11::CYP4F…</td><td>&quot;CYP4F11::CYP4F…</td><td>&quot;19:15914762-19…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>2</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;SLC49A3::ATP5M…</td><td>&quot;SLC49A3::ATP5M…</td><td>&quot;4:686532-4:673…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>3</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;AZGP1::GJC3__7…</td><td>&quot;AZGP1::GJC3&quot;</td><td>&quot;7:99971746-7:9…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;in-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>4</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;MALRD1::PLXDC2…</td><td>&quot;MALRD1::PLXDC2…</td><td>&quot;10:19790186-10…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;UTR&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>4</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;AC004066.1::PP…</td><td>&quot;AC004066.1::PP…</td><td>&quot;4:105552193-4:…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;exonic(no-know…</td><td>&quot;CDS(truncated)…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>4</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;NPEPPS::TBC1D3…</td><td>&quot;NPEPPS::TBC1D3…</td><td>&quot;17:47592545-17…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS(complete)&quot;</td><td>&quot;UTR&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>4</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;PRH1::AC007450…</td><td>&quot;PRH1::AC007450…</td><td>&quot;12:11171422-12…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;UTR&quot;</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>5</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;AZGP1::GJC3__7…</td><td>&quot;AZGP1::GJC3&quot;</td><td>&quot;7:99971746-7:9…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;in-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>5</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;CTBS::GNG5__1:…</td><td>&quot;CTBS::GNG5&quot;</td><td>&quot;1:84563257-1:8…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;in-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>5</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;RPS19::AC06793…</td><td>&quot;RPS19::AC06793…</td><td>&quot;19:41872769-8:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;UTR&quot;</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>990</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;GTPBP3::AC0977…</td><td>&quot;GTPBP3::AC0977…</td><td>&quot;19:17339599-2:…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>990</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;KANSL1::ARL17B…</td><td>&quot;KANSL1::ARL17B…</td><td>&quot;17:46094560-17…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>990</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;KANSL1::ARL17A…</td><td>&quot;KANSL1::ARL17A…</td><td>&quot;17:46094560-17…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>990</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;COX4I1::AC0269…</td><td>&quot;COX4I1::AC0269…</td><td>&quot;16:85801278-17…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;intronic&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;KMT5B::METTL15…</td><td>&quot;KMT5B::METTL15…</td><td>&quot;11:68213477-11…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;UTR&quot;</td><td>&quot;UTR&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;NPEPPS::TBC1D3…</td><td>&quot;NPEPPS::TBC1D3…</td><td>&quot;17:47592545-17…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS(complete)&quot;</td><td>&quot;UTR&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;CTBS::GNG5__1:…</td><td>&quot;CTBS::GNG5&quot;</td><td>&quot;1:84563257-1:8…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;in-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;MRPS30-DT::LIN…</td><td>&quot;MRPS30-DT::LIN…</td><td>&quot;5:44808642-5:4…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;exonic(no-know…</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;NBEA::CR382287…</td><td>&quot;NBEA::CR382287…</td><td>&quot;13:35070852-21…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;HACL1::COLQ__3…</td><td>&quot;HACL1::COLQ&quot;</td><td>&quot;3:15563358-3:1…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;PRCP::TASOR2__…</td><td>&quot;PRCP::TASOR2&quot;</td><td>&quot;11:82900235-10…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>992</td><td>&quot;FusionCatcher&quot;</td></tr></tbody></table></div>



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
<small>shape: (80_829, 11)</small><table border="1" class="dataframe"><thead><tr><th>fusionTranscriptID</th><th>fusionGeneID</th><th>breakpointPair</th><th>strand1</th><th>strand2</th><th>site1</th><th>site2</th><th>type</th><th>confidence</th><th>sampleID</th><th>toolID</th></tr><tr><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>i64</td><td>cat</td></tr></thead><tbody><tr><td>&quot;TRMT11::SMG6__…</td><td>&quot;TRMT11::SMG6&quot;</td><td>&quot;6:125986622-17…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;translocation&quot;</td><td>&quot;high&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STAG3::MEF2C-A…</td><td>&quot;STAG3::MEF2C-A…</td><td>&quot;7:100189570-5:…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation/…</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAPK13::C1QL1_…</td><td>&quot;MAPK13::C1QL1&quot;</td><td>&quot;6:36132629-17:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation/…</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1_…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAPK13::NMT1__…</td><td>&quot;MAPK13::NMT1&quot;</td><td>&quot;6:36132629-17:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;TMEM63A::SRP9(…</td><td>&quot;TMEM63A::SRP9(…</td><td>&quot;1:225853629-1:…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;intergenic&quot;</td><td>&quot;deletion/read-…</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCMT1::AL13609…</td><td>&quot;PCMT1::AL13609…</td><td>&quot;6:149773169-6:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;intergenic&quot;</td><td>&quot;deletion&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;RCOR3::AL59104…</td><td>&quot;RCOR3::AL59104…</td><td>&quot;1:211274262-1:…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;exon&quot;</td><td>&quot;inversion/3&#x27;-3…</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;TRAF5::C4BPB__…</td><td>&quot;TRAF5::C4BPB&quot;</td><td>&quot;1:211326889-1:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;CDS/splice-sit…</td><td>&quot;duplication&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;SIPA1L2::FMN2_…</td><td>&quot;SIPA1L2::FMN2&quot;</td><td>&quot;1:232415494-1:…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;inversion&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;ATP5ME::GTF2IR…</td><td>&quot;ATP5ME::GTF2IR…</td><td>&quot;4:673303-7:748…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;translocation&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PRPF4B::SLC22A…</td><td>&quot;PRPF4B::SLC22A…</td><td>&quot;6:4037570-6:35…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;intergenic&quot;</td><td>&quot;duplication&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;RPS19::AC06793…</td><td>&quot;RPS19::AC06793…</td><td>&quot;19:41872769-8:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;UTR&quot;</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>990</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;GTPBP3::AC0977…</td><td>&quot;GTPBP3::AC0977…</td><td>&quot;19:17339599-2:…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>990</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;KANSL1::ARL17B…</td><td>&quot;KANSL1::ARL17B…</td><td>&quot;17:46094560-17…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>990</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;KANSL1::ARL17A…</td><td>&quot;KANSL1::ARL17A…</td><td>&quot;17:46094560-17…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>990</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;COX4I1::AC0269…</td><td>&quot;COX4I1::AC0269…</td><td>&quot;16:85801278-17…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;intronic&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;KMT5B::METTL15…</td><td>&quot;KMT5B::METTL15…</td><td>&quot;11:68213477-11…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;UTR&quot;</td><td>&quot;UTR&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;NPEPPS::TBC1D3…</td><td>&quot;NPEPPS::TBC1D3…</td><td>&quot;17:47592545-17…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS(complete)&quot;</td><td>&quot;UTR&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;CTBS::GNG5__1:…</td><td>&quot;CTBS::GNG5&quot;</td><td>&quot;1:84563257-1:8…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;in-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;MRPS30-DT::LIN…</td><td>&quot;MRPS30-DT::LIN…</td><td>&quot;5:44808642-5:4…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;exonic(no-know…</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;NBEA::CR382287…</td><td>&quot;NBEA::CR382287…</td><td>&quot;13:35070852-21…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;HACL1::COLQ__3…</td><td>&quot;HACL1::COLQ&quot;</td><td>&quot;3:15563358-3:1…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>991</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;PRCP::TASOR2__…</td><td>&quot;PRCP::TASOR2&quot;</td><td>&quot;11:82900235-10…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>992</td><td>&quot;FusionCatcher&quot;</td></tr></tbody></table></div>



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
<small>shape: (80_829, 11)</small><table border="1" class="dataframe"><thead><tr><th>fusionTranscriptID</th><th>fusionGeneID</th><th>breakpointPair</th><th>strand1</th><th>strand2</th><th>site1</th><th>site2</th><th>type</th><th>confidence</th><th>sampleID</th><th>toolID</th></tr><tr><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>i64</td><td>cat</td></tr></thead><tbody><tr><td>&quot;TRMT11::SMG6__…</td><td>&quot;TRMT11::SMG6&quot;</td><td>&quot;6:125986622-17…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;translocation&quot;</td><td>&quot;high&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STAG3::MEF2C-A…</td><td>&quot;STAG3::MEF2C-A…</td><td>&quot;7:100189570-5:…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation/…</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAPK13::C1QL1_…</td><td>&quot;MAPK13::C1QL1&quot;</td><td>&quot;6:36132629-17:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation/…</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1_…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAPK13::NMT1__…</td><td>&quot;MAPK13::NMT1&quot;</td><td>&quot;6:36132629-17:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;TMEM63A::SRP9(…</td><td>&quot;TMEM63A::SRP9(…</td><td>&quot;1:225853629-1:…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;intergenic&quot;</td><td>&quot;deletion/read-…</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCMT1::AL13609…</td><td>&quot;PCMT1::AL13609…</td><td>&quot;6:149773169-6:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;intergenic&quot;</td><td>&quot;deletion&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;RCOR3::AL59104…</td><td>&quot;RCOR3::AL59104…</td><td>&quot;1:211274262-1:…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;exon&quot;</td><td>&quot;inversion/3&#x27;-3…</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;TRAF5::C4BPB__…</td><td>&quot;TRAF5::C4BPB&quot;</td><td>&quot;1:211326889-1:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;CDS/splice-sit…</td><td>&quot;duplication&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;SIPA1L2::FMN2_…</td><td>&quot;SIPA1L2::FMN2&quot;</td><td>&quot;1:232415494-1:…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;inversion&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;ATP5ME::GTF2IR…</td><td>&quot;ATP5ME::GTF2IR…</td><td>&quot;4:673303-7:748…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;translocation&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PRPF4B::SLC22A…</td><td>&quot;PRPF4B::SLC22A…</td><td>&quot;6:4037570-6:35…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;intergenic&quot;</td><td>&quot;duplication&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;5:141357621-5:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;5:141357621-5:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;5:141357621-5:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;PCDHGA4::PCDHG…</td><td>&quot;5:141357621-5:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;RSL24D1::RNU6-…</td><td>&quot;RSL24D1::RNU6-…</td><td>&quot;15:55196840-10…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intergenic&quot;</td><td>&quot;translocation&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;AL353899.1(156…</td><td>&quot;AL353899.1(156…</td><td>&quot;1:157453420-1:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;intergenic&quot;</td><td>&quot;3&#x27;UTR&quot;</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;DENND5B::AC087…</td><td>&quot;DENND5B::AC087…</td><td>&quot;12:31479608-12…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intergenic&quot;</td><td>&quot;inversion&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;LINC01145::AC2…</td><td>&quot;LINC01145::AC2…</td><td>&quot;1:145201150-1:…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;exon&quot;</td><td>&quot;exon&quot;</td><td>&quot;duplication/5&#x27;…</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;NET1::RNF169__…</td><td>&quot;NET1::RNF169&quot;</td><td>&quot;10:5412820-11:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;translocation&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAN2C1::SIN3A_…</td><td>&quot;MAN2C1::SIN3A&quot;</td><td>&quot;15:75366522-15…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;duplication&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;LINC02224::MRP…</td><td>&quot;LINC02224::MRP…</td><td>&quot;5:44658462-5:4…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;exon/splice-si…</td><td>&quot;exon/splice-si…</td><td>&quot;duplication&quot;</td><td>&quot;low&quot;</td><td>992</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PRCP::TASOR2__…</td><td>&quot;PRCP::TASOR2&quot;</td><td>&quot;11:82900235-10…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>992</td><td>&quot;FusionCatcher&quot;</td></tr></tbody></table></div>




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
<small>shape: (46_254, 2)</small><table border="1" class="dataframe"><thead><tr><th>fusionTranscriptID</th><th>count</th></tr><tr><td>cat</td><td>u32</td></tr></thead><tbody><tr><td>&quot;CTBS::GNG5__1:…</td><td>661</td></tr><tr><td>&quot;AZGP1::GJC3__7…</td><td>608</td></tr><tr><td>&quot;NPEPPS::TBC1D3…</td><td>608</td></tr><tr><td>&quot;TMED7::TICAM2_…</td><td>428</td></tr><tr><td>&quot;SIDT2::TAGLN__…</td><td>412</td></tr><tr><td>&quot;TVP23C::CDRT4_…</td><td>400</td></tr><tr><td>&quot;RBM14::RBM4__1…</td><td>328</td></tr><tr><td>&quot;AC092807.3::DD…</td><td>308</td></tr><tr><td>&quot;RBM14::RBM4__1…</td><td>304</td></tr><tr><td>&quot;STX16::NPEPL1_…</td><td>291</td></tr><tr><td>&quot;SMG1::NPIPB5__…</td><td>286</td></tr><tr><td>&quot;TVP23C::CDRT4_…</td><td>272</td></tr><tr><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;TACC2::C10ORF1…</td><td>1</td></tr><tr><td>&quot;TACC2::C10ORF1…</td><td>1</td></tr><tr><td>&quot;PPP1R1B::EIF3D…</td><td>1</td></tr><tr><td>&quot;SDC1::GPAT3__2…</td><td>1</td></tr><tr><td>&quot;SPOP::FLJ40194…</td><td>1</td></tr><tr><td>&quot;SPOP::FLJ40194…</td><td>1</td></tr><tr><td>&quot;BPTF::GNA13__1…</td><td>1</td></tr><tr><td>&quot;AL021546.1::DY…</td><td>1</td></tr><tr><td>&quot;CNOT1::C16ORF7…</td><td>1</td></tr><tr><td>&quot;RPS19::AC06793…</td><td>1</td></tr><tr><td>&quot;GTPBP3::AC0977…</td><td>1</td></tr><tr><td>&quot;COX4I1::AC0269…</td><td>1</td></tr></tbody></table></div>



# Filter Out non-TNBCs

We can filter out rows corresponding to `sampleID` more than 172, because these are not TNBC samples.




```python
tnbc_only_df = joined_df.filter(pl.col("sampleID") < 173)
```


```python
tnbc_only_df
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (15_264, 11)</small><table border="1" class="dataframe"><thead><tr><th>fusionTranscriptID</th><th>fusionGeneID</th><th>breakpointPair</th><th>strand1</th><th>strand2</th><th>site1</th><th>site2</th><th>type</th><th>confidence</th><th>sampleID</th><th>toolID</th></tr><tr><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>i64</td><td>cat</td></tr></thead><tbody><tr><td>&quot;TRMT11::SMG6__…</td><td>&quot;TRMT11::SMG6&quot;</td><td>&quot;6:125986622-17…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;translocation&quot;</td><td>&quot;high&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STAG3::MEF2C-A…</td><td>&quot;STAG3::MEF2C-A…</td><td>&quot;7:100189570-5:…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation/…</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAPK13::C1QL1_…</td><td>&quot;MAPK13::C1QL1&quot;</td><td>&quot;6:36132629-17:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation/…</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1_…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAPK13::NMT1__…</td><td>&quot;MAPK13::NMT1&quot;</td><td>&quot;6:36132629-17:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;TMEM63A::SRP9(…</td><td>&quot;TMEM63A::SRP9(…</td><td>&quot;1:225853629-1:…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;intergenic&quot;</td><td>&quot;deletion/read-…</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCMT1::AL13609…</td><td>&quot;PCMT1::AL13609…</td><td>&quot;6:149773169-6:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;intergenic&quot;</td><td>&quot;deletion&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;RCOR3::AL59104…</td><td>&quot;RCOR3::AL59104…</td><td>&quot;1:211274262-1:…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;exon&quot;</td><td>&quot;inversion/3&#x27;-3…</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;TRAF5::C4BPB__…</td><td>&quot;TRAF5::C4BPB&quot;</td><td>&quot;1:211326889-1:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;CDS/splice-sit…</td><td>&quot;duplication&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;SIPA1L2::FMN2_…</td><td>&quot;SIPA1L2::FMN2&quot;</td><td>&quot;1:232415494-1:…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;inversion&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;ATP5ME::GTF2IR…</td><td>&quot;ATP5ME::GTF2IR…</td><td>&quot;4:673303-7:748…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;translocation&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PRPF4B::SLC22A…</td><td>&quot;PRPF4B::SLC22A…</td><td>&quot;6:4037570-6:35…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;intergenic&quot;</td><td>&quot;duplication&quot;</td><td>&quot;high&quot;</td><td>2</td><td>&quot;Arriba&quot;</td></tr><tr><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;UTP4::TANGO6__…</td><td>&quot;UTP4::TANGO6&quot;</td><td>&quot;16:69133618-16…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;in-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>171</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;CASC15::CDKAL1…</td><td>&quot;CASC15::CDKAL1…</td><td>&quot;6:21666645-6:2…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;exonic(no-know…</td><td>&quot;CDS(truncated)…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>171</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;MAPK8::RBP3__1…</td><td>&quot;MAPK8::RBP3&quot;</td><td>&quot;10:48431270-10…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>171</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;MFGE8::HAPLN3_…</td><td>&quot;MFGE8::HAPLN3&quot;</td><td>&quot;15:88901551-15…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>171</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;NRDC::LINC0156…</td><td>&quot;NRDC::LINC0156…</td><td>&quot;1:51834017-1:5…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;intronic&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>171</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;NRDC::LINC0156…</td><td>&quot;NRDC::LINC0156…</td><td>&quot;1:51834017-1:5…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>171</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;PBXIP1::PMVK__…</td><td>&quot;PBXIP1::PMVK&quot;</td><td>&quot;1:154945572-1:…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;in-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>171</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;SIDT2::TAGLN__…</td><td>&quot;SIDT2::TAGLN&quot;</td><td>&quot;11:117195915-1…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;UTR&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>171</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;CTNNBIP1::CLST…</td><td>&quot;CTNNBIP1::CLST…</td><td>&quot;1:9871187-1:97…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;in-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>171</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;NPEPPS::TBC1D3…</td><td>&quot;NPEPPS::TBC1D3…</td><td>&quot;17:47592545-17…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS(complete)&quot;</td><td>&quot;UTR&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>171</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;SYT8::TNNI2__1…</td><td>&quot;SYT8::TNNI2&quot;</td><td>&quot;11:1836861-11:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;UTR&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>171</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;AZGP1::GJC3__7…</td><td>&quot;AZGP1::GJC3&quot;</td><td>&quot;7:99971746-7:9…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;in-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>172</td><td>&quot;FusionCatcher&quot;</td></tr></tbody></table></div>



### Counting Sharedness by Breakpoint Pair

To find rows where the `breakpointPair` value is repeated, use Polars' `filter` with a window function


```python
duplicated_df = joined_df.filter(
    pl.col('breakpointPair').count().over('breakpointPair') > 1
).sort('breakpointPair')
```


```python
duplicated_df
```


```python
## do the same to TNBC only df
shared_ft_df = tnbc_only_df.filter(
    pl.col('breakpointPair').count().over('breakpointPair') > 1
).sort('breakpointPair')
```


```python
shared_ft_df
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (7_375, 11)</small><table border="1" class="dataframe"><thead><tr><th>fusionTranscriptID</th><th>fusionGeneID</th><th>breakpointPair</th><th>strand1</th><th>strand2</th><th>site1</th><th>site2</th><th>type</th><th>confidence</th><th>sampleID</th><th>toolID</th></tr><tr><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>i64</td><td>cat</td></tr></thead><tbody><tr><td>&quot;MAPK13::C1QL1_…</td><td>&quot;MAPK13::C1QL1&quot;</td><td>&quot;6:36132629-17:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation/…</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;MAPK13::NMT1__…</td><td>&quot;MAPK13::NMT1&quot;</td><td>&quot;6:36132629-17:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;intron&quot;</td><td>&quot;translocation&quot;</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1_…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>1</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1_…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>4</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1_…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>9</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1_…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;medium&quot;</td><td>10</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1_…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>18</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1_…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>22</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1_…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>37</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1_…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>40</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1_…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>42</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;STX16::NPEPL1_…</td><td>&quot;STX16::NPEPL1&quot;</td><td>&quot;20:58673711-20…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>54</td><td>&quot;Arriba&quot;</td></tr><tr><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;IGKC::GTF2I__2…</td><td>&quot;IGKC::GTF2I&quot;</td><td>&quot;2:88857161-7:7…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;UTR&quot;</td><td>&quot;UTR&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>107</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;IGKC::GTF2I__2…</td><td>&quot;IGKC::GTF2I&quot;</td><td>&quot;2:88857161-7:7…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;UTR&quot;</td><td>&quot;UTR&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>120</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;IGK@::AP3D1__2…</td><td>&quot;IGK@::AP3D1&quot;</td><td>&quot;2:88947305-19:…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;---&quot;</td><td>&quot;intronic&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>107</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;IGK@::AP3D1__2…</td><td>&quot;IGK@::AP3D1&quot;</td><td>&quot;2:88947305-19:…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;---&quot;</td><td>&quot;intronic&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>118</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;SAMD5::SASH1__…</td><td>&quot;SAMD5::SASH1&quot;</td><td>&quot;6:147509387-6:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;in-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>112</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;SAMD5::SASH1__…</td><td>&quot;SAMD5::SASH1&quot;</td><td>&quot;6:147509387-6:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;in-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>130</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;AP003066.1::AC…</td><td>&quot;AP003066.1::AC…</td><td>&quot;11:96590474-7:…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;exonic(no-know…</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>116</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;AP003066.1::AC…</td><td>&quot;AP003066.1::AC…</td><td>&quot;11:96590474-7:…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;exonic(no-know…</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>153</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;IGK@::ATP11B__…</td><td>&quot;IGK@::ATP11B&quot;</td><td>&quot;2:88857243-3:1…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;---&quot;</td><td>&quot;intronic&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>120</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;IGK@::ATP11B__…</td><td>&quot;IGK@::ATP11B&quot;</td><td>&quot;2:88857243-3:1…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;---&quot;</td><td>&quot;intronic&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>126</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;AC007402.2::FG…</td><td>&quot;AC007402.2::FG…</td><td>&quot;2:51442042-1:5…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;exonic(no-know…</td><td>&quot;CDS(truncated)…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>141</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;AC007402.2::FG…</td><td>&quot;AC007402.2::FG…</td><td>&quot;2:51442042-1:5…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;exonic(no-know…</td><td>&quot;CDS(truncated)…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>164</td><td>&quot;FusionCatcher&quot;</td></tr></tbody></table></div>




```python
duplicated_df.write_csv("duplicated-manual.tsv", separator="\t")
```


```python
summary_df = joined_df.group_by('breakpointPair').agg([
    pl.n_unique('sampleID').alias('unique_samples')
]).filter(
    pl.col('unique_samples') > 1
).sort('unique_samples', descending=False)


```


```python
summary_df
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (3_196, 2)</small><table border="1" class="dataframe"><thead><tr><th>breakpointPair</th><th>unique_samples</th></tr><tr><td>cat</td><td>u32</td></tr></thead><tbody><tr><td>&quot;8:69640613-19:…</td><td>2</td></tr><tr><td>&quot;18:50002040-17…</td><td>2</td></tr><tr><td>&quot;6:106970265-2:…</td><td>2</td></tr><tr><td>&quot;6:73524875-20:…</td><td>2</td></tr><tr><td>&quot;2:88944680-2:1…</td><td>2</td></tr><tr><td>&quot;9:120659374-9:…</td><td>2</td></tr><tr><td>&quot;7:103099569-2:…</td><td>2</td></tr><tr><td>&quot;8:98026188-12:…</td><td>2</td></tr><tr><td>&quot;15:75890483-15…</td><td>2</td></tr><tr><td>&quot;12:76074203-12…</td><td>2</td></tr><tr><td>&quot;2:88882339-19:…</td><td>2</td></tr><tr><td>&quot;7:102665693-7:…</td><td>2</td></tr><tr><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;17:15545785-17…</td><td>272</td></tr><tr><td>&quot;16:18858170-16…</td><td>286</td></tr><tr><td>&quot;20:58673711-20…</td><td>291</td></tr><tr><td>&quot;11:66617057-11…</td><td>304</td></tr><tr><td>&quot;1:85496166-1:8…</td><td>308</td></tr><tr><td>&quot;11:66617057-11…</td><td>328</td></tr><tr><td>&quot;17:15503098-17…</td><td>400</td></tr><tr><td>&quot;11:117195915-1…</td><td>412</td></tr><tr><td>&quot;5:115616318-5:…</td><td>428</td></tr><tr><td>&quot;7:99971746-7:9…</td><td>608</td></tr><tr><td>&quot;17:47592545-17…</td><td>608</td></tr><tr><td>&quot;1:84563257-1:8…</td><td>661</td></tr></tbody></table></div>




```python
# rename column breakpointPair to breakpointID
summary_df = summary_df.rename({"breakpointPair": "breakpointID"})
summary_df
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (3_196, 2)</small><table border="1" class="dataframe"><thead><tr><th>breakpointID</th><th>unique_samples</th></tr><tr><td>cat</td><td>u32</td></tr></thead><tbody><tr><td>&quot;8:69640613-19:…</td><td>2</td></tr><tr><td>&quot;18:50002040-17…</td><td>2</td></tr><tr><td>&quot;6:106970265-2:…</td><td>2</td></tr><tr><td>&quot;6:73524875-20:…</td><td>2</td></tr><tr><td>&quot;2:88944680-2:1…</td><td>2</td></tr><tr><td>&quot;9:120659374-9:…</td><td>2</td></tr><tr><td>&quot;7:103099569-2:…</td><td>2</td></tr><tr><td>&quot;8:98026188-12:…</td><td>2</td></tr><tr><td>&quot;15:75890483-15…</td><td>2</td></tr><tr><td>&quot;12:76074203-12…</td><td>2</td></tr><tr><td>&quot;2:88882339-19:…</td><td>2</td></tr><tr><td>&quot;7:102665693-7:…</td><td>2</td></tr><tr><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;17:15545785-17…</td><td>272</td></tr><tr><td>&quot;16:18858170-16…</td><td>286</td></tr><tr><td>&quot;20:58673711-20…</td><td>291</td></tr><tr><td>&quot;11:66617057-11…</td><td>304</td></tr><tr><td>&quot;1:85496166-1:8…</td><td>308</td></tr><tr><td>&quot;11:66617057-11…</td><td>328</td></tr><tr><td>&quot;17:15503098-17…</td><td>400</td></tr><tr><td>&quot;11:117195915-1…</td><td>412</td></tr><tr><td>&quot;5:115616318-5:…</td><td>428</td></tr><tr><td>&quot;7:99971746-7:9…</td><td>608</td></tr><tr><td>&quot;17:47592545-17…</td><td>608</td></tr><tr><td>&quot;1:84563257-1:8…</td><td>661</td></tr></tbody></table></div>




```python
summary_df.write_csv("summary_agg.tsv", separator='\t')
```


```python
# now load TCGANormals

# load up Arriba and FusionCatcher merged dataframes lazily
arriba_norms_mdf = pl.scan_parquet('output/TCGANormals/Arriba-TCGANormals-fusiontranscript-raw-list.parquet')
fc_norms_mdf = pl.scan_parquet('output/TCGANormals/FusionCatcher-TCGANormals-fusiontranscript-raw-list.parquet')
```


```python
joined_norms_df = pl.concat(
    [
        arriba_norms_mdf.collect(),
        fc_norms_mdf.collect()
    ]
)

joined_norms_df
```

    sys:1: CategoricalRemappingWarning: Local categoricals have different encodings, expensive re-encoding is done to perform this merge operation. Consider using a StringCache or an Enum type if the categories are known in advance





<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (14_561, 11)</small><table border="1" class="dataframe"><thead><tr><th>fusionTranscriptID</th><th>fusionGeneID</th><th>breakpointPair</th><th>strand1</th><th>strand2</th><th>site1</th><th>site2</th><th>type</th><th>confidence</th><th>sampleID</th><th>toolID</th></tr><tr><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>cat</td><td>str</td><td>cat</td></tr></thead><tbody><tr><td>&quot;TRPM7::SPPL2A_…</td><td>&quot;TRPM7::SPPL2A&quot;</td><td>&quot;15:50686531-15…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;duplication&quot;</td><td>&quot;high&quot;</td><td>&quot;TCGA-A7-A0CE&quot;</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;AC084756.2::SP…</td><td>&quot;AC084756.2::SP…</td><td>&quot;15:50686531-15…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;duplication&quot;</td><td>&quot;high&quot;</td><td>&quot;TCGA-A7-A0CE&quot;</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;BOLA2B::SMG1P5…</td><td>&quot;BOLA2B::SMG1P5…</td><td>&quot;16:30193358-16…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;exon&quot;</td><td>&quot;duplication&quot;</td><td>&quot;high&quot;</td><td>&quot;TCGA-A7-A0CE&quot;</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;FBXO25::SEPTIN…</td><td>&quot;FBXO25::SEPTIN…</td><td>&quot;8:435707-7:557…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;CDS/splice-sit…</td><td>&quot;translocation&quot;</td><td>&quot;high&quot;</td><td>&quot;TCGA-A7-A0CE&quot;</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;TVP23C::CDRT4_…</td><td>&quot;TVP23C::CDRT4&quot;</td><td>&quot;17:15503098-17…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS&quot;</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;medium&quot;</td><td>&quot;TCGA-A7-A0CE&quot;</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;TVP23C::CDRT4_…</td><td>&quot;TVP23C::CDRT4&quot;</td><td>&quot;17:15540433-17…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;5&#x27;UTR/splice-s…</td><td>&quot;deletion/read-…</td><td>&quot;low&quot;</td><td>&quot;TCGA-A7-A0CE&quot;</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCDHGB2::PCDHG…</td><td>&quot;PCDHGB2::PCDHG…</td><td>&quot;5:141362556-5:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;deletion/read-…</td><td>&quot;medium&quot;</td><td>&quot;TCGA-A7-A0CE&quot;</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCDHGB2::PCDHG…</td><td>&quot;PCDHGB2::PCDHG…</td><td>&quot;5:141362556-5:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;deletion/read-…</td><td>&quot;medium&quot;</td><td>&quot;TCGA-A7-A0CE&quot;</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCDHGB2::PCDHG…</td><td>&quot;PCDHGB2::PCDHG…</td><td>&quot;5:141362556-5:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;deletion/read-…</td><td>&quot;medium&quot;</td><td>&quot;TCGA-A7-A0CE&quot;</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCDHGB2::PCDHG…</td><td>&quot;PCDHGB2::PCDHG…</td><td>&quot;5:141362556-5:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;deletion/read-…</td><td>&quot;medium&quot;</td><td>&quot;TCGA-A7-A0CE&quot;</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCDHGB2::PCDHG…</td><td>&quot;PCDHGB2::PCDHG…</td><td>&quot;5:141362556-5:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;deletion/read-…</td><td>&quot;medium&quot;</td><td>&quot;TCGA-A7-A0CE&quot;</td><td>&quot;Arriba&quot;</td></tr><tr><td>&quot;PCDHGB2::PCDHG…</td><td>&quot;PCDHGB2::PCDHG…</td><td>&quot;5:141362556-5:…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS&quot;</td><td>&quot;CDS/splice-sit…</td><td>&quot;deletion/read-…</td><td>&quot;medium&quot;</td><td>&quot;TCGA-A7-A0CE&quot;</td><td>&quot;Arriba&quot;</td></tr><tr><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;NSF::LRRC37A3_…</td><td>&quot;NSF::LRRC37A3&quot;</td><td>&quot;17:46704854-17…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;UTR&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;TCGA-GI-A2C9&quot;</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;SERF1A::SMN1__…</td><td>&quot;SERF1A::SMN1&quot;</td><td>&quot;5:70917048-5:7…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;out-of-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;TCGA-GI-A2C9&quot;</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;C6ORF47::BAG6_…</td><td>&quot;C6ORF47::BAG6&quot;</td><td>&quot;6:31660321-6:3…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;UTR&quot;</td><td>&quot;UTR&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;TCGA-GI-A2C9&quot;</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;SLC29A1::HSP90…</td><td>&quot;SLC29A1::HSP90…</td><td>&quot;6:44232428-6:4…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;in-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;TCGA-GI-A2C9&quot;</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;AC018362.1::XA…</td><td>&quot;AC018362.1::XA…</td><td>&quot;15:42532388-X:…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;exonic(no-know…</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;TCGA-GI-A2C9&quot;</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;DCN::AHNAK__12…</td><td>&quot;DCN::AHNAK&quot;</td><td>&quot;12:91180302-11…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;UTR&quot;</td><td>&quot;intronic&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;TCGA-GI-A2C9&quot;</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;AIMP2::STAG3L5…</td><td>&quot;AIMP2::STAG3L5…</td><td>&quot;7:6009378-7:10…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;TCGA-GI-A2C9&quot;</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;F11R::RBM4__1:…</td><td>&quot;F11R::RBM4&quot;</td><td>&quot;1:161019399-11…</td><td>&quot;-&quot;</td><td>&quot;+&quot;</td><td>&quot;intronic&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;TCGA-GI-A2C9&quot;</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;FBXO25::FAM157…</td><td>&quot;FBXO25::FAM157…</td><td>&quot;8:435707-9:138…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;CDS(truncated)…</td><td>&quot;exonic(no-know…</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;TCGA-GI-A2C9&quot;</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;MAGT1::STAC2__…</td><td>&quot;MAGT1::STAC2&quot;</td><td>&quot;X:77891455-17:…</td><td>&quot;-&quot;</td><td>&quot;-&quot;</td><td>&quot;intronic&quot;</td><td>&quot;UTR&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;TCGA-GI-A2C9&quot;</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;PTPRF::IGK@__1…</td><td>&quot;PTPRF::IGK@&quot;</td><td>&quot;1:43622869-2:8…</td><td>&quot;+&quot;</td><td>&quot;-&quot;</td><td>&quot;UTR&quot;</td><td>&quot;---&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;TCGA-GI-A2C9&quot;</td><td>&quot;FusionCatcher&quot;</td></tr><tr><td>&quot;EEF1AKMT3::TSF…</td><td>&quot;EEF1AKMT3::TSF…</td><td>&quot;12:57773128-12…</td><td>&quot;+&quot;</td><td>&quot;+&quot;</td><td>&quot;in-frame&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;.&quot;</td><td>&quot;TCGA-GI-A2C9&quot;</td><td>&quot;FusionCatcher&quot;</td></tr></tbody></table></div>




```python
summary_norms_df = joined_norms_df.group_by('breakpointPair').agg([
    pl.n_unique('sampleID').alias('unique_samples')
]).filter(
    pl.col('unique_samples') > 1
).sort('unique_samples', descending=False)
```


```python
# rename column breakpointPair to breakpointID
summary_norms_df = summary_norms_df.rename({"breakpointPair": "breakpointID"})
summary_norms_df
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (647, 2)</small><table border="1" class="dataframe"><thead><tr><th>breakpointID</th><th>unique_samples</th></tr><tr><td>cat</td><td>u32</td></tr></thead><tbody><tr><td>&quot;16:2097868-16:…</td><td>2</td></tr><tr><td>&quot;2:197503128-2:…</td><td>2</td></tr><tr><td>&quot;11:65295990-11…</td><td>2</td></tr><tr><td>&quot;16:2097728-16:…</td><td>2</td></tr><tr><td>&quot;5:158892090-1:…</td><td>2</td></tr><tr><td>&quot;5:81814204-5:8…</td><td>2</td></tr><tr><td>&quot;2:86741003-2:8…</td><td>2</td></tr><tr><td>&quot;3:68998011-6:4…</td><td>2</td></tr><tr><td>&quot;13:24017701-13…</td><td>2</td></tr><tr><td>&quot;16:2097728-16:…</td><td>2</td></tr><tr><td>&quot;3:186858134-5:…</td><td>2</td></tr><tr><td>&quot;12:26801-15:10…</td><td>2</td></tr><tr><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;17:15503098-17…</td><td>82</td></tr><tr><td>&quot;5:141376507-5:…</td><td>83</td></tr><tr><td>&quot;5:141357621-5:…</td><td>84</td></tr><tr><td>&quot;5:141366751-5:…</td><td>85</td></tr><tr><td>&quot;5:115616318-5:…</td><td>85</td></tr><tr><td>&quot;16:18858170-16…</td><td>90</td></tr><tr><td>&quot;5:141362556-5:…</td><td>90</td></tr><tr><td>&quot;1:85496166-1:8…</td><td>91</td></tr><tr><td>&quot;17:47592545-17…</td><td>101</td></tr><tr><td>&quot;11:117195915-1…</td><td>102</td></tr><tr><td>&quot;1:84563257-1:8…</td><td>102</td></tr><tr><td>&quot;7:99971746-7:9…</td><td>112</td></tr></tbody></table></div>




```python
summary_norms_df.write_csv("summary-norms-counts.tsv", separator="\t")
```


```python
## do the same to TNBC 

summary_tnbc_df = tnbc_only_df.group_by('breakpointPair').agg([
    pl.n_unique('sampleID').alias('unique_samples')
]).filter(
    pl.col('unique_samples') > 1
).sort('unique_samples', descending=False)
```


```python
summary_tnbc_df
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (623, 2)</small><table border="1" class="dataframe"><thead><tr><th>breakpointPair</th><th>unique_samples</th></tr><tr><td>cat</td><td>u32</td></tr></thead><tbody><tr><td>&quot;4:10068163-5:6…</td><td>2</td></tr><tr><td>&quot;6:73489274-19:…</td><td>2</td></tr><tr><td>&quot;14:102000717-6…</td><td>2</td></tr><tr><td>&quot;3:69100069-6:1…</td><td>2</td></tr><tr><td>&quot;5:151684611-19…</td><td>2</td></tr><tr><td>&quot;1:119428516-1:…</td><td>2</td></tr><tr><td>&quot;15:64821912-1:…</td><td>2</td></tr><tr><td>&quot;9:97290902-9:9…</td><td>2</td></tr><tr><td>&quot;10:71841710-3:…</td><td>2</td></tr><tr><td>&quot;1:1442852-1:14…</td><td>2</td></tr><tr><td>&quot;9:100442356-9:…</td><td>2</td></tr><tr><td>&quot;16:78726461-8:…</td><td>2</td></tr><tr><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;5:141362556-5:…</td><td>40</td></tr><tr><td>&quot;11:66617057-11…</td><td>48</td></tr><tr><td>&quot;1:85496166-1:8…</td><td>49</td></tr><tr><td>&quot;17:15545785-17…</td><td>50</td></tr><tr><td>&quot;11:66617057-11…</td><td>52</td></tr><tr><td>&quot;16:18858170-16…</td><td>55</td></tr><tr><td>&quot;17:15503098-17…</td><td>59</td></tr><tr><td>&quot;11:117195915-1…</td><td>62</td></tr><tr><td>&quot;5:115616318-5:…</td><td>73</td></tr><tr><td>&quot;7:99971746-7:9…</td><td>80</td></tr><tr><td>&quot;17:47592545-17…</td><td>98</td></tr><tr><td>&quot;1:84563257-1:8…</td><td>110</td></tr></tbody></table></div>




```python
# rename column breakpointPair to breakpointID
summary_tnbc_df = summary_tnbc_df.rename({"breakpointPair": "breakpointID"})
summary_tnbc_df
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (623, 2)</small><table border="1" class="dataframe"><thead><tr><th>breakpointID</th><th>unique_samples</th></tr><tr><td>cat</td><td>u32</td></tr></thead><tbody><tr><td>&quot;4:10068163-5:6…</td><td>2</td></tr><tr><td>&quot;6:73489274-19:…</td><td>2</td></tr><tr><td>&quot;14:102000717-6…</td><td>2</td></tr><tr><td>&quot;3:69100069-6:1…</td><td>2</td></tr><tr><td>&quot;5:151684611-19…</td><td>2</td></tr><tr><td>&quot;1:119428516-1:…</td><td>2</td></tr><tr><td>&quot;15:64821912-1:…</td><td>2</td></tr><tr><td>&quot;9:97290902-9:9…</td><td>2</td></tr><tr><td>&quot;10:71841710-3:…</td><td>2</td></tr><tr><td>&quot;1:1442852-1:14…</td><td>2</td></tr><tr><td>&quot;9:100442356-9:…</td><td>2</td></tr><tr><td>&quot;16:78726461-8:…</td><td>2</td></tr><tr><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;5:141362556-5:…</td><td>40</td></tr><tr><td>&quot;11:66617057-11…</td><td>48</td></tr><tr><td>&quot;1:85496166-1:8…</td><td>49</td></tr><tr><td>&quot;17:15545785-17…</td><td>50</td></tr><tr><td>&quot;11:66617057-11…</td><td>52</td></tr><tr><td>&quot;16:18858170-16…</td><td>55</td></tr><tr><td>&quot;17:15503098-17…</td><td>59</td></tr><tr><td>&quot;11:117195915-1…</td><td>62</td></tr><tr><td>&quot;5:115616318-5:…</td><td>73</td></tr><tr><td>&quot;7:99971746-7:9…</td><td>80</td></tr><tr><td>&quot;17:47592545-17…</td><td>98</td></tr><tr><td>&quot;1:84563257-1:8…</td><td>110</td></tr></tbody></table></div>




```python
# filter summary_df by breakpointID value in the summary_norms_df: if it exists in summary_norms_df, filter out the row in the summary_df

# Filter df1 to keep only rows where the breakpointID is not present in df2
norm_filtered_summary_tnbc_df = summary_tnbc_df.filter(~pl.col('breakpointID').is_in(summary_norms_df['breakpointID']))
```

    /var/folders/cj/48kxs2k94p7fmnmdbb42jxv00000gn/T/ipykernel_69903/430244449.py:4: CategoricalRemappingWarning: Local categoricals have different encodings, expensive re-encoding is done to perform this merge operation. Consider using a StringCache or an Enum type if the categories are known in advance
      norm_filtered_summary_tnbc_df = summary_tnbc_df.filter(~pl.col('breakpointID').is_in(summary_norms_df['breakpointID']))



```python
norm_filtered_summary_tnbc_df
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (474, 2)</small><table border="1" class="dataframe"><thead><tr><th>breakpointID</th><th>unique_samples</th></tr><tr><td>cat</td><td>u32</td></tr></thead><tbody><tr><td>&quot;4:10068163-5:6…</td><td>2</td></tr><tr><td>&quot;6:73489274-19:…</td><td>2</td></tr><tr><td>&quot;14:102000717-6…</td><td>2</td></tr><tr><td>&quot;3:69100069-6:1…</td><td>2</td></tr><tr><td>&quot;5:151684611-19…</td><td>2</td></tr><tr><td>&quot;15:64821912-1:…</td><td>2</td></tr><tr><td>&quot;10:71841710-3:…</td><td>2</td></tr><tr><td>&quot;1:1442852-1:14…</td><td>2</td></tr><tr><td>&quot;9:100442356-9:…</td><td>2</td></tr><tr><td>&quot;16:78726461-8:…</td><td>2</td></tr><tr><td>&quot;7:56021119-7:2…</td><td>2</td></tr><tr><td>&quot;10:42752102-22…</td><td>2</td></tr><tr><td>&hellip;</td><td>&hellip;</td></tr><tr><td>&quot;6:73510899-X:1…</td><td>14</td></tr><tr><td>&quot;2:89149998-8:5…</td><td>14</td></tr><tr><td>&quot;12:11918726-12…</td><td>14</td></tr><tr><td>&quot;19:1274440-19:…</td><td>14</td></tr><tr><td>&quot;15:88905757-15…</td><td>15</td></tr><tr><td>&quot;13:45378671-12…</td><td>15</td></tr><tr><td>&quot;6:31120881-6:2…</td><td>16</td></tr><tr><td>&quot;13:45383235-12…</td><td>18</td></tr><tr><td>&quot;22:36190233-13…</td><td>20</td></tr><tr><td>&quot;11:122538653-6…</td><td>20</td></tr><tr><td>&quot;6:79312368-6:7…</td><td>24</td></tr><tr><td>&quot;5:141352669-5:…</td><td>34</td></tr></tbody></table></div>




```python
norm_filtered_summary_tnbc_df.write_csv("norm-filt-summary-TNBC-only-counts.tsv", separator="\t")
```


```python
import matplotlib.pyplot as plt

# Count the number of unique breakpointIDs for each unique_samples value
sharedness_counts = (
    norm_filtered_summary_df
    .group_by('unique_samples')
    .agg(
        pl.col('breakpointID').n_unique()
        .alias('unique_bp_count')
    )
    .sort('unique_samples')
)

# Prepare the data for plotting
x = sharedness_counts['unique_samples'].to_list()
y = sharedness_counts['unique_bp_count'].to_list()

# Create the bar plot
plt.figure(figsize=(10, 8), dpi=300)
plt.bar(x, y)

# Set labels and title
plt.xlabel('Sharedness Score (Number of Patients A Unique FT Is Observed)')
plt.ylabel('Count of Unique FTs')
plt.title('Frequency of Unique Tumor-Specific FTs by Sharedness Score')

# Rotate x-axis labels for better readability
plt.xticks(rotation=90)

plt.show()
```


    
![png](README_files/README_36_0.png)
    



```python
# Count the number of unique breakpointIDs for each unique_samples value
sharedness_score = (
    norm_filtered_summary_tnbc_df
    .group_by('unique_samples')
    .agg(
        pl.col('breakpointID').n_unique()
        .alias('unique_ft_count')
    )
    .sort('unique_samples')
)
sharedness_score
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (19, 2)</small><table border="1" class="dataframe"><thead><tr><th>unique_samples</th><th>unique_ft_count</th></tr><tr><td>u32</td><td>u32</td></tr></thead><tbody><tr><td>2</td><td>289</td></tr><tr><td>3</td><td>76</td></tr><tr><td>4</td><td>38</td></tr><tr><td>5</td><td>20</td></tr><tr><td>6</td><td>15</td></tr><tr><td>7</td><td>9</td></tr><tr><td>8</td><td>5</td></tr><tr><td>9</td><td>2</td></tr><tr><td>10</td><td>3</td></tr><tr><td>11</td><td>2</td></tr><tr><td>12</td><td>1</td></tr><tr><td>13</td><td>1</td></tr><tr><td>14</td><td>5</td></tr><tr><td>15</td><td>2</td></tr><tr><td>16</td><td>1</td></tr><tr><td>18</td><td>1</td></tr><tr><td>20</td><td>2</td></tr><tr><td>24</td><td>1</td></tr><tr><td>34</td><td>1</td></tr></tbody></table></div>




```python
sharedness_score = sharedness_score.with_columns([
        pl.col('unique_samples').cast(pl.Utf8)
    ])
sharedness_score
```




<div><style>
.dataframe > thead > tr,
.dataframe > tbody > tr {
  text-align: right;
  white-space: pre-wrap;
}
</style>
<small>shape: (19, 2)</small><table border="1" class="dataframe"><thead><tr><th>unique_samples</th><th>unique_ft_count</th></tr><tr><td>str</td><td>u32</td></tr></thead><tbody><tr><td>&quot;2&quot;</td><td>289</td></tr><tr><td>&quot;3&quot;</td><td>76</td></tr><tr><td>&quot;4&quot;</td><td>38</td></tr><tr><td>&quot;5&quot;</td><td>20</td></tr><tr><td>&quot;6&quot;</td><td>15</td></tr><tr><td>&quot;7&quot;</td><td>9</td></tr><tr><td>&quot;8&quot;</td><td>5</td></tr><tr><td>&quot;9&quot;</td><td>2</td></tr><tr><td>&quot;10&quot;</td><td>3</td></tr><tr><td>&quot;11&quot;</td><td>2</td></tr><tr><td>&quot;12&quot;</td><td>1</td></tr><tr><td>&quot;13&quot;</td><td>1</td></tr><tr><td>&quot;14&quot;</td><td>5</td></tr><tr><td>&quot;15&quot;</td><td>2</td></tr><tr><td>&quot;16&quot;</td><td>1</td></tr><tr><td>&quot;18&quot;</td><td>1</td></tr><tr><td>&quot;20&quot;</td><td>2</td></tr><tr><td>&quot;24&quot;</td><td>1</td></tr><tr><td>&quot;34&quot;</td><td>1</td></tr></tbody></table></div>




```python
import matplotlib.pyplot as plt
import seaborn as sns
import polars as pl

# Create the bar plot
plt.figure(figsize=(10, 8), dpi=300)
sns.barplot(x=sharedness_score['unique_samples'],y=sharedness_score['unique_ft_count'], color='steelblue')

# Add value labels on top of the bars
for i, v in enumerate(sharedness_score['unique_ft_count']):
    plt.text(i, v, str(v), color='black', ha='center', fontweight='bold', fontsize=8)
    

# Set labels and title
plt.xlabel('Sharedness Score (Number of Patients A Unique FT Is Observed)')
plt.ylabel('Count of Unique FTs')
plt.title('Frequency of Unique Tumor-Specific FTs by Sharedness Score')

# Rotate x-axis labels for better readability
plt.xticks(rotation=90)

plt.show()
```


    
![png](README_files/README_39_0.png)
    

