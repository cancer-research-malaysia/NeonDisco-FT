## SAC Proposal (Project 1)

### Introduction
1. Cancer is driven by genetic abnormalities, specifically nonsynonymous DNA mutations, such as single nucleotide variants (*SNVs*) or small insertions and deletions (*indels*), and gene fusions due to chromosomal rearrangements or translocations (Vogelstein et al., 2013)

2. These genetic aberrations can translate into tumor-specific neoantigens – a class of altered, foreign protein products that are only expressed in one particular tumor cell and absent in normal, healthy cells

3. The resulting neoantigens can be highly immunogenic and may induce an antitumor adaptive immune response as the immune system will see them as 'nonself' (Gubin et al., 2015)

4. Neoantigens can originate from both mutation-dependent and mutation-independent molecular processes (Fig 1)[assets/katsikis.png], but   

4. Cancer neoantigens have canonically been studied in the context of SNVs in the protein-coding exome, but noncanonical neoantigens from genetic alterations that produce fusion genes and transcriptional aberrations that dysregulate normal alternative splicing events (ASE) remain largely unexplored (Capietto et al., 2022)

5. Tumor mutational burden demonstrates positive correlation with response to immunotherapy (Samstein et al., 2019), but many cancer types have low mutational burden, so 

6. Therefore, there is a


### Knowledge Gap

The idea of personalized neoantigen-based vaccine is exciting but logistically daunting and focusing research efforts on discovering shared immunogenic tumor-specific neoantigens that can be used to aid immunothrapy for a specific cancer occurring within a specific population of patinets may be more realistic and achievable.

A neoantigen database, specifically one that collates neoantigen data from sources beyond SN mutations and indels of coding regions of the genome, would expand the space of neoantigen repertoire to inform future neoantigen-based therapeutic development. We also attempt to contextualize genomic data obtained from local cohort of patients in Malaysia to create a more direct therapeutic impact and contribute back to the local community, both in general and in cancer clinical context.



### Main Objectives

A) Identification of shared, public neoantigen candidates found amAsian populations owith predicted high immunogenicity for translation into clinical therapy

B) Harnessing mutation-independent neoantigen sources to expand search of therapeutically relevant neoantigens for vaccine targets

### Proposed Methodology

* We have designed a neoantigen prediction framework (Fig 2) from which a Nextflow pipeline could be designed to identify potential neoantigens derived from alternative tumor-specific sources of genomic and transcriptomic alterations as described previously. 

The Nextflow pipeline is highly modular and would allow different prediction modules to be "plugged and played" as desired, based on a specific neoantigen-deriving aberrant source. 

#### ***Gene Fusions***

lorem ipsum

#### ***Aberrant ASE***

lorem ipsum

#### ***