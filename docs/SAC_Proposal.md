## SAC Proposal (Project 1)

### Introduction
1. Cancer is driven by genetic abnormalities, specifically nonsynonymous DNA mutations, such as single nucleotide variants (*SNVs*) or small insertions and deletions (*indels*), and gene fusions due to chromosomal rearrangements or translocations (Vogelstein et al., 2013)

2. These genetic aberrations can translate into tumor-specific neoantigens – a class of altered, foreign protein products that are only expressed in one particular tumor cell and absent in normal, healthy cells

3. The resulting neoantigens can be highly immunogenic and may induce an antitumor adaptive immune response as the immune system will see them as 'nonself' (Gubin et al., 2015)

4. Neoantigens can be generalized into mutation-dependent and mutation-independent neoantigens (**Figure 1**) that could be produced through various molecular mechanisms at different levels of regulation (genetic, transcriptional, translational, post-translational).

![Figure 1: Categories of neoantigens (adapted from Katsikis et al., 2023)](assets/katsikis.png)
<em>Figure 1: Categories of neoantigens (adapted from Katsikis et al., 2023)</em>

5. Cancer neoantigens have canonically been studied in the context of SNVs in the protein-coding exome, but noncanonical neoantigens from genetic alterations that produce fusion genes and transcriptional aberrations that dysregulate normal alternative splicing events (ASE) remain largely unexplored (Capietto et al., 2022)

6. These alternative sources of neoantigens can theoretically produce tumor-specific neopeptides (**Figure 2**) that can bind to MHC molecules to mount an adaptive immune response in patients, especially in tumor types with low tumor mutational burden, a proxy metric that tends to be biased towards SNV-derived neoantigens.

![Figure 2: Theoretical sources of neoantigens (adapted from Bräunlein & Krackhardt, 2017)](assets/braunlein&krackhardt.jpg)
<em>Figure 2: Theoretical sources of neoantigens (adapted from Bräunlein & Krackhardt, 2017)</em>


### Significance 

* Personalized neoantigen-based vaccine development targets a set of both unique, *private* neoantigens alongside recurrent, *public* neoantigens. Such vaccine targets are tailored to individual patients and thus impose prohibitive logistical and financial contraints to widespread clinical accessibility (Pearlman et al., 2021)

* Focusing research efforts on discovering shared (public) tumor-specific neoantigens that can be targeted via an *of-the-shelf* immunotherapy may be more realistic and may speed up clinical translation

* Compiling a neoantigen database that covers neoantigens from sources beyond single base mutations and indels of coding regions of the genome would expand the space of neoantigen repertoire to inform future neoantigen-based therapeutic development

* Contextualizing neoantigen discovery with a specific Asian population would also address genomic data inequity that is currently inherent in genomic precision medicine field that is predominantly Western or European-centric (Tawfik et al., 2023), hopefully leading to a more direct therapeutic impact for local Asian community in Malaysia

**Project Objective** 
> **Identification of highly immunogenic and shared, public neoantigen candidates in Asian populations by harnessing neoantigen-producing sources beyond SNVs to expand neoantigen prediction space for a faster and more comprehensive clinical translation**

### Proposed Methodology

* We have designed a Nextflow pipeline following a generalized neoantigen prediction framework (**Figure 3**) to identify potential neoantigens derived from alternative tumor-specific sources of genomic and transcriptomic alterations as described previously. 

* The Nextflow pipeline is highly modular and would allow different prediction modules to be "plugged and played" as desired, based on a specific neoantigen-deriving aberrant source. Two aberrant neoantigen-producing mechanisms have been chosen as our primary focus in establishing an Asian cancer neoantigen database.

	1. ***Gene fusion neopeptides***: They represent the aberrant output of gene fusions caused by chromosomal genetic rearrangements via translocations, deletions, or inversions

	2. ***Alternatively spliced neopeptides***: We will focus on *aberrant* splicing events that produce neopeptides that are absent in normal cells

#### ***Gene Fusions***

As a preliminary step, we have optimized a gene fusion transcript calling using two published bioinformatics tools – Arriba and FusionCatcher.

#### ***Aberrant ASE***

lorem ipsum

