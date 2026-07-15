# Genomic insights into the karyotypic radiation of a narrow endemic holocentric plant Carex helodes

This repository contains the computational workflows and scripts used in the analyses described in the manuscript: Genomic insights into the karyotypic radiation of a narrow endemic holocentric plant Carex helodes Link (doi: 10.1098/rstb.2015.0615)

**"Chromosomal rearrangements in *Carex helodes*: genome assembly, comparative genomics, and evolutionary analyses."**

The repository is organized according to the main methodological sections of the manuscript. Each directory contains the scripts used for a specific analytical component, together with a dedicated README describing the workflow and execution order.

---

# Repository structure

```
.
├── genome_assembly/
├── characterization_rearrangements/
├── gene_divergence/
├── gene_ontology/
├── graphs/
└── README.md
```

---

# Contents

## genome_assembly/

Pipeline for chromosome-scale genome assembly and comparative genome analyses.

Main analyses include:

- PacBio HiFi assembly with Hifiasm
- Duplicate purging with purge_dups
- Hi-C scaffolding using Arima pipeline and YaHS
- Manual curation in JuiceBox
- Quality assessment (BUSCO, Merqury, gfastats)
- Contaminant filtering
- Repeat library construction and genome masking
- Synteny analyses using GENESPACE
- Reciprocal Hi-C cross-mapping
- Structural variant validation using pbsv
- GC content analyses
- Final assembly quality control
- CR validation
- CR heterozygosity

This directory corresponds primarily to manuscript sections:

- **Genome assembly**
- **De novo repeat discovery and annotation, Gene annotation, synteny analysis, CRs validation and inversions heterozygosity and chromosomal rearrangement identification**

---

## characterization_rearrangements/

Scripts used to characterize genomic properties associated with conserved and rearranged genomic regions.

Analyses include:

- GC content
- Gene density
- Repeat content
- Repeat family enrichment
- Linear and generalized linear models
- Permutation tests using regioneR
- Breakpoint enrichment analyses
- Telomeric versus interstitial comparisons

This directory corresponds to manuscript section:

- **Characterization of conserved versus rearranged genomic regions**

---

## gene_divergence/

Pipeline for estimating sequence divergence and selective pressure associated with chromosomal rearrangements.

Main analyses include:

- Ortholog extraction
- CDS retrieval
- MAFFT alignments
- Codon alignments with PAL2NAL
- dN/dS estimation using PAML (codeml)
- Regional enrichment analyses
- Comparative analyses between rearranged and collinear regions

Corresponds to manuscript section:

- **Divergence and selection in rearrangements**

---

## gene_ontology/

Gene Ontology enrichment analyses of genes located within rearranged genomic regions.

Analyses include:

- GO enrichment with topGO
- EggNOG annotation processing
- Biological Process enrichment
- GO semantic similarity clustering
- Result summarization and visualization

Corresponds to manuscript section:

- **Gene ontology within rearranged regions**

---

## graphs/

Scripts used to generate the figures presented in the manuscript and supplementary material.

Includes:

- Genomic Variables autocorrelation analyses


---
