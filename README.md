# Breakpoint Enrichment Analysis

This repository contains R scripts for analyzing and comparing genomic features (gene density, GC content, repeat elements) between conserved (collinear) and rearranged genomic regions, alongside breakpoint enrichment analysis.

## 📂 Repository Contents

* **`enrichment_analysis.R`** / **`enrichment_analysis_chr.R`**
  * Performs permutation-based enrichment analysis of genomic features across the whole genome and specific chromosomes using `regioneR`.
* **`enrichment_analysis_chr_area_depend.R`**
  * Runs section-specific (subtelomeric vs. interstitial) breakpoint enrichment analyses.
* **`transform_robj_table_chr.R`** / **`transform_robj_table_chr_area_depend.R`**
  * Data tidying scripts to transform R data objects into structured tables for statistical modeling and downstream visualization.

---

## 🔬 Methodology Overview

###  Breakpoint Enrichment Analysis
* Statistical enrichment of genomic features at breakpoint regions is evaluated using **1000 permutations** via the `regioneR` package.
* **Evaluation Statistics:** `numOverlaps` for categorical features; `meanInRegions` for GC content.
* **Spatial Controls:** Telomeric regions (terminal 10kb) are masked. Analysis is contextualized by partitioning chromosomes into subtelomeric (terminal 10%) and interstitial sections based on GC-content stabilization.

---

## 🛠️ Key R Dependencies

* **Genomics:** `GenomicRanges`, `regioneR`
* **Statistics & Modeling:** `lme4`, `glmmTMB`, `MASS`
* **Visualization & Utilities:** `ggcorrplot`, `stats`
