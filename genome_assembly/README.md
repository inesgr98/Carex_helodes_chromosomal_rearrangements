
# Carex helodes Genome Assembly & Comparative Genomics Pipeline

This repository contains the complete suite of SLURM (`.sbatch`) and R scripts used for the genome assembly, *de novo* repeat annotation, and validation of large chromosomal rearrangements (CRs) in two high-quality, chromosome-scale assemblies of *Carex helodes*: **Aznalcóllar (AZN)** and **Monchique (MON)**.

## 📋 Table of Contents
1. [Workflow Overview](#-workflow-overview)
2. [Data Summary](#-data-summary)
3. [Repository Structure & Execution Order](#-repository-structure--execution-order)
4. [Detailed Component Breakdown](#-detailed-component-breakdown)
   - [Phase A: Assembly & Quality Control](#phase-a-assembly--quality-control)
   - [Phase B: Hi-C Scaffolding & Cross-Mapping](#phase-b-hi-c-scaffolding--cross-mapping)
   - [Phase C: De Novo Repeat Discovery & Masking](#phase-c-de-novo-repeat-discovery--masking)
   - [Phase D: Synteny & Structural Variant Validation](#phase-d-synteny--structural-variant-validation)

---

## 🧬 Workflow Overview

The computational pipeline is structured into four integrated phases designed to process PacBio HiFi long-reads and Omni-C/Hi-C data, generate highly accurate assemblies, isolate a custom pan-genome transposable element (TE) library, and validate structural evolutions (inversions, translocations, and fusions).


```

[ PacBio HiFi Reads ] ──> A1_hifiasm ──> A2_purge_dups ──> Chromosome Scaffolding
│
[ Transposable Elements ] <── C11_TE_library <── C0_cleanup <─────┤
│                                                      │
[ Host Gene Filtering ] ──> C14 to C16 ──> RepeatMasker <─────────┤
▼
[ Structural Validation ] <── genespace.R <── pbsv_pipeline <── B1/B2/B3 Hi-C Cross-Mapping

```

---

## 📊 Data Summary

The scripts in this repository were tuned to process the following genomic inputs for *Carex helodes*:
* **AZN (Aznalcóllar):** 60 Gb PacBio HiFi data (~42.4× coverage) | 145 Gb Hi-C data (392M reads).
* **MON (Monchique):** 72 Gb PacBio HiFi data (~49.2× coverage) | 204 Gb Hi-C data (569M reads).
* **Outputs:** * `AZN`: 346.74 Mb assembled into **35 chromosome-scale pseudomolecules** (Scaffold N50: 10.12 Mb).
  * `MON`: 348.51 Mb assembled into **36 chromosome-scale pseudomolecules** (Scaffold N50: 9.69 Mb).

---

## 🗂️ Repository Structure & Execution Order

Scripts are prefix-coded to indicate their logical sequence in the workflow:

| Step | Script Name | Language / Engine | Core Task / Target |
| :--- | :--- | :--- | :--- |
| **A1** | `A1_hifiasm.sbatch` | SLURM / `hifiasm` | Primary *de novo* long-read assembly |
| **A2** | `A2_histogram.sbatch` | SLURM / `GenomeScope` | K-mer based heterozygosity/profile analysis |
| **A2** | `A2_purge_dups_np.sbatch` | SLURM / `purge_dups` | Purging haplotypic duplications |
| **A2** | `A2_purge_dups_np_4C.sbatch` | SLURM / `purge_dups` | Optimized duplicate purging for 4C/structural contexts |
| **B1** | `B1_hiC_trimming_np.sbatch` | SLURM / `fastp`/`trimmomatic` | Pre-processing and cleaning raw Hi-C reads |
| **B1** | `B1_hiC_index.sbatch` | SLURM / `bwa index` | Indexing assemblies for Hi-C read mapping |
| **B1** | `B1_hiC_mapping.sbatch` | SLURM / `bwa mem` | Native alignment of Hi-C pairs to self-reference |
| **B2** | `B2_juicer_bf_mc.sbatch` | SLURM / `Juicer Tools` | Generating Hi-C matrices *Before Manual Curation* |
| **B3** | `B3_juicer_af_mc.sbatch` | SLURM / `Juicer Tools` | Compiling refined Hi-C contacts *After Manual Curation* |
| **B4** | `B4_decontaminate.sbatch` | SLURM / `BlobTools`/`Kraken` | Screening and removal of exogenous/microbial contigs |
| **C0** | `C0_remove_small_contigs.sbatch` | SLURM / `awk` / Python | Truncating short unanchored contigs for clean TE runs |
| **C11**| `C11_TE_library.sbatch` | SLURM / `EarlGrey` | Constructing individual *de novo* repeat libraries |
| **C11**| `C11_TE_library_wR.sbatch` | SLURM / `EarlGrey` | Integrating TREP and PlantSat reference libraries |
| **C12**| `C12_TE_lib_concadenated_wR.sbatch` | SLURM / `EarlGrey` | Compiling unified cross-assembly concatenated libraries |
| **C13**| `C13_cluster_library_80.sbatch` | SLURM / `CD-HIT` | Clustering consensus sequences based on the 80/80/80 rule |
| **C14**| `C14_prepare_host_genes_for_filtering.sbatch` | SLURM / `BLAST+` | Database construction of functional plant protein-coding genes |
| **C15**| `C15_filter_host_genes_r80.sbatch` | SLURM / `BLASTN` | Eliminating host genes from the transposable element library |
| **C16**| `C16_finish_host_genes.sbatch` | SLURM / Custom Script | Generating the final common *Carex helodes* RE library |
| **-** | `unmask_genome.sbatch` | SLURM / `RepeatMasker` | Production-level soft/hard masking of final genomes |
| **-** | `hiC_mapping_crossed_refMON.sbatch` | SLURM / `bwa mem` | Reciprocal cross-mapping: AZN Hi-C reads onto MON assembly |
| **-** | `hiC_mapping_crossed_refAZN.sbatch` | SLURM / `bwa mem` | Reciprocal cross-mapping: MON Hi-C reads onto AZN assembly |
| **-** | `index_assemblies.sbatch` | SLURM / `samtools` | Creating spatial index structures (`.fai`) for mapping support |
| **-** | `pbsv_refAZN` | SLURM / PacBio `pbsv` | Calling structural variants and structural signature discovery |
| **-** | `bed_extract_reords.sbatch` | SLURM / `bedtools` | Extracting coverage intervals ± 5000 bp around predicted breakpoints |
| **-** | `genespace.R` | R Environment / `GENESPACE` | Comparative synteny mapping and ancestral block tracking |
| **-** | `syntheny.sbatch` | SLURM / `MCScanX` / `minimap2` | Micro/macro-synteny visualization extraction pipelines |
| **-** | `syntheny_reord_4C.sbatch` | SLURM / Custom script | Reordering syntenic data blocks for clean 4C validation plots |
| **-** | `gc_content.sbatch` | SLURM / `bedtools` | Positional calculation of GC skew across chromosomal pseudomolecules |
| **-** | `final_QC.sbatch` | SLURM / `BUSCO` | Final validation of gene space completeness (Poales vs Embryophyta) |

---

## 🛠️ Detailed Component Breakdown

### Phase A: Assembly & Quality Control
* **`A1_hifiasm.sbatch`** computes primary configurations using PacBio HiFi inputs. 
* **`A2_histogram.sbatch`** verifies the low heterozygosity profiles calculated by GenomeScope 1.0 (**0.23% for AZN** and **0.26% for MON**).
* **`A2_purge_dups_np*.sbatch`** eliminates alternative haplotypic overlaps, establishing precise genome architectures of **346.74 Mb** (AZN) and **348.51 Mb** (MON).

### Phase B: Hi-C Scaffolding & Cross-Mapping
* **`B1_hiC_*.sbatch`** maps sequence pairings natively to place scaffolds into continuous pseudo-chromosomes.
* **`B2_juicer_bf_mc.sbatch`** and **`B3_juicer_af_mc.sbatch`** handle structural tracking pre- and post- manual configuration via Juicebox.
* **`hiC_mapping_crossed_refMON.sbatch`** and **`hiC_mapping_crossed_refAZN.sbatch`** provide the data for reciprocal cross-mapping. By mapping AZN reads onto the MON reference and vice versa, structural discordances can be directly visualized to identify true evolutionary rearrangements instead of assembly errors.

### Phase C: De Novo Repeat Discovery & Masking
To minimize the fraction of unclassified elements, a customized plant repetitive element workflow was deployed:
1. **Library Gathering (`C11`):** `EarlGrey` parses elements using a custom base built with the **TREP database** (Transposable Elements) and **PlantSat** (Satellites).
2. **Cross-Annotation (`C12`):** Repeats identified from the Poales library in AZN are dynamically run against MON to establish a standardized, shared pan-species *Carex helodes* repeat database.
3. **Host Gene Filtering (`C14` - `C16`):** Elements overlapping actual functional plant coding frames are isolated and stripped away (`C15_filter_host_genes_r80.sbatch`) to prevent false-positive masking.
4. **Production Run (`unmask_genome.sbatch`):** Coordinates final `RepeatMasker v4.1.5` execution, identifying that **50.20% (AZN)** and **50.38% (MON)** of the genomes consist of repetitive DNA.

### Phase D: Synteny & Structural Variant Validation
* **`genespace.R`** evaluates micro-synteny block structural relationships, identifying **3 explicit inversions** (`Iinv_chr6`, `Iinv_chr8`, `Iinv_chr32`), **4 translocations** (`Ttrans_chr5`, `Ttrans_chr9`, `Ttrans_chr16`, `Ttrans_chr24`), and **1 fusion event** where ancestral segments `MON-chr17` and `MON-chr19` fused into the largest modern AZN chromosome.
* **`pbsv_refAZN`** triggers structural call parameters (`discover` and `call`) using long split-reads.
* **`bed_extract_reords.sbatch`** inspects structural breakends (BNDs) within a window of ± 5000 bp to provide confirmation statistics (e.g., confirming the heterozygous inversion in `chr32` for MON where 76 out of 130 split-reads support the variant).

---



```

```
