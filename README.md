# coh_PMM_paper

## Project Overview

This repository contains code and analysis pipelines for the **Plasmablastic Multiple Myeloma (PMM)** project, focusing on the characterization of **L1 repeat family** activity using RNA-seq and whole-genome bisulfite sequencing (WGBS) data.

PMM is a high-risk subtype of multiple myeloma characterized by immature plasma cells (plasmablasts) in the bone marrow, showing aggressive behavior and reduced survival. These cells display distinct morphological features, including a large nucleus, diffuse chromatin, and a prominent nucleolus.

---

## Data Summary

- **Samples:**
  - **RNA-seq:** PMM3, PMM4, PMM6, PMM7, PMM11, PMM14, PMM15, PMM16; BMPC1, BMPC2, BMPC3
  - **WGBS:** PMM1, PMM2, PMM3, PMM4, PMM6, PMM7, PMM9(*), PMM11, PMM12, PMM13, PMM14, PMM15, PMM16, PMM17, PMM18; B1_rest, B2_rest, B3_rest

- **Genome build:** hg38 for final analyses

---

## Pipeline Summary

### RNA-seq Analysis

1. **Preprocessing**
   - Alignment: `STAR`
   - Indexing & sorting: `samtools`
   - Duplicate removal: `picard`

2. **Gene counts**
   - Identification: `STAR`

3. **Strand separation**
   - Separate + and – strand reads: `samtools`

4. **Visualization**
   - Convert bedGraph to BigWig: UCSC tools

5. **Transposable Element (TE) transcripts**
   - Called separately on Watson and Crick strands, with activation scores calculated using an in-house pipeline

---

### WGBS Analysis

1. **Preprocessing and CpG DNAm calling**
   - QC & trimming: `trimGalore`
   - Alignment: `bismark` (bowtie1)
   - Deduplication: `picard`
   - Methylation extraction: `bismark_methylation_extractor`

2. **Region calling**
   - QC & trimming: `trimGalore`
   - Alignment: `abismal`
   - Deduplication & methylation extraction: `DNMTools`
   - Identify PMD (partially methylated domains) using `DNMTools`

---

## Repository Structure

```
/pipeline_scripts    # Analysis pipelines
/figure_scripts      # code for figures in the paper
```

---

## How to Run

git clone https://github.com/QianhuiWan/coh_PMM_paper
change path in the script to the path of your data 
bash /pipeline_scripts/*.sh

---

## License

[MIT]

---

## Contact

For questions, please contact Qianhui Wan at qwan@coh.org




