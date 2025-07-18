# coh_PMM_paper

_Analysis codes for L1 repeat activity in Plasmablastic Multiple Myeloma (PMM)_

## Project Overview

This repository contains code and analysis pipelines for the **Plasmablastic Multiple Myeloma (PMM)** project, focusing on the characterization of **L1 repeat family** activity using RNA-seq and whole-genome bisulfite sequencing (WGBS) data.

PMM is a high-risk subtype of multiple myeloma characterized by immature plasma cells (plasmablasts) in the bone marrow, showing aggressive behavior and reduced survival. These immature plasma cells display distinct morphological features, including a large nucleus, diffuse chromatin, and a prominent nucleolus.

---

## Data Summary

- **Samples:**
  - **RNA-seq:** PMM3, PMM4, PMM6, PMM7, PMM11, PMM14, PMM15, PMM16; BMPC1, BMPC2, BMPC3
  - **WGBS:** PMM1, PMM2, PMM3, PMM4, PMM6, PMM7, PMM9, PMM11, PMM12, PMM13, PMM14, PMM15, PMM16, PMM17, PMM18; B1_rest, B2_rest, B3_rest

- **Genome build:** hg38 (for final analyses)

---

## Pipeline Summary

### RNA-seq Analysis

1. **Preprocessing**
   - Alignment: `STAR`
   - Indexing & sorting: `samtools`
   - Duplicate removal: `picard`

2. **Gene counts**
   - Quantification: `STAR`

3. **Strand separation**
   - Separate Watson and Crick strand reads: `samtools`

4. **Visualization**
   - Convert bam to BigWig: `bamCoverage`

5. **Transposable Element (TE) transcripts**
   - Called separately on Watson and Crick strands, with activation scores calculated using an in-house pipeline

---

### WGBS Analysis

1. **Preprocessing and CpG DNA methylation calling**
   - Quality control & trimming: `trimGalore`
   - Alignment: `bismark` (bowtie1)
   - Deduplication: `picard`
   - Methylation extraction: `bismark_methylation_extractor`
   - Coverage filtering: `methylKit`

2. **Region calling**
   - Quality control & trimming: `trimGalore`
   - Alignment: `abismal`
   - Deduplication & methylation extraction: `DNMTools`
   - Identification of PMD (partially methylated domains): `DNMTools`

3. **INS analysis**
   - Alignment: `methylGrapher`
   - Deduplication & methylation extraction: `methylGrapher`
   - INS supporting reads detection: in-house pipeline

---

## Repository Structure

```
/pipeline_scripts    # Analysis pipelines
/figure_scripts      # Code for generating paper figures
```

---

## How to Run

```bash
git clone https://github.com/QianhuiWan/coh_PMM_paper
# Change paths in the scripts to match your local data directory
# Install all required packages specified in the scripts
# Then run the pipeline scripts using `bash`:
bash /pipeline_scripts/*.sh

# Figures are generated using R version 4.4.1
# Run the figure scripts with `Rscript` after correcting the data paths
Rscript /figure_scripts/*.R
```

---

## License

[MIT License]

---

## Contact

For questions, please contact **Qianhui Wan** at [qwan@coh.org](mailto:qwan@coh.org).


