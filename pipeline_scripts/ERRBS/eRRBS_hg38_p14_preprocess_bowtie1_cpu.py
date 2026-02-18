# this round: use Bismark bowtie1 for eRRBS data
# this script was adapted from WGBS_hg38_p14_preprocess_bowtie1_cpu.py for eRRBS
# CHANGES:
# 1. trim_galore uses --rrbs
# 2. NO deduplication (skip Picard MarkDuplicates)
# 3. Direct mapping output to methylation extractor

import sys
import os

config = sys.argv[1]

# software directory
Bismark = "/home/qwan/miniconda3/envs/bismark_env/bin/bismark"
bismark_methylation_extractor = "/home/qwan/miniconda3/envs/bismark_env/bin/bismark_methylation_extractor"
# Bismark_deduplicate = "/home/qwan/miniconda3/envs/coh_PMM/bin/deduplicate_bismark"
trim_galore = "/home/qwan/miniconda3/envs/coh_PMM/bin/trim_galore"
dnmtools = "/home/qwan/miniconda3/envs/coh_PMM/bin/dnmtools"
java = "/home/qwan/miniconda3/envs/coh_PMM/bin/java"
Samtools = "/home/qwan/miniconda3/envs/coh_PMM/bin/samtools"
wigToBigWig = "/home/qwan/miniconda3/envs/coh_PMM/bin/wigToBigWig"
bigWigToBedGraph = "/home/qwan/miniconda3/envs/coh_PMM/bin/bigWigToBedGraph"

# hg38 bismark index
hg38_p14_bismark = "/scratch/qwan/projects/hg38_p14/bismark_0_19_1_bowtie1_hg38_p14_filter"

with open(config) as f:
    for line in f:
        Rpath, outputPath, sample = line.strip().split()
        if not os.path.exists(f"{outputPath}/{sample}"):
            os.mkdir(f"{outputPath}/{sample}")
        if not os.path.exists(f"{outputPath}/{sample}/01.QC"):
            os.mkdir(f"{outputPath}/{sample}/01.QC")
        if not os.path.exists(f"{outputPath}/{sample}/02.MAP"):
            os.mkdir(f"{outputPath}/{sample}/02.MAP")
        if not os.path.exists(f"{outputPath}/{sample}/03.BismarkPosMethylation"):
            os.mkdir(f"{outputPath}/{sample}/03.BismarkPosMethylation")
        with open(f"{outputPath}/{sample}/{sample}_eRRBSPreprocess.sh", "w") as f:
            f.write(f"""#!/bin/bash
#SBATCH --job-name={sample}_eRRBS_bowtie1_hg38_p14
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=qwan@coh.org      # Where to send mail
#SBATCH -n 16                         # Number of cores
#SBATCH -N 1                          # Min - Max Nodes
#SBATCH -p compute                    # gpu queue
#SBATCH --mem=150G                    # Amount of memory in GB
#SBATCH --time=60:00:00               # Time limit hrs:min:sec
#SBATCH --output={outputPath}/{sample}/{sample}_eRRBS_bowtie1_hg38_p14_%j.log   # Standard output and error log

# jobs on gemini
# module spider trimgalore/0.6.5

## 1. data clean #######################################################################################################
cd {outputPath}/{sample}/01.QC
ln -s {Rpath}/{sample}_R1.fastq.gz {sample}_R1.fq.gz
ln -s {Rpath}/{sample}_R2.fastq.gz {sample}_R2.fq.gz

source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/coh_PMM

# eRRBS trimming: Use --rrbs
# This automatically handles MspI site adaptation and 2bp removal if needed.
# {trim_galore} --rrbs --paired --basename {sample} --cores 4 \\
# {sample}_R1.fq.gz {sample}_R2.fq.gz

conda deactivate

## 2. data mapping #####################################################################################################
cd {outputPath}/{sample}/02.MAP/
ln -s {outputPath}/{sample}/01.QC/{sample}*_1.fq.gz .
ln -s {outputPath}/{sample}/01.QC/{sample}*_2.fq.gz .

source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/bismark_env

## note path adjustment if needed, but keeping original structure
# {Bismark} --bowtie1 --path_to_bowtie /home/qwan/miniconda3/envs/coh_PMM/bin/ \\
# --non_bs_mm --unmapped --phred33-quals --parallel 12 \\
# --genome_folder {hg38_p14_bismark} \\
# -1 {sample}_val_1.fq.gz -2 {sample}_val_2.fq.gz -o . 

conda deactivate

# about 12-24h for mapping

## 3. remove duplicates (SKIPPED FOR eRRBS) ############################################################################
# RRBS data should NOT be deduplicated by position as many reads share start sites (MspI).
# We proceed directly to extraction using the mapped BAM.

## 4. DNA methylation extract #############################################################################################
cd {outputPath}/{sample}/03.BismarkPosMethylation/
# Link directly to the mapped BAM (not deduplicated)
# The Bismark output usually follows pattern: {sample}_val_1_bismark_bt2_pe.bam or similar.
# For bowtie1 it is usually: {sample}_val_1_bismark_pe.bam (check exact output name!)
# Assuming logical naming convention from mapping step output.

rm {sample}_val_1_bismark_pe.bam
ln -s {outputPath}/{sample}/02.MAP/{sample}_val_1_bismark_pe.bam .

# ### sort bam by name for methylation extractor
# {Samtools} sort -n -@ 14 -O bam -o {sample}_name_sorted.bam {sample}_val_1_bismark_pe.bam

### Extract DNAm CpG sites
# {bismark_methylation_extractor} -p --no_overlap --comprehensive --merge_non_CpG \\
# --gzip --bedGraph --counts --buffer_size 20G --cytosine_report --zero_based --parallel 16 \\
# --genome_folder {hg38_p14_bismark} \\
# {sample}_name_sorted.bam \\
# --samtools_path {Samtools} 


### bedGraph to bigWig format
# Make CpG percent bedGraph (0-100)
# Uses the coverage file from Bismark: {sample}_name_sorted.bismark.cov.gz
# zcat {sample}_name_sorted.bismark.cov.gz | \\
#   awk 'BEGIN{{OFS="\\t"}} {{s=$2-1; if(s<0)s=0; print $1,s,$3,$4}}' > {sample}_CpG_percent.bedGraph

# sort -k1,1 -k2,2n {sample}_CpG_percent.bedGraph > {sample}_CpG_percent_sorted.bedGraph

{wigToBigWig} {sample}_CpG_percent_sorted.bedGraph \\
  /scratch/qwan/projects/hg38_p14/filtered_hg38_p14.chrom.sizes \\
  {sample}_CpG_percent.bw

# rm {sample}_CpG_percent.bedGraph
conda deactivate
""")
