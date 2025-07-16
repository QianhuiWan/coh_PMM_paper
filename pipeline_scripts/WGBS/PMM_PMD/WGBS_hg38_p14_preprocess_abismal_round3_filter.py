# aim to optimise abismal for more L1 repeats
# round3: optimise trim_galore parameters for increase abismal performance
# this round: add filter, filter out patches from hg38.p14, only keep chr1:22 + chrX + chrY + chrM

import sys
import os

config = sys.argv[1]

# software directory
bismark_methylation_extractor = "/home/qwan/miniconda3/envs/coh/bin/bismark_methylation_extractor"
Bismark_deduplicate = "/home/qwan/miniconda3/envs/coh/bin/deduplicate_bismark"
trim_galore = "/home/qwan/miniconda3/envs/coh/bin/trim_galore"
abismal = "/home/qwan/miniconda3/envs/coh/bin/abismal"
dnmtools = "/home/qwan/miniconda3/envs/coh/bin/dnmtools"
java = "/home/qwan/miniconda3/envs/coh/bin/java"
Samtools = "/home/qwan/miniconda3/envs/coh/bin/samtools"
wigToBigWig = "/home/qwan/miniconda3/envs/coh/bin/wigToBigWig"
bigWigToBedGraph = "/home/qwan/miniconda3/envs/coh/bin/bigWigToBedGraph"

# hg38 abismal index
hg38_p14_abismal = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/abismal_hg38_p14_filter" \
                   "/hg38_p14_abismalidx"

with open(config) as f:
    for line in f:
        Rpath, outputPath, sample = line.strip().split()
        if not os.path.exists(f"{outputPath}/{sample}"):
            os.mkdir(f"{outputPath}/{sample}")
        if not os.path.exists(f"{outputPath}/{sample}/01.QC"):
            os.mkdir(f"{outputPath}/{sample}/01.QC")
        if not os.path.exists(f"{outputPath}/{sample}/02.MAP"):
            os.mkdir(f"{outputPath}/{sample}/02.MAP")
        if not os.path.exists(f"{outputPath}/{sample}/03.DNMtoolsPosMethylation"):
            os.mkdir(f"{outputPath}/{sample}/03.DNMtoolsPosMethylation")
        with open(f"{outputPath}/{sample}/{sample}_wgbsPreprocess.sh", "w") as f:
            f.write(f"""#!/bin/bash
#SBATCH --job-name={sample}_abismal_wgbs_hg38_p14_filter 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 16
#SBATCH -N 1-4
#SBATCH -p all
#SBATCH --mem=150G
#SBATCH --time=48:00:00 # Time limit hrs:min:sec
#SBATCH --output={sample}_abismal_wgbs_hg38_p14_filter_%j.log

source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/coh

## 1. data clean #######################################################################################################
cd {outputPath}/{sample}/01.QC
ln -s {Rpath}/{sample}_R1.fastq.gz {sample}_R1.fq.gz
ln -s {Rpath}/{sample}_R2.fastq.gz {sample}_R2.fq.gz

{trim_galore} --paired -q 0 --length 0 \
--basename {sample} --cores 4 \
{sample}_R1.fq.gz {sample}_R2.fq.gz

# DNMtools recommended settings: trim_galore --paired -q 0 --length 0 human_esc_1.fastq human_esc_2.fastq

## 2. data mapping #####################################################################################################
cd {outputPath}/{sample}/02.MAP/
ln -s {outputPath}/{sample}/01.QC/{sample}*_1.fq.gz .
ln -s {outputPath}/{sample}/01.QC/{sample}*_2.fq.gz .

{abismal} -B -v -i  {hg38_p14_abismal} \
-o ./{sample}_pe.bam -t 16 \
-stats ./{sample}_abismal_stat.yaml \
{sample}_val_1.fq.gz {sample}_val_2.fq.gz

### format pe.bam ##########
{dnmtools} format -f abismal -t 8 -B {sample}_pe.bam {sample}_format_pe.bam

## 3. remove duplicates ################################################################################################
### Use dnmtools to remove duplicates##########
{Samtools} sort -@ 8 -o sorted_{sample}_pe.bam {sample}_format_pe.bam

{dnmtools} uniq -t 8 -S duplicate_removal_stats.txt sorted_{sample}_pe.bam nr_sorted_{sample}_pe.bam


## 4. DNAm extract #####################################################################################################
cd {outputPath}/{sample}/03.DNMtoolsPosMethylation/
ln -s {outputPath}/{sample}/02.MAP/nr_sorted_{sample}_pe.bam .

### Extract BS convertion rate##########
{dnmtools} bsrate -c /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_hg38_p14.fa \
-t 2 -o {sample}.bsrate nr_sorted_{sample}_pe.bam

### Extract all DNAm for C sites##########
{dnmtools} counts -c /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_hg38_p14.fa \
-t 2 -z -o {sample}.meth.gz nr_sorted_{sample}_pe.bam

### Extract all CpG sites##########
{dnmtools} sym -t 2 -z -o {sample}_CpG.meth.gz {sample}.meth.gz

### summary stats##########
{dnmtools} levels -o {sample}.levels {sample}.meth.gz

### CpG.meth.gz file to meth.bw file and meth.bedGraph##########
pigz -d -p 8 -c {sample}_CpG.meth.gz > {sample}_CpG.meth

awk '{{print $1,$2,$2+1,$5}}' {sample}_CpG.meth > {sample}_CpG_meth.bedGraph

{wigToBigWig} {sample}_CpG_meth.bedGraph \
/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_hg38_p14.chrom.sizes \
{sample}_CpG_meth.bw

### CpG.meth file to cov.bw file and cov.bedGraph##########
awk '{{print $1,$2,$2+1,$6}}' {sample}_CpG.meth > {sample}_CpG_cov.bedGraph

{wigToBigWig} {sample}_CpG_cov.bedGraph \
/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_hg38_p14.chrom.sizes \
{sample}_CpG_cov.bw

rm {sample}_CpG.meth

conda deactivate
""")
