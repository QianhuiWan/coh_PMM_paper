# this round: use Bismark bowtie1
# this script was updated based on WGBS_hg38_p14_preprocess_round2.py (bismark bowtie2)
# move DNAm data and run on Gemini cpu

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
        with open(f"{outputPath}/{sample}/{sample}_wgbsPreprocess.sh", "w") as f:
            f.write(f"""#!/bin/bash
#SBATCH --job-name={sample}_bowtie1_wgbs_hg38_p14
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=qwan@coh.org      # Where to send mail
#SBATCH -n 16                         # Number of cores
#SBATCH -N 1                          # Min - Max Nodes
#SBATCH -p compute                    # gpu queue. Replace gpu-a100 with gpu-v100 to run on V100 cards
#SBATCH --mem=150G                    # Amount of memory in GB
#SBATCH --time=60:00:00               # Time limit hrs:min:sec
#SBATCH --output={sample}_bowtie1_wgbs_hg38_p14_%j.log   # Standard output and error log

# jobs on gemini
# module spider trimgalore/0.6.5

## 1. data clean #######################################################################################################
cd {outputPath}/{sample}/01.QC
ln -s {Rpath}/{sample}_R1.fastq.gz {sample}_R1.fq.gz
ln -s {Rpath}/{sample}_R2.fastq.gz {sample}_R2.fq.gz

source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/coh_PMM

# {trim_galore} --paired --clip_R1 8 --clip_R2 8 \
# --three_prime_clip_R1 9 --three_prime_clip_R2 9 \
# --basename {sample} --cores 4 \
# {sample}_R1.fq.gz {sample}_R2.fq.gz

conda deactivate

## --trim1: 1bp off every read from its 3' end, which may be needed for bowtie1, 
## --three_prime_clip_R1 1 --three_prime_clip_R2 1, this effectively achieves the same behavior as --trim1
## about 2h for trimming

## 2. data mapping #####################################################################################################
cd {outputPath}/{sample}/02.MAP/
ln -s {outputPath}/{sample}/01.QC/{sample}*_1.fq.gz .
ln -s {outputPath}/{sample}/01.QC/{sample}*_2.fq.gz .

source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/bismark_env

# {Bismark} --bowtie1 --path_to_bowtie /home/qwan/miniconda3/envs/coh_PMM/bin/ \
# --non_bs_mm --unmapped --phred33-quals --parallel 12 \
# --genome_folder {hg38_p14_bismark} \
# -1 {sample}_val_1.fq.gz -2 {sample}_val_2.fq.gz -o . 

conda deactivate

# about 12-24h for mapping

## 3. remove duplicates ################################################################################################
source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/coh_PMM

### Use picard to remove duplicates, keep the same with previous analysis from Joo
# {java} -Djava.io.tmpdir=/scratch/qwan/temp \
# -jar /home/qwan/miniconda3/envs/coh_PMM/share/picard-2.27.5-0/picard.jar \
# SortSam -INPUT {sample}_*_pe.bam -OUTPUT sorted_{sample}_pe.bam -SORT_ORDER coordinate

# {java} -Djava.io.tmpdir=/scratch/qwan/temp \
# -jar /home/qwan/miniconda3/envs/coh_PMM/share/picard-2.27.5-0/picard.jar \
# MarkDuplicates -INPUT sorted_{sample}_pe.bam -OUTPUT nr_sorted_{sample}_pe.bam \
# -REMOVE_DUPLICATES true \
# -METRICS_FILE ./picard-stats.txt

## 4. DNA methylation extract #############################################################################################
cd {outputPath}/{sample}/03.BismarkPosMethylation/
ln -s {outputPath}/{sample}/02.MAP/nr_sorted_{sample}_pe.bam .

# ### sort deduplicated bam e.g. 18_S5-name-sorted.bam
# {Samtools} sort -n -@ 14 -O bam -o {sample}_name_sorted.bam nr_sorted_{sample}_pe.bam

### Extract DNAm CpG sites
# {bismark_methylation_extractor} -p --no_overlap --comprehensive --merge_non_CpG \
# --gzip --bedGraph --counts --buffer_size 20G --cytosine_report --zero_based --parallel 16 \
# --genome_folder {hg38_p14_bismark} \
# {sample}_name_sorted.bam \
# --samtools_path {Samtools} 

### bedGraph to bigWig format
pigz -d -p 8 -c {sample}_name_sorted.bedGraph.gz > {sample}_name_sorted.bedGraph

{wigToBigWig} {sample}_name_sorted.bedGraph \
/scratch/qwan/projects/hg38_p14/filtered_hg38_p14.chrom.sizes \
{sample}_name_sorted.bw

# rm {sample}_name_sorted.bedGraph
conda deactivate
""")

