#!/bin/bash
#SBATCH --job-name=plot_L1_bindKZFPs_R2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --output=plot_L1_bindKZFPs_R2_%j.log

source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/coh

module load R/RStudio_R-4.4.1

outdir=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs

########## 293T chip-exo #################
## R1: use MACS score which is -log10(p-value) from chip-exo bed files to 
## plot KZFPs binding at full length L1

# Rscript \
# /home/qwan/githubRepo/coh_PMM/RNAseq_scripts/L1_ORF1_5UTR_predict/use_encode_chipSeq/L1_all_bin_plot_L1bindKZFPs_R1.R


# R2: only look at L1PA2
Rscript \
/home/qwan/githubRepo/coh_PMM/RNAseq_scripts/L1_ORF1_5UTR_predict/use_encode_chipSeq/L1_all_bin_plot_L1bindKZFPs_R2.R


conda deactivate


