#!/bin/bash
#SBATCH --job-name=deltaDNAm_againstHEPG2_Wave_hg38
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=deltaDNAm_againstHEPG2_Wave_hg38_%j.log

source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/coh

module load R/RStudio_R-4.4.1

Rscript \
/home/qwan/githubRepo/coh_PMM/WGBS_scripts/wgbs_stats/deltaDNAm_waveScorePlot.R

conda deactivate

