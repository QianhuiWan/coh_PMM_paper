#!/bin/bash
#SBATCH --job-name=norm_R3_methylKit
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=192G
#SBATCH --time=24:00:00 # Time limit hrs:min:sec
#SBATCH --output=norm_R3_methylKit_%j.log

# Using SLURM_JOB_ID in the script
echo "This job's ID is $SLURM_JOB_ID"

module load R/RStudio_R-4.4.1

# round 1: only save rds
# Rscript /home/qwan/githubRepo/coh_PMM/WGBS_scripts/chr_plots_deeptools/fig2D_deeptools_prepare_normBW_R2.R

# round2: save output CpG.txt files to methylRaw_output, and use multi-cores
# Rscript /home/qwan/githubRepo/coh_PMM/WGBS_scripts/chr_plots_deeptools/fig2D_deeptools_prepare_normBW_R2_testMC.R
# 
# mv /home/qwan/githubRepo/coh_PMM/WGBS_scripts/chr_plots_deeptools/Rplots.pdf \
# /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/methylKit_output/plots/"Rplots_wgbs_${SLURM_JOB_ID}.pdf"

# round3: save output CpG.txt files to methylRaw_output, and use multi-cores
# Rscript /home/qwan/githubRepo/coh_PMM/WGBS_scripts/chr_plots_deeptools/fig2D_deeptools_prepare_normBW_R3_MC_cov1.R
# 
# mv /home/qwan/githubRepo/coh_PMM/WGBS_scripts/chr_plots_deeptools/Rplots.pdf \
# /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/methylKit_output/plots/"Rplots_wgbs_${SLURM_JOB_ID}.pdf"

# round4: save output CpG.txt files to methylRaw_output, and use multi-cores
Rscript /home/qwan/githubRepo/coh_PMM/WGBS_scripts/redo_DMR/step1_prepare_normBW_MC.R

mv /home/qwan/githubRepo/coh_PMM/WGBS_scripts/redo_DMR/Rplots.pdf \
/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/methylKit_output/plots/"Rplots_wgbs_${SLURM_JOB_ID}.pdf"

