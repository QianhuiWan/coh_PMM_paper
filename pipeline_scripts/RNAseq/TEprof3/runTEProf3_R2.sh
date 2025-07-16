#!/bin/bash
#SBATCH --job-name=teprof3_R2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 36 # Number of cores
#SBATCH -N 1-4 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=teprof3_R2_%j.log

# R2: use --quanreadlength 50

source /home/qwan/.bashrc
module load Mamba/24.3.0-0
mamba activate /home/qwan/miniconda3/envs/coh

# Define the directory:
teprof3Dir="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/TEProf3_TEdetect_R2"
# mkdir -p ${teprof3Dir}
cd ${teprof3Dir}

# ln -s /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/*/*_sorted_nr_sorted.bam .
# ln -s /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/*/*_sorted_nr_sorted.bam.bai .
# ln -s /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/*/*_SJ.out.tab .

# get sample_manifest file
# rm sample_manifest.txt
# find . -maxdepth 1 -name "*.bam" | while read file ; do xbase=$(basename $file); sample_name=${xbase/_sorted_nr_sorted.bam/} ; echo -e "${sample_name}\tshort\t${xbase}" >> sample_manifest.txt; done ;
# find . -maxdepth 1 -name "*.tab" | while read file ; do xbase=$(basename $file); sample_name=${xbase/_SJ.out.tab/} ; echo -e "${sample_name}\tSJ\t${xbase}" >> sample_manifest.txt; done ;
# find . -maxdepth 1 -name "*.gtf" | while read file ; do xbase=$(basename $file); sample_name=${xbase/.gtf/} ; echo -e "${sample_name}\tgtf\t${xbase}" >> sample_manifest.txt; done

# reset all folders if error log
teprof3 -f sample_manifest.txt --reset

# run teprof3:
# Assemble step (-as) \ # run Process assemble step (-ps) \ # run Filter transcripts (-fs) \\ 
# run Mega assembly (--tacothread) \ # run Quantification (-qs)
teprof3 -f sample_manifest.txt -s 11 -am 1 -as 11 --assemblethread 4 --assemblestrand 1 --assemblejunctionread 2 \
  -ps 11 -pt 0.5 -ptn 100 \
  -fs 11 -fm 1 --filterintronretention 3 \
  --filtermonoexontpm 1 --filterdownstreammate 2 --filterratio 0.5 \
  --tacothread 11 \
  -qs 11 --quansamplenumbercon 100 --quanmode 1 --quanreadlength 50

mamba deactivate





