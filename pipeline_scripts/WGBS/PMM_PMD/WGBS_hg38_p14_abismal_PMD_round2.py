
# current PMD pipeline: round2 analysis for PMDs etc.
# Use filtered hg38.p14, patches filtered out
# also use default entropy and entropy_F in this round:
# - default entropy is for methylated entropy sum(p(M)*log2(p(M)),
# - entropy_F -F also considered un-methylated entropy: sum(p(U)*log2(p(U)) when prob(M)<0.5

import sys
import os

config = sys.argv[1]
#currentpath = os.getcwd()

# software directory
abismal = "/home/qwan/miniconda3/envs/coh/bin/abismal"
dnmtools = "/home/qwan/miniconda3/envs/coh/bin/dnmtools"
Samtools = "/home/qwan/miniconda3/envs/coh/bin/samtools"

with open(config) as f:
    for line in f:
        Rpath, outputPath, sample = line.strip().split()
        if not os.path.exists(f"{outputPath}/{sample}/04.HMR"):
            os.mkdir(f"{outputPath}/{sample}/04.HMR")
        if not os.path.exists(f"{outputPath}/{sample}/05.HyperMR"):
            os.mkdir(f"{outputPath}/{sample}/05.HyperMR")
        if not os.path.exists(f"{outputPath}/{sample}/06.entropy_F"):
            os.mkdir(f"{outputPath}/{sample}/06.entropy_F")
        if not os.path.exists(f"{outputPath}/{sample}/07.PMD"):
            os.mkdir(f"{outputPath}/{sample}/07.PMD")
        if not os.path.exists(f"{outputPath}/{sample}/08.Hydroxymethylation"):
            os.mkdir(f"{outputPath}/{sample}/08.Hydroxymethylation")
        with open(f"{outputPath}/{sample}/{sample}_wgbsManalyses.sh", "w") as f:
            f.write(f"""#!/bin/bash
#SBATCH --job-name={sample}_DNAmDomainAnalysis 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org 
#SBATCH -n 8
#SBATCH -N 1-4
#SBATCH -p all
#SBATCH --mem=50G
#SBATCH --time=8:00:00 # Time limit hrs:min:sec 
#SBATCH --output={sample}_DNAmDomainAnalysis_%j.log

source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/coh

## 5. HMR hypo-methylated regions and PMR partially methylated regions##################################################
### we can start with the sample_CpG.meth.gz for getting HMRs
# cd {outputPath}/{sample}/04.HMR/
# ln -s {outputPath}/{sample}/03.DNMtoolsPosMethylation/{sample}_CpG.meth.gz .

### get hypo-methylated regions
# {dnmtools} hmr -v -p hmr_params.txt -o {sample}_CpG.hmr {sample}_CpG.meth.gz

### get partially-methylated regions
# {dnmtools} hmr -partial -v -p pmr_params.txt -o {sample}_CpG.pmr {sample}_CpG.meth.gz

## 6. HyperMR ##########################################################################################################
# cd {outputPath}/{sample}/05.HyperMR/
# ln -s {outputPath}/{sample}/03.DNMtoolsPosMethylation/{sample}_CpG.meth.gz .

### get hyper-methylated regions
# pigz -d -p 8 -c {sample}_CpG.meth.gz > {sample}_CpG.meth
# awk '{{$5=1-$5; print $0}}' {sample}_CpG.meth > {sample}_CpG_inverted.meth

# {dnmtools} hmr -v -p hmr_params.txt -o {sample}_CpG_inverted.hmr {sample}_CpG_inverted.meth

# rm {sample}_CpG.meth {sample}_CpG_inverted.meth

### get hyper-methylated domains
# {dnmtools} hypermr -v -p hypermr_params.txt -o {sample}_CpG.hypermr {sample}_CpG.meth.gz

## 7. entropy ##########################################################################################################
# cd {outputPath}/{sample}/06.entropy_F/
# ln -s {outputPath}/{sample}/03.DNMtoolsPosMethylation/{sample}_CpG.meth.gz .
# ln -s {outputPath}/{sample}/02.MAP/nr_sorted_{sample}_pe.bam .
# 
# {dnmtools} states -c /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_hg38_p14.fa \
# -o {sample}.epiread nr_sorted_{sample}_pe.bam
# 
# LC_ALL=C sort -k1,1 -k2,2g {sample}.epiread -o {sample}_sorted.epiread
# 
# {dnmtools} entropy -w 5 -F -v -o {sample}_entropy.meth \
# /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_chroms \
# {sample}_sorted.epiread

## 8. PMD ##############################################################################################################
# cd {outputPath}/{sample}/07.PMD/
# ln -s {outputPath}/{sample}/03.DNMtoolsPosMethylation/{sample}_CpG.meth.gz .

## PMD for single PMM sample and B cells
# {dnmtools} pmd -v -p pmd_params.txt -i 1000 -o {sample}_pmd.bed {sample}_CpG.meth.gz

## PMD for all samples (this part only need to run once)
# mkdir -p {outputPath}/pmd_all_samples
# cd {outputPath}/pmd_all_samples
# {dnmtools} pmd -v -p {outputPath}/pmd_all_samples/pmd_params.txt \
# -i 1000 -o {outputPath}/pmd_all_samples/all_samples_pmd.bed \
# {outputPath}/*/03.DNMtoolsPosMethylation/*_CpG.meth.gz

## PMD for all PMM samples (this part also only need to run once)
mkdir -p {outputPath}/pmd_allPMM_samples
cd {outputPath}/pmd_allPMM_samples
{dnmtools} pmd -v -p {outputPath}/pmd_allPMM_samples/pmd_params.txt \
-i 1000 -o {outputPath}/pmd_allPMM_samples/allPMM_samples_pmd.bed \
{outputPath}/PMM*/03.DNMtoolsPosMethylation/*_CpG.meth.gz

conda deactivate
""")
