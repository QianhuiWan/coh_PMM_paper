#!/bin/bash
#SBATCH --job-name=aws_upload_R1_WGBS
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=32G
#SBATCH --time=24:00:00 # Time limit hrs:min:sec
#SBATCH --output=aws_upload_R1_WGBS_%j.log

source /home/qwan/.bashrc
module load Mamba/24.3.0-0
mamba activate /home/qwan/miniconda3/envs/coh

#   "CD19-GSM2171212" = "BMPC4",
#   "CD19-GSM2171213" = "BMPC5",
#   "CD19-GSM2171214" = "BMPC6",
#   "CD19-GSM2171215" = "BMPC7"


cp /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/bowtie1_bw_perc/PMM*.bw /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_DNAm_hubDir/DNAm/
cp /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/bowtie1_bw_perc/B*_rest*.bw /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_DNAm_hubDir/DNAm/

cp /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/eRRBS_hg38_p14_bowtie1_filter/*/03.BismarkPosMethylation/*_CpG_percent.bw /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_DNAm_hubDir/DNAm/


cp /home/qwan/githubRepo/coh_PMM/upload_data/bw_tracks/pmm_DNAm/hub.txt /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_DNAm_hubDir
cp /home/qwan/githubRepo/coh_PMM/upload_data/bw_tracks/pmm_DNAm/genomes.txt /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_DNAm_hubDir
cp /home/qwan/githubRepo/coh_PMM/upload_data/bw_tracks/pmm_DNAm/trackDb.txt /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_DNAm_hubDir/DNAm/


aws s3 cp /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_DNAm_hubDir/hub.txt s3://qwan/PMM_DNAm/hub.txt
aws s3 cp /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_DNAm_hubDir/genomes.txt s3://qwan/PMM_DNAm/genomes.txt
aws s3 cp /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_DNAm_hubDir/DNAm/ s3://qwan/PMM_DNAm/DNAm --recursive

# https://genome.ucsc.edu/cgi-bin/hgHubConnect?hubCheckUrl=https%3A%2F%2Fqwan.s3.amazonaws.com%2FPMM%2Fhub.txt&hgsid=3570381813_OmI6IaYUairDunBy8sqduyUSOlqE#hubDeveloper
hubCheck https://qwan.s3.amazonaws.com/PMM_addBMPC/hub.txt
hubCheck https://qwan.s3.amazonaws.com/PMM_DNAm/hub.txt


# rsync -avPL /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/*/PMM*_ReadsPerGene.out.tab /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/GEO_qwan_2025/processed_data_GEO/geo_submission_RNAseq

# md5sum /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/GEO_qwan_2025/processed_data_GEO/geo_submission_RNAseq/* > /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/GEO_qwan_2025/processed_data_GEO/md5sum_RNAseq.txt

# check files:
aws s3 ls s3://qwan --recursive

mamba deactivate
