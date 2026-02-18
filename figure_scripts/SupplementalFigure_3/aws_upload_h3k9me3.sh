#!/bin/bash
#SBATCH --job-name=aws_upload_R1_RNAseq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=32G
#SBATCH --time=24:00:00 # Time limit hrs:min:sec
#SBATCH --output=aws_upload_R1_RNAseq_%j.log

source /home/qwan/.bashrc
module load Mamba/24.3.0-0
mamba activate /home/qwan/miniconda3/envs/coh

rsync -avP /Users/qwan/githubRepo/coh_PMM/data/chromHMM/GSE151556_RAW/GSM4585183_MM1S_H3K9me3_treat_pileup.bw qwan@apollo-acc.coh.org:/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/histone_hubDir/histone/
rsync -avP /Users/qwan/githubRepo/coh_PMM/data/chromHMM/GSE151556_RAW/GSM4585177_KMS11_H3K9me3_treat_pileup.bw qwan@apollo-acc.coh.org:/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/histone_hubDir/histone/


cp /home/qwan/githubRepo/coh_PMM/upload_data/bw_tracks/mm_cellline_h3k9me3/hub.txt /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/histone_hubDir
cp /home/qwan/githubRepo/coh_PMM/upload_data/bw_tracks/mm_cellline_h3k9me3/genomes.txt /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/histone_hubDir
cp /home/qwan/githubRepo/coh_PMM/upload_data/bw_tracks/mm_cellline_h3k9me3/trackDb.txt /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/histone_hubDir/histone/


aws s3 cp /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/histone_hubDir/hub.txt s3://qwan/PMM_histone/hub.txt
aws s3 cp /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/histone_hubDir/genomes.txt s3://qwan/PMM_histone/genomes.txt
aws s3 cp /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/histone_hubDir/histone/ s3://qwan/PMM_histone/histone --recursive

# https://genome.ucsc.edu/cgi-bin/hgHubConnect?hubCheckUrl=https%3A%2F%2Fqwan.s3.amazonaws.com%2FPMM%2Fhub.txt&hgsid=3570381813_OmI6IaYUairDunBy8sqduyUSOlqE#hubDeveloper
hubCheck https://qwan.s3.amazonaws.com/PMM_histone/hub.txt


# rsync -avPL /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/*/PMM*_ReadsPerGene.out.tab /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/GEO_qwan_2025/processed_data_GEO/geo_submission_RNAseq

# md5sum /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/GEO_qwan_2025/processed_data_GEO/geo_submission_RNAseq/* > /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/GEO_qwan_2025/processed_data_GEO/md5sum_RNAseq.txt

# check files:
aws s3 ls s3://qwan --recursive

mamba deactivate
