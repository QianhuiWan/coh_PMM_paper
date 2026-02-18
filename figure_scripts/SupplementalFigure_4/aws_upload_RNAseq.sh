#!/bin/bash
#SBATCH --job-name=aws_upload_R1_RNAseq_addBMPC
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=32G
#SBATCH --time=24:00:00 # Time limit hrs:min:sec
#SBATCH --output=aws_upload_R1_RNAseq_addBMPC_%j.log

source /home/qwan/.bashrc
module load Mamba/24.3.0-0
mamba activate /home/qwan/miniconda3/envs/coh

#   "CD19-GSM2171212" = "BMPC4",
#   "CD19-GSM2171213" = "BMPC5",
#   "CD19-GSM2171214" = "BMPC6",
#   "CD19-GSM2171215" = "BMPC7"


ln -s /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/CD19-GSM2171212_1/bam_split/*_rpkm.bw /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/rpkm_bw/BMPC4_1_rpkm.bw
ln -s /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/CD19-GSM2171212_2/bam_split/*_rpkm.bw /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/rpkm_bw/BMPC4_2_rpkm.bw
ln -s /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/CD19-GSM2171212_3/bam_split/*_rpkm.bw /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/rpkm_bw/BMPC4_3_rpkm.bw

ln -s /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/CD19-GSM2171213_1/bam_split/*_rpkm.bw /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/rpkm_bw/BMPC5_1_rpkm.bw
ln -s /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/CD19-GSM2171213_2/bam_split/*_rpkm.bw /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/rpkm_bw/BMPC5_2_rpkm.bw
ln -s /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/CD19-GSM2171213_3/bam_split/*_rpkm.bw /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/rpkm_bw/BMPC5_3_rpkm.bw

outdir="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/rpkm_bw"
basedir="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1"

declare -A map=(
  [GSM2171212]=BMPC4
  [GSM2171213]=BMPC5
  [GSM2171214]=BMPC6
  [GSM2171215]=BMPC7
)

for gsm in GSM2171212 GSM2171213 GSM2171214 GSM2171215; do
  sample="${map[$gsm]}"
  if [ -z "$sample" ]; then
    echo "ERROR: no sample mapping for $gsm" >&2
    continue
  fi

  for rep in 1 2 3; do
    src_glob="$basedir/CD19-${gsm}_${rep}/bam_split/"*_rpkm.bw
    files=( $src_glob )

    if [ "${#files[@]}" -ne 1 ]; then
      echo "ERROR: expected 1 match, got ${#files[@]}: $src_glob" >&2
      continue
    fi

    ln -sfn "${files[0]}" "$outdir/${sample}_${rep}_rpkm.bw"
    echo "linked -> $outdir/${sample}_${rep}_rpkm.bw"
  done
done

cp /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/rpkm_bw/*.bw /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_BMPC_hubDir/RNA/
cp /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/rpkm_bw/*.bw /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_BMPC_hubDir/RNA/


cp /home/qwan/githubRepo/coh_PMM/upload_data/bw_tracks/pmm_addBMPC_RNA/hub.txt /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_BMPC_hubDir
cp /home/qwan/githubRepo/coh_PMM/upload_data/bw_tracks/pmm_addBMPC_RNA/genomes.txt /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_BMPC_hubDir
cp /home/qwan/githubRepo/coh_PMM/upload_data/bw_tracks/pmm_addBMPC_RNA/trackDb.txt /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_BMPC_hubDir/RNA/


aws s3 cp /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_BMPC_hubDir/hub.txt s3://qwan/PMM_addBMPC/hub.txt
aws s3 cp /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_BMPC_hubDir/genomes.txt s3://qwan/PMM_addBMPC/genomes.txt
aws s3 cp /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/pmm_BMPC_hubDir/RNA/ s3://qwan/PMM_addBMPC/RNA --recursive

# https://genome.ucsc.edu/cgi-bin/hgHubConnect?hubCheckUrl=https%3A%2F%2Fqwan.s3.amazonaws.com%2FPMM%2Fhub.txt&hgsid=3570381813_OmI6IaYUairDunBy8sqduyUSOlqE#hubDeveloper
hubCheck https://qwan.s3.amazonaws.com/PMM_addBMPC/hub.txt


# rsync -avPL /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/*/PMM*_ReadsPerGene.out.tab /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/GEO_qwan_2025/processed_data_GEO/geo_submission_RNAseq

# md5sum /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/GEO_qwan_2025/processed_data_GEO/geo_submission_RNAseq/* > /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/GEO_qwan_2025/processed_data_GEO/md5sum_RNAseq.txt

# check files:
aws s3 ls s3://qwan --recursive

mamba deactivate
