#!/bin/sh

# In this round, we calculate gene level counts, use exons
# match with s1_checkStrandness_R4 and _getTxCounts_R5

# Check if the input file is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <input_file>"
  exit 1
fi

## variables
samples=$1

# Check if the input file exists
if [ ! -f "$samples" ]; then
  echo "Input file not found!"
  exit 1
fi

# read each line from the $1 input sampleNames.txt file and create .sh script for each sample
while IFS= read -r sample_name; do
  echo "Creating bash script for sample: $sample_name"
  ## data path
  datapath_MM=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/allFastqs_MMRF_rna
  s1_outdir=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s1_checkStrandness_R4/${sample_name}
  s2_outdir=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s2_getBW/${sample_name}
  s3_outdir=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s3_countTx_R5/${sample_name}
  mkdir -p /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2
  ## make directory for each sample
  mkdir -p /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/${sample_name}
  outdir=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/${sample_name}
  ## software
  featureCounts=/home/qwan/miniconda3/envs/coh/bin/featureCounts
  # annotations
  hg38_transcriptGTF=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/teAnno_round3/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf
  #hg38_mRNAexonsSAF=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14_saf/hg38_p14_mRNA.saf
  ## write .sh script for each sample
  {
    echo -e "#!/bin/bash
#SBATCH --job-name=${sample_name}_getGeneCounts_R2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 10 # Number of cores
#SBATCH -N 1-4 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=30G
#SBATCH --time=06:00:00
#SBATCH --output=${sample_name}_getGeneCounts_R2_%j.log\n

source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/coh

# cp data from s1 output
ln -s ${s1_outdir}/${sample_name}_sorted_nr_sorted.bam ${outdir}/${sample_name}_sorted_nr_sorted.bam
ln -s ${s1_outdir}/${sample_name}_sorted_nr_sorted.bam.bai ${outdir}/${sample_name}_sorted_nr_sorted.bam.bai

#### gene level counts #################################################################################################
##### -s 0 :un-stranded
${featureCounts} -B -p --countReadPairs -Q 30 -s 0 -T 8 -a ${hg38_transcriptGTF} \
-o ${outdir}/${sample_name}_gene_counts.txt \
${outdir}/${sample_name}_sorted_nr_sorted.bam

conda deactivate"
  } > ${outdir}/${sample_name}_rnaGetGeneCounts.sh
#  cd ${outdir}
#  sbatch ${sample_name}_rnaGetGeneCounts.sh
done < "${samples}"

echo "All sample script files created successfully."
