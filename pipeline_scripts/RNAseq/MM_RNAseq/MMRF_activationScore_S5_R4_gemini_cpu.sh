#!/bin/sh

# In this round, we get up&downstream counts for TE and random regions with featureCounts instead of pysam
# s5_activationScore_R4 match with s1_checkStrandness_R4 geneCount_S4_R2 getBW_S2_R2 and _getTxCounts_R5

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
  ### Define the default BAM path
  s1_outdir=/scratch/qwan/projects/PMM/RNA_MMRF/s1_checkStrandness_R4/${sample_name}
  ### Check if the directory contains any *_sorted_nr_sorted.bam files
  if ls "${s1_outdir}"/*_sorted_nr_sorted.bam > /dev/null 2>&1; then
    : # Do nothing, keep the default bamPath
  else
    # Switch to the alternative path
    s1_outdir=/scratch/qwan/projects/PMM/RNA_MMRF/s1_checkStrandness_R4_apollo/${sample_name}
  fi
  # Print the final bamPath
  echo "Final bamPath for ${sample_name}: ${s1_outdir}"

  mkdir -p /scratch/qwan/projects/PMM/RNA_MMRF/s5_activationScore_R4
  ## make directory for each sample
  mkdir -p /scratch/qwan/projects/PMM/RNA_MMRF/s5_activationScore_R4/${sample_name}
  outdir=/scratch/qwan/projects/PMM/RNA_MMRF/s5_activationScore_R4/${sample_name}
  ## software
  featureCounts=/home/qwan/miniconda3/envs/coh_PMM/bin/featureCounts
  samtools=/home/qwan/miniconda3/envs/coh_PMM/bin/samtools
  # annotations
  #hg38_transcriptGTF=/scratch/qwan/projects/PMM/RNA_MMRF/annoFiles/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf
  filtered_rmsk_hg38_start100bp_SAF=/scratch/qwan/projects/PMM/RNA_MMRF/annoFiles/filtered_rmsk_hg38_start100bp.saf
  filtered_rmsk_hg38_end100bp_SAF=/scratch/qwan/projects/PMM/RNA_MMRF/annoFiles/filtered_rmsk_hg38_end100bp.saf
  #filtered_random_hg38_start100bp_SAF=/scratch/qwan/projects/PMM/RNA_MMRF/annoFiles/filtered_random_hg38_start100bp.saf
  #filtered_random_hg38_end100bp_SAF=/scratch/qwan/projects/PMM/RNA_MMRF/annoFiles/filtered_random_hg38_end100bp.saf
  filtered_random_hg38_start100bp_SAF=/scratch/qwan/projects/PMM/RNA_MMRF/annoFiles/filtered_random5683690_hg38_start100bp.saf
  filtered_random_hg38_end100bp_SAF=/scratch/qwan/projects/PMM/RNA_MMRF/annoFiles/filtered_random5683690_hg38_end100bp.saf

  ## write .sh script for each sample
  {
    echo -e "#!/bin/bash
#SBATCH --job-name=${sample_name}_region100bpCounts_S5R4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 10 # Number of cores
#SBATCH -N 1 # Min - Max Nodes
#SBATCH -p compute
#SBATCH --mem=40G
#SBATCH --time=7:00:00
#SBATCH --output=${outdir}/${sample_name}_region100bpCounts_S5R4_%j.log

source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/coh_PMM

# cp data from s1 output
ln -s ${s1_outdir}/${sample_name}_sorted_nr_sorted.bam ${outdir}/${sample_name}_sorted_nr_sorted.bam
ln -s ${s1_outdir}/${sample_name}_sorted_nr_sorted.bam.bai ${outdir}/${sample_name}_sorted_nr_sorted.bam.bai

#### region level counts ###############################################################################################
##### get rmsk_hg38_start100bp_counts
##### -s 0 :un-stranded
${featureCounts} -F SAF -B -p --countReadPairs -Q 30 -s 0 -T 8 -a ${filtered_rmsk_hg38_start100bp_SAF} \
-o ${outdir}/${sample_name}_rmsk_hg38_start100bp_counts.txt \
${outdir}/${sample_name}_sorted_nr_sorted.bam

##### get rmsk_hg38_end100bp_counts
${featureCounts} -F SAF -B -p --countReadPairs -Q 30 -s 0 -T 8 -a ${filtered_rmsk_hg38_end100bp_SAF} \
-o ${outdir}/${sample_name}_rmsk_hg38_end100bp_counts.txt \
${outdir}/${sample_name}_sorted_nr_sorted.bam

##### get random_hg38_start100bp_counts
${featureCounts} -F SAF -B -p --countReadPairs -Q 30 -s 0 -T 8 -a ${filtered_random_hg38_start100bp_SAF} \
-o ${outdir}/${sample_name}_random_hg38_start100bp_counts.txt \
${outdir}/${sample_name}_sorted_nr_sorted.bam

##### get random_hg38_end100bp_counts
${featureCounts} -F SAF -B -p --countReadPairs -Q 30 -s 0 -T 8 -a ${filtered_random_hg38_end100bp_SAF} \
-o ${outdir}/${sample_name}_random_hg38_end100bp_counts.txt \
${outdir}/${sample_name}_sorted_nr_sorted.bam

conda deactivate"
  } > ${outdir}/${sample_name}_rna100bpCounts.sh
#  cd ${outdir}
#  sbatch ${sample_name}_rna100bpCounts.sh
done < "${samples}"

echo "All sample script files created successfully."
