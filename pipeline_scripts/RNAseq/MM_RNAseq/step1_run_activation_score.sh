#!/bin/bash

# get activated L1 for each sample
# we use s5_activationScore_R4_gemini as input
# output dir: s6_activationScore_R3

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
  # datapath_input=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s5_activationScore_R4_gemini
  # mkdir -p /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/
  mkdir -p /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_redo
  mkdir -p /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_redo/${sample_name}
  outdir=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_redo/${sample_name}
  {
    echo -e "#!/bin/bash
#SBATCH --job-name=${sample_name}_redo_get_activeL1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=32G
#SBATCH --time=00:40:00
#SBATCH --output=${outdir}/${sample_name}_redo_get_activeL1_%j.log\n

module load R/RStudio_R-4.4.1

Rscript \
/home/qwan/githubRepo/coh_PMM/RNAseq_scripts/MMRF_analysis/redo_DE/step1_activation_score.R \
${sample_name}

"
  } > ${outdir}/${sample_name}_FilterActiveL1s.sh
#  cd ${outdir}
#  sbatch ${sample_name}_rnaPreprocess.sh
done < ${samples}

echo "All sample script files created successfully."




