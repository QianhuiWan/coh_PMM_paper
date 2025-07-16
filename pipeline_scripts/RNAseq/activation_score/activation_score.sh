#!/bin/bash

activation_score_script() {
  local sample_name=$1
  local crick_bam=$2
  local watson_bam=$3
  local out_directory=$4

  # check input
  if [[ -z "$sample_name" || -z "$crick_bam" || -z "$watson_bam" || -z "$out_directory" ]]; then
    echo "Usage: generate_fc_script <sample_name> <crick_bam> <watson_bam> <out_directory>"
    return 1
  fi

  if [[ ! -f "$crick_bam" ]]; then
    echo "Crick BAM file not found: $crick_bam"
    return 1
  fi

  if [[ ! -f "$watson_bam" ]]; then
    echo "Watson BAM file not found: $watson_bam"
    return 1
  fi

  # path settings
  local outdir="${out_directory}"
  mkdir -p "$outdir"

  local featureCounts=/home/qwan/miniconda3/envs/coh/bin/featureCounts
  local SAF_DIR=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo
  local SAF_start=${SAF_DIR}/filtered_rmsk_hg38_start100bp.saf
  local SAF_end=${SAF_DIR}/filtered_rmsk_hg38_end100bp.saf
  local SAF_random_start=${SAF_DIR}/filtered_random5683690_hg38_start100bp.saf
  local SAF_random_end=${SAF_DIR}/filtered_random5683690_hg38_end100bp.saf

  local script_path="${outdir}/${sample_name}_activationScore.sh"

  # heredoc
  cat <<EOF > "${script_path}"
#!/bin/bash
#SBATCH --job-name=${sample_name}_activationScore
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=${outdir}/${sample_name}_activationScore_%j.log

source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/coh

# test1: rm --countReadPairs
# test2: add -O

# # === featureCounts: rmsk_hg38_start100bp ===
${featureCounts} -F SAF -B -p -Q 30 -s 0 -T 8 -a ${SAF_start} \\
  -o ${outdir}/${sample_name}_rmsk_hg38_start100bp_counts.txt \\
  ${watson_bam} ${crick_bam}

# # === featureCounts: rmsk_hg38_end100bp ===
${featureCounts} -F SAF -B -p -Q 30 -s 0 -T 8 -a ${SAF_end} \\
  -o ${outdir}/${sample_name}_rmsk_hg38_end100bp_counts.txt \\
  ${watson_bam} ${crick_bam}

# # === featureCounts: random_hg38_start100bp ===
${featureCounts} -F SAF -B -p -Q 30 -s 0 -T 8 -a ${SAF_random_start} \\
  -o ${outdir}/${sample_name}_random_hg38_start100bp_counts.txt \\
  ${watson_bam} ${crick_bam}

# # === featureCounts: random_hg38_end100bp ===
${featureCounts} -F SAF -B -p -Q 30 -s 0 -T 8 -a ${SAF_random_end} \\
  -o ${outdir}/${sample_name}_random_hg38_end100bp_counts.txt \\
  ${watson_bam} ${crick_bam}

# === R analysis for L1 scores ===
module load R/RStudio_R-4.4.1
mkdir -p ${outdir}/filtered_L1
Rscript \
/home/qwan/githubRepo/coh_PMM/RNAseq_scripts/benchmark/activationScore_functions/activation_score_S2_fil_L1.R \
${sample_name} ${outdir} ${outdir}/filtered_L1

echo "activation socre finished"

conda deactivate
EOF

  echo "Script created and submitted: $script_path"
  sbatch ${script_path}
  return 0
}


