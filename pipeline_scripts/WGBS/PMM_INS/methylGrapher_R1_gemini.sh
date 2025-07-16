#!/bin/bash
## variables
samples=$1
# read each line from the $1 input sampleNames.txt file and create .sh script for each sample
while IFS= read -r sample_name; do
  echo "Creating bash script for sample: $sample_name"
  # Define the directory:
  datapath_PMM=/scratch/qwan/projects/PMM/rawData/allFastqs_wgbs
  mkdir -p /scratch/qwan/projects/PMM/DNAm/DNAm_hg38_p14_methylGrapher
  mkdir -p /scratch/qwan/projects/PMM/DNAm/DNAm_hg38_p14_methylGrapher/${sample_name}
  outdir=/scratch/qwan/projects/PMM/DNAm/DNAm_hg38_p14_methylGrapher/${sample_name}
  mkdir -p /scratch/qwan/projects/hg38_p14/methylGrapher_CHM13_index
  index_dir=/scratch/qwan/projects/hg38_p14/methylGrapher_CHM13_index
  # pkg path:
  methylGrapher=/home/qwan/miniconda3/envs/methylGrapher/bin/methylGrapher
  # prepaire data
  # get lambda sequence:
  lambda_phage_fa=/scratch/qwan/projects/hg38_p14/wgbs_CHM13_ref/lambda/lambda_phage.fasta
  # get pangenome T2T reference
  CHM13_Graph_dir=/scratch/qwan/projects/hg38_p14/wgbs_CHM13_ref/gfa
  CHM13_Graph_gz=/scratch/qwan/projects/hg38_p14/wgbs_CHM13_ref/gfa/hprc-v1.1-mc-chm13.gfa.gz
{
    echo -e "#!/bin/bash
#SBATCH --job-name=methylGrapher_R1_${sample_name}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 12 # Number of cores
#SBATCH -N 1 # Min - Max Nodes
#SBATCH -p compute
#SBATCH --mem=192G
#SBATCH --time=36:00:00
#SBATCH --output=methylGrapher_R1_${sample_name}_%j.log


# miniconda3 on gemini
source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/coh_PMM

# step 1 -- index (spikein: -lp ${lambda_phage_fa})
# pigz -d -p 8 -c ${CHM13_Graph_gz} > ${CHM13_Graph_dir}/hprc_v1_1_mc_chm13.gfa

# generate index 288G RAM needed
# methylGrapher PrepareGenome -gfa ${CHM13_Graph_dir}/hprc_v1_1_mc_chm13.gfa \\
#   -lp ${lambda_phage_fa} \\
#   -prefix ${index_dir}/CHM13 -t 12

# step 2 -- main: it automatically executes PrepareLibrary, Align and MethylCall in a sequence.
# pigz -d -p 8 -c ${datapath_PMM}/${sample_name}_R1.fastq.gz > ${datapath_PMM}/${sample_name}_R1.fastq
# pigz -d -p 8 -c ${datapath_PMM}/${sample_name}_R2.fastq.gz > ${datapath_PMM}/${sample_name}_R2.fastq

methylGrapher Main -t 12 \\
  -work_dir ${outdir} \\
  -index_prefix ${index_dir}/CHM13 \\
  -fq1 ${datapath_PMM}/${sample_name}_R1.fastq -fq2 ${datapath_PMM}/${sample_name}_R2.fastq

# step3 -- align
# methylGrapher Align -t 12 \\
#   -work_dir ${outdir} \\
#   -index_prefix ${index_dir}/CHM13

# step4 -- MethylCall
# methylGrapher MethylCall  -t 8 \\
#   -work_dir ${outdir} \\
#   -index_prefix ${index_dir}/CHM13

# step5 -- ConversionRate
methylGrapher ConversionRate -work_dir ${outdir} \\
  -index_prefix ${index_dir}/CHM13

conda deactivate"
  } > ${outdir}/${sample_name}_methylGrapher.sh
#  cd ${outdir}
#  sbatch ${sample_name}_methylGrapher.sh
done < ${samples}

echo "All sample script files created successfully."

