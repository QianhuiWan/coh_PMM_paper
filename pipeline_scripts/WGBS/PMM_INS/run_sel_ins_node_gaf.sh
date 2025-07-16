#!/bin/bash

dataPath=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_methylGrapher_gemini
mkdir -p /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_methylGrapher_gemini/ins_supporting_reads
outDir=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_methylGrapher_gemini/ins_supporting_reads


for gaf in ${dataPath}/*/alignment.gaf; do
  sample=$(basename "$(dirname "${gaf}")")  # Get folder name, not file name
  echo "Generating script for sample: ${sample}"
  mkdir -p ${outDir}/${sample}
  script="${outDir}/${sample}/${sample}_check_ins_nodes.sh"

  cat <<EOF > "$script"
#!/bin/bash
#SBATCH --job-name=check_${sample}_INS
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=${outDir}/${sample}/check_${sample}_INS_%j.log

source /home/qwan/.bashrc
module load Mamba/24.3.0-0
mamba activate /home/qwan/miniconda3/envs/coh

# python3 /home/qwan/githubRepo/coh_PMM/WGBS_scripts/methylGrapher_T2T/panGenome_hg38_INS_flank_DNAm/select_ins_node_gaf.py \\
#   --gaf ${gaf} \\
#   --node_map /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/methylGrapher_ref/pangenomeInsertion/pan_ins_CHM13graph_Length50_insID_nodes_map.txt \\
#   --out ${outDir}/${sample}/${sample}_ins_reads_nodes.txt \\
#   --sample ${sample}

# count total unique reads in GAF
cut -f1 ${gaf} | sort | uniq | wc -l > ${outDir}/${sample}/${sample}_total_reads.txt

EOF

  chmod +x "$script"
done









