#!/bin/bash

## variables
samples=$1

# read each line from the $1 input sampleNames.txt file and create .sh script for each sample
while IFS= read -r sample_name; do
  echo "Creating bash script for sample: ${sample_name}"
  short_sample_name=${sample_name/sample/S}
  echo "Short sample name is: ${short_sample_name}"
  # Define the directory:
  input_path=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/results_simulated/results

  mkdir -p /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/results_simulated/results/mergedFiles/${sample_name}
  mkdir -p /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/results_simulated/results/tetx_detect/${sample_name}

  mergedFile_path=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/results_simulated/results/mergedFiles/${sample_name}
  tetx_output_path=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/results_simulated/results/tetx_detect/${sample_name}
  teprof3_star_path=${tetx_output_path}/teprof3/star_output
  # software:
  SalmonTE=/home/qwan/nonCondaPkg/SalmonTE/SalmonTE.py
  STAR=/home/qwan/miniconda3/envs/coh/bin/STAR
  java=/home/qwan/miniconda3/envs/coh/bin/java
  teprof3=/home/qwan/nonCondaPkg/TEProf3/bin/teprof3

  cat <<EOF > "${tetx_output_path}/${sample_name}_testTetxMethods.sh"
#!/bin/bash
#SBATCH --job-name=TEtx_benchmark_simulatedData_S5_6_${sample_name}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 8 # Number of cores
#SBATCH -N 1 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=64G
#SBATCH --time=12:00:00 #24h
#SBATCH --output=TEtx_benchmark_simulatedData_S5_6_${sample_name}_%j.log

# ==============================
# TE Transcript Identification Benchmarking (Simulated Data)
# Author: Qianhui Wan
# Date: 2025 March
# ==============================

source /home/qwan/.bashrc
module load Mamba/24.3.0-0
mamba activate /home/qwan/miniconda3/envs/coh

# ------------------------------
# Step 5: merge L1 and L2
# ------------------------------
echo "Step 5: Merge files from 2 allels i.e. L1 and L2, for each sample..."

## for bam files -- add RG lines, read group information:
# samtools addreplacerg -r "ID:L1\tSM:${sample_name}\tLB:simulatedLib_beers2\tPL:ILLUMINA_beers2" -o ${mergedFile_path}/L1_rg.bam ${input_path}/${short_sample_name}_L1.bam
# samtools addreplacerg -r "ID:L2\tSM:${sample_name}\tLB:simulatedLib_beers2\tPL:ILLUMINA_beers2" -o ${mergedFile_path}/L2_rg.bam ${input_path}/${short_sample_name}_L2.bam

# samtools merge -@ 8 ${mergedFile_path}/merged.bam ${mergedFile_path}/L1_rg.bam ${mergedFile_path}/L2_rg.bam
# samtools sort -@ 8 -o ${mergedFile_path}/merged_sort.bam ${mergedFile_path}/merged.bam
# samtools index -@ 8 ${mergedFile_path}/merged_sort.bam

## for fastq files -- cat them together:
# cat ${input_path}/${short_sample_name}_L1_R1.fastq ${input_path}/${short_sample_name}_L2_R1.fastq > ${mergedFile_path}/${sample_name}_R1.fastq

# cat ${input_path}/${short_sample_name}_L1_R2.fastq ${input_path}/${short_sample_name}_L2_R2.fastq > ${mergedFile_path}/${sample_name}_R2.fastq


# ------------------------------
# Step 6: Run TE Identification Tools
# ------------------------------
echo "Step 6: Running TE Identification Tools..."

# 6.1 SalmonTE
echo "Running SalmonTE..."
## regenerate output if error log
# rm -r ${tetx_output_path}/salmonTE/*

# mkdir -p ${tetx_output_path}/salmonTE
# ln -s ${mergedFile_path}/${sample_name}_R1.fastq ${tetx_output_path}/salmonTE/${sample_name}_R1.fastq
# ln -s ${mergedFile_path}/${sample_name}_R2.fastq ${tetx_output_path}/salmonTE/${sample_name}_R2.fastq

## add coh path the first, temporarily
# export PATH="/home/qwan/miniconda3/envs/coh/bin:\$PATH"

# python3 ${SalmonTE} quant --reference=hs \
# --outpath=${tetx_output_path}/salmonTE \
# --num_threads=8 \
# ${tetx_output_path}/salmonTE

# 6.2 TEprof3
echo "Running TEprof3..."
## regenerate output if error log
# rm -r ${tetx_output_path}/teprof3/*

## make teprof3_star_path:
# mkdir -p ${tetx_output_path}/teprof3
# mkdir -p ${tetx_output_path}/teprof3/star_output

## 6.2.1 star alignment
# ln -s ${mergedFile_path}/${sample_name}_R1.fastq ${teprof3_star_path}/${sample_name}_R1.fastq
# ln -s ${mergedFile_path}/${sample_name}_R2.fastq ${teprof3_star_path}/${sample_name}_R2.fastq

# ${STAR} --genomeDir /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/STAR_hg38_p14_geneCodeGTF_filter \
# --readFilesIn ${teprof3_star_path}/${sample_name}_R1.fastq ${teprof3_star_path}/${sample_name}_R2.fastq \
# --runThreadN 8 \
# --twopassMode Basic \
# --outFilterMultimapNmax 20 \
# --alignSJoverhangMin 8 \
# --alignSJDBoverhangMin 1 \
# --outFilterMismatchNmax 999 \
# --outFilterMismatchNoverLmax 0.1 \
# --alignIntronMin 20 \
# --alignIntronMax 1000000 \
# --alignMatesGapMax 1000000 \
# --outFilterType BySJout \
# --outFilterScoreMinOverLread 0.33 \
# --outFilterMatchNminOverLread 0.33 \
# --limitSjdbInsertNsj 1200000 \
# --outFileNamePrefix ${teprof3_star_path}/${sample_name}_ \
# --outSAMstrandField intronMotif \
# --outFilterIntronMotifs None \
# --alignSoftClipAtReferenceEnds Yes \
# --quantMode TranscriptomeSAM GeneCounts \
# --outSAMtype BAM Unsorted \
# --outSAMunmapped Within \
# --genomeLoad NoSharedMemory \
# --chimSegmentMin 15 \
# --chimJunctionOverhangMin 15 \
# --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
# --chimOutJunctionFormat 1 \
# --chimMainSegmentMultNmax 1 \
# --outSAMattributes NH HI AS nM NM ch

# --outFilterMultimapNmax

### sort by coordinates, Remove duplicates and sort
# samtools sort -@ 8 -O bam -o ${teprof3_star_path}/${sample_name}_sorted.bam ${teprof3_star_path}/${sample_name}_Aligned.out.bam
# samtools index ${teprof3_star_path}/${sample_name}_sorted.bam
# rm ${teprof3_star_path}/${sample_name}_Aligned.out.bam

# ${java} -Djava.io.tmpdir=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/temp \
# -jar /home/qwan/miniconda3/envs/coh/share/picard-2.27.5-0/picard.jar MarkDuplicates \
# -INPUT ${teprof3_star_path}/${sample_name}_sorted.bam \
# -OUTPUT ${teprof3_star_path}/${sample_name}_nr_sorted.bam \
# -REMOVE_DUPLICATES true -METRICS_FILE ${teprof3_star_path}/${sample_name}_picardStats.txt

# samtools sort -@ 8 -O bam -o ${teprof3_star_path}/${sample_name}_sorted_nr_sorted.bam ${teprof3_star_path}/${sample_name}_nr_sorted.bam
# samtools index -@ 8 ${teprof3_star_path}/${sample_name}_sorted_nr_sorted.bam

# ------------------------------
# Step 6.2.2: split waston and crick strand bam files
# ------------------------------
echo "Step 6.2.2: Split waston and crick strand bam files..."

# samtools view -@ 8 -h -f 99 -S -b ${teprof3_star_path}/${sample_name}_sorted_nr_sorted.bam > ${teprof3_star_path}/crick1.bam
# samtools view -@ 8 -h -f 147 -S -b ${teprof3_star_path}/${sample_name}_sorted_nr_sorted.bam > ${teprof3_star_path}/crick2.bam
# samtools view -@ 8 -h -f 83 -S -b ${teprof3_star_path}/${sample_name}_sorted_nr_sorted.bam > ${teprof3_star_path}/watson1.bam
# samtools view -@ 8 -h -f 163 -S -b ${teprof3_star_path}/${sample_name}_sorted_nr_sorted.bam > ${teprof3_star_path}/watson2.bam

# samtools sort -@ 8 -o ${teprof3_star_path}/crick1_sorted.bam ${teprof3_star_path}/crick1.bam
# samtools sort -@ 8 -o ${teprof3_star_path}/crick2_sorted.bam ${teprof3_star_path}/crick2.bam
# samtools sort -@ 8 -o ${teprof3_star_path}/watson1_sorted.bam ${teprof3_star_path}/watson1.bam
# samtools sort -@ 8 -o ${teprof3_star_path}/watson2_sorted.bam ${teprof3_star_path}/watson2.bam

# samtools merge -f ${teprof3_star_path}/${sample_name}_crick.bam \
# ${teprof3_star_path}/crick1_sorted.bam ${teprof3_star_path}/crick2_sorted.bam

# samtools merge -f ${teprof3_star_path}/${sample_name}_watson.bam \
# ${teprof3_star_path}/watson1_sorted.bam ${teprof3_star_path}/watson2_sorted.bam

# samtools index ${teprof3_star_path}/${sample_name}_crick.bam
# samtools index ${teprof3_star_path}/${sample_name}_watson.bam

# rm ${teprof3_star_path}/crick1.bam ${teprof3_star_path}/crick2.bam ${teprof3_star_path}/watson1.bam ${teprof3_star_path}/watson2.bam
# rm ${teprof3_star_path}/crick1_sorted.bam ${teprof3_star_path}/crick2_sorted.bam \
# ${teprof3_star_path}/watson1_sorted.bam ${teprof3_star_path}/watson2_sorted.bam

## 6.2.3 TEprof3 detection
### regenerate output if error log
# rm -r ${tetx_output_path}/teprof3/teprof3_output/*
### create sub-directory
# mkdir -p ${tetx_output_path}/teprof3/teprof3_output

# cd ${tetx_output_path}/teprof3/teprof3_output

# reset all folders if error log
# ${teprof3} -f sample_manifest.txt --reset

# ln -s ${teprof3_star_path}/${sample_name}_sorted_nr_sorted.bam .
# ln -s ${teprof3_star_path}/${sample_name}_sorted_nr_sorted.bam.bai .
# ln -s ${teprof3_star_path}/${sample_name}_SJ.out.tab .

# get sample_manifest file
# rm sample_manifest.txt
# find . -maxdepth 1 -name "*.bam" | while read file ; do xbase=\$(basename \$file); sample_name=\${xbase/_sorted_nr_sorted.bam/} ; echo -e "${sample_name}\tshort\t\${xbase}" >> sample_manifest.txt; done ;
# find . -maxdepth 1 -name "*.tab" | while read file ; do xbase=\$(basename \$file); sample_name=\${xbase/_SJ.out.tab/} ; echo -e "${sample_name}\tSJ\t\${xbase}" >> sample_manifest.txt; done ;
# find . -maxdepth 1 -name "*.gtf" | while read file ; do xbase=\$(basename \$file); sample_name=\${xbase/.gtf/} ; echo -e "${sample_name}\tgtf\t\${xbase}" >> sample_manifest.txt; done

## add coh env path the first, temporarily
# export PATH="/home/qwan/miniconda3/envs/taco_env/bin:/home/qwan/miniconda3/envs/coh/bin:\$PATH"

# ${teprof3} -f sample_manifest.txt -s 1 -am 1 -as 1 --assemblethread 4 --assemblestrand 1 --assemblejunctionread 2 \
# -ps 1 -pt 0.5 -ptn 100 \
# -fs 1 -fm 1 --filterintronretention 3 \
# --filtermonoexontpm 1 --filterdownstreammate 2 --filterratio 0.5 \
# --tacothread 1 \
# -qs 1 --quansamplenumbercon 100 --quanmode 1 --quanreadlength 50

# 6.3 Activation score
echo "Running Activation score..."
## regenerate output if error log
# rm -r ${tetx_output_path}/activation_score/*
# mkdir -p ${tetx_output_path}/activation_score
# cd ${tetx_output_path}/activation_score
# source ~/githubRepo/coh_PMM/RNAseq_scripts/benchmark/activationScore_functions/activation_score.sh
# activation_score_script ${sample_name} \
# ${teprof3_star_path}/${sample_name}_crick.bam ${teprof3_star_path}/${sample_name}_watson.bam \
# ${tetx_output_path}/activation_score

# for test1: rm --countReadPairs, not work 
# mkdir -p ${tetx_output_path}/activation_score_countReadsTest
# cd ${tetx_output_path}/activation_score_countReadsTest

# source ~/githubRepo/coh_PMM/RNAseq_scripts/benchmark/activationScore_functions/activation_score.sh
# activation_score_script ${sample_name} \
# ${teprof3_star_path}/${sample_name}_crick.bam ${teprof3_star_path}/${sample_name}_watson.bam \
# ${tetx_output_path}/activation_score_countReadsTest

# for test2: add -O add --countReadPairs
mkdir -p ${tetx_output_path}/activation_score_addOtest
cd ${tetx_output_path}/activation_score_addOtest

source ~/githubRepo/coh_PMM/RNAseq_scripts/benchmark/activationScore_functions/activation_score_o.sh
activation_score_script ${sample_name} \
${teprof3_star_path}/${sample_name}_crick.bam ${teprof3_star_path}/${sample_name}_watson.bam \
${tetx_output_path}/activation_score_addOtest

mamba deactivate

EOF
echo "Script created: ${tetx_output_path}/${sample_name}_testTetxMethods.sh"
#  cd ${tetx_output_path}
#  sbatch ${sample_name}_testTetxMethods.sh
done < ${samples}


echo "All sample script files created successfully."





