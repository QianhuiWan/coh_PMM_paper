#!/bin/bash

# use TCGA 2passSTAR alignment and use hg38.p14 without patches
# S1 round2: remove multimapping in S1 Round2 to reduce simple repeat alignment
# S1 round3: round2+ remove low quality and short reads in `fastq` step to reduce simple repeat alignment
# S1 round4: round2 + remove fastq step, use original fq for alignment

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
  ## we can soft-link all data to a data folder, so we can use it more conveniently
  datapath_MM=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/allFastqs_MMRF_rna
  mkdir -p /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s1_checkStrandness_R4
  mkdir -p /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s1_checkStrandness_R4/${sample_name}
  outdir=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s1_checkStrandness_R4/${sample_name}
  # software
  fastqc=/home/qwan/miniconda3/envs/coh/bin/fastqc
  fastp=/home/qwan/miniconda3/envs/coh/bin/fastp
  STAR=/home/qwan/miniconda3/envs/coh/bin/STAR
  allo=/home/qwan/miniconda3/envs/coh/bin/allo
  samtools=/home/qwan/miniconda3/envs/coh/bin/samtools
  picard=/home/qwan/miniconda3/envs/coh/bin/picard
  java=/home/qwan/miniconda3/envs/coh/bin/java
  bedtools=/home/qwan/miniconda3/envs/coh/bin/bedtools
  wigToBigWig=/home/qwan/miniconda3/envs/coh/bin/wigToBigWig
  bamCoverage=/home/qwan/miniconda3/envs/coh/bin/bamCoverage
  {
    echo -e "#!/bin/bash
#SBATCH --job-name=${sample_name}_checkStrandness_R4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 12 # Number of cores
#SBATCH -N 1  # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=${outdir}/${sample_name}_checkStrandness_R4_%j.log\n

source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/coh

#### fastqc for fastq files
mkdir -p ${outdir}/fastqc_out
${fastqc} -t 8 -o ${outdir}/fastqc_out \
${datapath_MM}/${sample_name}_1*.gz ${datapath_MM}/${sample_name}_2*.gz


### Align reads with STAR:

${STAR} --genomeDir /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/STAR_hg38_p14_geneCodeGTF_filter \
--readFilesIn ${datapath_MM}/${sample_name}_1*.gz ${datapath_MM}/${sample_name}_2*.gz \
--readFilesCommand zcat \
--runThreadN 8 \
--twopassMode Basic \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFilterType BySJout \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNminOverLread 0.33 \
--limitSjdbInsertNsj 1200000 \
--outFileNamePrefix ${outdir}/${sample_name}_ \
--outSAMstrandField intronMotif \
--outFilterIntronMotifs None \
--alignSoftClipAtReferenceEnds Yes \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--genomeLoad NoSharedMemory \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
--chimOutJunctionFormat 1 \
--chimMainSegmentMultNmax 1 \
--outSAMattributes NH HI AS nM NM ch

#--outSAMattrRGline ID:sample1 SM:sample1 PL:ILLUMINA LB:lib1 PU:unit1; not sure lib info. so leave it out
#--outFilterMismatchNmax 10, not use in R1 and R2 because:
# other para. covered this filtering, the default value is unlimited
#--outFilterMultimapNmax 100, for repeat, not use in this round
#--winAnchorMultimapNmax 100, for repeat, not use in this round
# --outSAMmultNmax 25, for allo, not use in this round

### sort by coordinates, Remove duplicates and sort
${samtools} sort -@ 8 -O bam -o ${outdir}/${sample_name}_sorted.bam ${outdir}/${sample_name}_Aligned.out.bam
${samtools} index -@ 8 ${outdir}/${sample_name}_sorted.bam
rm ${outdir}/${sample_name}_Aligned.out.bam

${java} -Djava.io.tmpdir=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/temp \
-jar /home/qwan/miniconda3/envs/coh/share/picard-2.27.5-0/picard.jar MarkDuplicates \
-INPUT ${outdir}/${sample_name}_sorted.bam -OUTPUT ${outdir}/${sample_name}_nr_sorted.bam \
-REMOVE_DUPLICATES true -METRICS_FILE ${outdir}/${sample_name}_picardStats.txt

${samtools} sort -@ 8 -O bam -o ${outdir}/${sample_name}_sorted_nr_sorted.bam ${outdir}/${sample_name}_nr_sorted.bam
${samtools} index -@ 8 ${outdir}/${sample_name}_sorted_nr_sorted.bam

# check for strand info. of RNAseq data here, if stranded, get separate bam for watson and crick
module load RSeQC/4.0.0
cd ${outdir}
infer_experiment.py -r /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/teAnno_round3/hg38_gencode_exon.bed \
-i ${outdir}/${sample_name}_sorted_nr_sorted.bam

conda deactivate"
  } > ${outdir}/${sample_name}_rnaPreprocess_s1.sh
#  cd ${outdir}
#  sbatch ${sample_name}_rnaPreprocess.sh
done < ${samples}

echo "All sample script files created successfully."
