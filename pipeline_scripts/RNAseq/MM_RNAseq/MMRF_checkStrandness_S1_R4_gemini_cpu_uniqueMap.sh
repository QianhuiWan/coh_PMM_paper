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
  datapath_MM=/scratch/qwan/projects/PMM/rawData/allFastqs_MMRF_rna
  mkdir -p /scratch/qwan/projects/PMM/RNA_MMRF/s1_checkStrandness_R4_uniqueMap
  mkdir -p /scratch/qwan/projects/PMM/RNA_MMRF/s1_checkStrandness_R4_uniqueMap/${sample_name}
  outdir=/scratch/qwan/projects/PMM/RNA_MMRF/s1_checkStrandness_R4_uniqueMap/${sample_name}
  # software
  fastqc=/home/qwan/miniconda3/envs/coh_PMM/bin/fastqc
  fastp=/home/qwan/miniconda3/envs/coh_PMM/bin/fastp
  STAR=/home/qwan/miniconda3/envs/coh_PMM/bin/STAR
  allo=/home/qwan/miniconda3/envs/coh_PMM/bin/allo
  samtools=/home/qwan/miniconda3/envs/coh_PMM/bin/samtools
  picard=/home/qwan/miniconda3/envs/coh_PMM/bin/picard
  java=/home/qwan/miniconda3/envs/coh_PMM/bin/java
  bedtools=/home/qwan/miniconda3/envs/coh_PMM/bin/bedtools
  wigToBigWig=/home/qwan/miniconda3/envs/coh_PMM/bin/wigToBigWig
  bamCoverage=/home/qwan/miniconda3/envs/coh_PMM/bin/bamCoverage
  {
    echo -e "#!/bin/bash
#SBATCH --job-name=${sample_name}_checkStrandness_R4_uniqueMap
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=qwan@coh.org      # Where to send mail
#SBATCH -n 10                         # Number of cores
#SBATCH -N 1                          # Min - Max Nodes
#SBATCH -p compute                    # gpu queue. Replace with gpu-v100/gpu-a100 to run on V100 cards
#SBATCH --mem=60G                     # Amount of memory in GB
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=${outdir}/${sample_name}_checkStrandness_R4_uniqueMap_%j.log

# miniconda3 on gemini
source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/coh_PMM

#### fastqc for fastq files
# mkdir -p ${outdir}/fastqc_out
# ${fastqc} -t 8 -o ${outdir}/fastqc_out \
# ${datapath_MM}/${sample_name}_1*.gz ${datapath_MM}/${sample_name}_2*.gz


### Align reads with STAR:

${STAR} --genomeDir /scratch/qwan/projects/PMM/RNA_MMRF/annoFiles/STAR_hg38_p14_geneCodeGTF_filter \
--readFilesIn ${datapath_MM}/${sample_name}_1*.gz ${datapath_MM}/${sample_name}_2*.gz \
--readFilesCommand zcat \
--runThreadN 8 \
--twopassMode Basic \
--outFilterMultimapNmax 1 \
--outSAMmultNmax 1 \
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
# --outFilterMultimapNmax 1; --outSAMmultNmax 1; keep unique mapped reads only

### sort by coordinates, Remove duplicates and sort
${samtools} sort -@ 8 -O bam -o ${outdir}/${sample_name}_sorted.bam ${outdir}/${sample_name}_Aligned.out.bam
${samtools} index -@ 8 ${outdir}/${sample_name}_sorted.bam
rm ${outdir}/${sample_name}_Aligned.out.bam

${java} -Djava.io.tmpdir=/scratch/qwan/temp \
-jar /home/qwan/miniconda3/envs/coh_PMM/share/picard-2.27.5-0/picard.jar MarkDuplicates \
-INPUT ${outdir}/${sample_name}_sorted.bam -OUTPUT ${outdir}/${sample_name}_nr_sorted.bam \
-REMOVE_DUPLICATES true -METRICS_FILE ${outdir}/${sample_name}_picardStats.txt

${samtools} sort -@ 8 -O bam -o ${outdir}/${sample_name}_sorted_nr_sorted.bam ${outdir}/${sample_name}_nr_sorted.bam
${samtools} index -@ 8 ${outdir}/${sample_name}_sorted_nr_sorted.bam

# check for strand info. of RNAseq data here, if stranded, get separate bam for watson and crick
#module load RSeQC/4.0.0
cd ${outdir}
infer_experiment.py -r /scratch/qwan/projects/PMM/RNA_MMRF/annoFiles/hg38_gencode_exon.bed \
-i ${outdir}/${sample_name}_sorted_nr_sorted.bam

### bam to bedGraph
${bedtools} genomecov -bg -ibam ${outdir}/${sample_name}_sorted_nr_sorted.bam \
-g /scratch/qwan/projects/PMM/RNA_MMRF/annoFiles/filtered_hg38_p14.chrom.sizes \
> ${outdir}/${sample_name}_readCount.bedGraph

### bedGraph to bigwig (read counts)
${wigToBigWig} ${outdir}/${sample_name}_readCount.bedGraph \
/scratch/qwan/projects/PMM/RNA_MMRF/annoFiles/filtered_hg38_p14.chrom.sizes \
${outdir}/${sample_name}_readCount.bw

### bigwig (cpm)
${bamCoverage} -b ${outdir}/${sample_name}_sorted_nr_sorted.bam \
-o ${outdir}/${sample_name}_cpm.bw --normalizeUsing CPM \
--numberOfProcessors 10 \
--effectiveGenomeSize 2913022398 --binSize 1

### bigwig (rpkm)
${bamCoverage} -b ${outdir}/${sample_name}_sorted_nr_sorted.bam \
-o ${outdir}/${sample_name}_rpkm.bw --normalizeUsing RPKM \
--numberOfProcessors 10 \
--effectiveGenomeSize 2913022398 --binSize 1

rm ${outdir}/*.bedGraph

if [ -f ${outdir}/${sample_name}_sorted_nr_sorted.bam ]; then
    rm ${datapath_MM}/${sample_name}_1*.gz ${datapath_MM}/${sample_name}_2*.gz
fi

echo \"Preprocessing FINISHED\"

conda deactivate"
  } > ${outdir}/${sample_name}_rnaPreprocess_s1.sh
#  cd ${outdir}
#  sbatch ${sample_name}_rnaPreprocess.sh
done < ${samples}

echo "All sample script files created successfully."




