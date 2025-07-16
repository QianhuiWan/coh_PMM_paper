#!/bin/bash
#SBATCH --job-name=TEtx_benchmark_simulatedData_S1to4_R2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 8 # Number of cores
#SBATCH -N 1 # Min - Max Nodes
#SBATCH -p all
#SBATCH --mem=192G
#SBATCH --time=72:00:00 #72
#SBATCH --output=TEtx_benchmark_simulatedData_S1to4_R2_%j.log

# ==============================
# TE Transcript Identification Benchmarking (Simulated Data)
# Author: Qianhui Wan
# Updated: 2025 March - Paired-end + RF stranded simulation
# use gene + TE annotation
# use this round: R2 as filnal script for data simulation
# ==============================

source /home/qwan/.bashrc
conda activate /home/qwan/miniconda3/envs/coh

echo "=== Conda Environment Info ==="
# conda list | head -n 10
which python

# === Paths ===
# REF_GENOME="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_hg38_p14.fa"
# REF_GENOME="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_hg38_p14.upper.fa"
REF_GENOME="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/data/example/HomoSapiens_GRCh38_Ensemblv99/HomoSapiens_GRCh38_Ensemblv99.oneline_seqs.fa"
REPEAT_TSV="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_repeatMasker_hg38_addTEid.tsv"
REAL_R1_FASTQ="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/PMM7/PMM7_R1.fastq"
REAL_R2_FASTQ="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/PMM7/PMM7_R2.fastq"

OUT_DIR="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2"
mkdir -p ${OUT_DIR} ${OUT_DIR}/data ${OUT_DIR}/simulated_TEreads/camparee_output ${OUT_DIR}/results_simulated

REPEAT_ANNOTATION=${OUT_DIR}/data/filtered_repeat_annotation.gtf

# ---------------------------------
# Step 1: Prepare Repeat Annotation
# ---------------------------------
echo "Step 1: Preparing Repeat Annotation..."

if [ ! -f "${REPEAT_ANNOTATION}" ]; then
    if [ -f "${REPEAT_TSV}" ]; then
        echo "Converting repeats.tsv to GTF format..."
        awk 'BEGIN{OFS="\t"} NR>1 {print $1, "RepeatMasker", "exon", $2, $3, ".", $4, ".", "gene_id \""$5"\"; transcript_id \""$5"\"; repName \""$6"\"; repClass \""$7"\"; repFamily \""$8"\";"}' ${REPEAT_TSV} > ${REPEAT_ANNOTATION}
    else
        echo "Error: repeats.tsv not found! Please provide a valid TE annotation file."
        exit 1
    fi
else
    echo "Repeat annotation already exists."
fi

# awk 'BEGIN {OFS="\t"}
#      NR > 1 { 
#        chrom = $1;
#        start = $2;
#        end = $3;
#        strand = $4;
#        te_id = $5;
#        repName = $6;
#        exonCount = 1;
#        exonStarts = start;
#        exonEnds = end;
#        transcriptID = te_id;
#        geneID = te_id;
#        geneSymbol = repName;
#        biotype = "repeat_element";
#        print chrom, strand, start, end, exonCount, exonStarts, exonEnds, transcriptID, geneID, geneSymbol, biotype;
#      }' ${REPEAT_TSV} > ${OUT_DIR}/data/repeat_annotation.camparee.txt

# cat /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/data/example/HomoSapiens_GRCh38_Ensemblv99/HomoSapiens_GRCh38_Ensemblv99.annotation.txt ${OUT_DIR}/data/repeat_annotation.camparee.txt > /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/data/example/HomoSapiens_GRCh38_Ensemblv99/HomoSapiens_GRCh38_Ensemblv99.annotation.addTE.txt

# sort -k1,1 -k3,3n /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/data/example/HomoSapiens_GRCh38_Ensemblv99/HomoSapiens_GRCh38_Ensemblv99.annotation.addTE.txt > /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/data/example/HomoSapiens_GRCh38_Ensemblv99/HomoSapiens_GRCh38_Ensemblv99.annotation.addTE.sorted.txt

# ------------------------------
# Step 2: Generate config YAML for CAMPAREE, and then generate RNA molecules
# ------------------------------
echo "Step 2: Generating config YAML for CAMPAREE..."

cat <<EOF > ${OUT_DIR}/simulated_TEreads/camparee_simulationConfig.yaml
setup:
    run_id: pmm7_sim_R2
    seed: 101
    scheduler_mode: serial
    default_scheduler_parameters:
        default_num_processors: 1
        default_memory_in_mb: 6000
    job_resub_limit: 3

resources:
    directory_path: /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/data/example/
    species_model: HomoSapiens_GRCh38_Ensemblv99
    star_genome_index_directory_name: star_index.genome
    reference_genome_filename: HomoSapiens_GRCh38_Ensemblv99.oneline_seqs.fa
    annotation_filename: HomoSapiens_GRCh38_Ensemblv99.annotation.addTE.sorted.txt
    chr_ploidy_filename: HomoSapiens_GRCh38_Ensemblv99.chr_ploidy.txt

input:
    fastq_directory_path: /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/allFastqs_rna/
    optional_inputs:

    data:
        'sample1':
            fastq_files:
                - PMM7_S4_R1_001.fastq.gz
                - PMM7_S4_R2_001.fastq.gz
            pooled: False
            gender: female
            molecule_count: 50000
            optional_inputs:
        'sample2':
            fastq_files:
                - PMM7_S4_R1_001.fastq.gz
                - PMM7_S4_R2_001.fastq.gz
            pooled: False
            gender: female
            molecule_count: 50000
            optional_inputs:
        'sample3':
            fastq_files:
                - PMM7_S4_R1_001.fastq.gz
                - PMM7_S4_R2_001.fastq.gz
            pooled: False
            gender: female
            molecule_count: 50000
            optional_inputs:

output:
    directory_path: ${OUT_DIR}/simulated_TEreads/camparee_output/
    type: packet
    override_sample_molecule_count: True
    default_molecule_count: 50000
    parameters:
        min_polyA_tail_length: 50
        max_polyA_tail_length: 250
    scheduler_parameters:
        memory_in_mb: 8000

steps:
    'genome_alignment.GenomeAlignmentStep':
        parameters:
            '--runThreadN': 4
            '--readFilesCommand': zcat
        scheduler_parameters:
            num_processors: 4
            memory_in_mb: 40000
    'genome_alignment.GenomeBamIndexStep':
    'intron_quant.IntronQuantificationStep':
        parameters:
            forward_read_is_sense: false
            flank_size: 1500
    'variants_finder.VariantsFinderStep':
        parameters:
            sort_by_entropy: false
            min_threshold: 0.03
        scheduler_parameters:
            memory_in_mb: 12000
    'variants_compilation.VariantsCompilationStep':
    'beagle.BeagleStep':
        parameters:
            nthreads: 1
    'genome_builder.GenomeBuilderStep':
        parameters:
            ignore_indels: false
            ignore_snps: false
    'update_annotation_for_genome.UpdateAnnotationForGenomeStep':
    'transcriptome_fasta_preparation.TranscriptomeFastaPreparationStep':
    'kallisto.KallistoIndexStep':
    'kallisto.KallistoQuantStep':
    'bowtie2.Bowtie2IndexStep':
        parameters:
            num_bowtie_threads: 7
        scheduler_parameters:
            num_processors: 7
            memory_in_mb: 40000
    'bowtie2.Bowtie2AlignStep':
        parameters:
            num_bowtie_threads: 7
        scheduler_parameters:
            num_processors: 7
            memory_in_mb: 40000
    'transcript_gene_quant.TranscriptGeneQuantificationStep':
    'allelic_imbalance_quant.AllelicImbalanceQuantificationStep':
        scheduler_parameters:
            memory_in_mb: 192000
EOF

echo "Launching CAMPAREE..."
## pmm7_sim_R2 is the round with outputs that are ok to use for beers2
## add coh env path the first, temporarily
export PATH="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/miniconda3/envs/beers2_env/bin:\$PATH"
# camparee -c ${OUT_DIR}/simulated_TEreads/comparee_simulationConfig.yaml -r pmm7_sim_R2 -m serial -s 52 -d

## pmm7_sim is only a test
# camparee -c ${OUT_DIR}/simulated_TEreads/comparee_simulationConfig.yaml -r pmm7_sim -m slurm -s 52 -d
echo "CAMPAREE step finished"

echo "then dump all pickles into a molecule_file.txt for each simulated sample"
# python3 /home/qwan/githubRepo/coh_PMM/RNAseq_scripts/benchmark/dump_pickles_to_txt.py \
#     -i /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/simulated_TEreads/camparee_output/run_pmm7_sim_R2/CAMPAREE/data/sample1 \
#     -o /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/simulated_TEreads/camparee_output/run_pmm7_sim_R2/CAMPAREE/data/sample1/molecule_file.txt

# python3 /home/qwan/githubRepo/coh_PMM/RNAseq_scripts/benchmark/dump_pickles_to_txt.py \
#     -i /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/simulated_TEreads/camparee_output/run_pmm7_sim_R2/CAMPAREE/data/sample2 \
#     -o /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/simulated_TEreads/camparee_output/run_pmm7_sim_R2/CAMPAREE/data/sample2/molecule_file.txt

# python3 /home/qwan/githubRepo/coh_PMM/RNAseq_scripts/benchmark/dump_pickles_to_txt.py \
#     -i /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/simulated_TEreads/camparee_output/run_pmm7_sim_R2/CAMPAREE/data/sample3 \
#     -o /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/simulated_TEreads/camparee_output/run_pmm7_sim_R2/CAMPAREE/data/sample3/molecule_file.txt
echo "finished dumping all pickles into a molecule_file.txt for each simulated sample"

# ------------------------------
# Step 3: Generate config YAML for BEERS2
# ------------------------------
echo "Step 3: Generating config YAML for BEERS2..."

cat <<EOF > ${OUT_DIR}/simulated_TEreads/beers2_simulationConfig.yaml
seed: 52

output:
    output_fastq: true
    output_sam: false
    output_bam: true
    full_logs: true

global_config:
    samples:
        '1':
            barcodes:
                i5: AGCGCTAG
                i7: AACCGCGG
        '2':
            barcodes:
                i5: GATATCGA
                i7: TTATAACC
        '3':
            barcodes:
                i5: TTCCGGAA
                i7: CCGGTTAA

    molecule_maker_parameters:
        min_polyA_tail_length: 50
        max_polyA_tail_length: 250

    resources:
        pre_i5_adapter: AATGATACGGCGACCACCGAGATCTACAC
        post_i5_adapter: ACACTCTTTCCCTACACGACGCTCTTCCGATCT
        pre_i7_adapter: GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
        post_i7_adapter: ATCTCGTATGCCGTCTTCTGCTTG
        reference_genome_fasta: ${REF_GENOME}

# probably need num_packets: 100 , to generate 20 million reads
library_prep_pipeline:
    input:
        directory_path: ${OUT_DIR}/simulated_TEreads/camparee_output/run_pmm7_sim_R2/CAMPAREE/data/
        from_distribution_data:
            '1':
                num_packets: 100
                num_molecules_per_packet: 200000
                sample_data_directory: ${OUT_DIR}/simulated_TEreads/camparee_output/run_pmm7_sim_R2/CAMPAREE/data/sample1/
            '2':
                num_packets: 100
                num_molecules_per_packet: 200000
                sample_data_directory: ${OUT_DIR}/simulated_TEreads/camparee_output/run_pmm7_sim_R2/CAMPAREE/data/sample2/
            '3':
                num_packets: 100
                num_molecules_per_packet: 200000
                sample_data_directory: ${OUT_DIR}/simulated_TEreads/camparee_output/run_pmm7_sim_R2/CAMPAREE/data/sample3/

    steps:
    -   step_name: polya_step.PolyAStep
        parameters:
            breakpoint_prob_per_base: 0.001
            max_retention_prob: 1.0
            min_retention_prob: 0.0
            min_polya_tail_length: 40
            length_retention_prob: 0.05

    -   step_name: fragment_step.FragmentStep
        parameters:
            method: uniform
            rate: 0.005
            runtime: 1
            min_frag_size: 20

    -   step_name: first_strand_synthesis_step.FirstStrandSynthesisStep
        parameters:
            perfect_priming: false
            position_probability_matrix:
                A: [0.50, 0.1, 0.40, 0.30, 0.25, 0.15]
                C: [0.20, 0.5, 0.3 , 0.25, 0.25, 0.15]
                G: [0.15, 0.1, 0.15, 0.25, 0.25, 0.20]
                T: [0.15, 0.3, 0.15, 0.20, 0.25, 0.50]
            primes_per_kb: 50

    -   step_name: second_strand_synthesis_step.SecondStrandSynthesisStep
        parameters:
            perfect_priming: false
            position_probability_matrix:
                A: [0.50, 0.1, 0.40, 0.30, 0.25, 0.15]
                C: [0.20, 0.5, 0.3 , 0.25, 0.25, 0.15]
                G: [0.15, 0.1, 0.15, 0.25, 0.25, 0.20]
                T: [0.15, 0.3, 0.15, 0.20, 0.25, 0.50]
            primes_per_kb: 50

    -   step_name: sizing_step.SizingStep
        parameters:
            min_length: 100
            max_length: 400
            select_all_start_length: 200
            select_all_end_length: 300

    -   step_name: adapter_ligation_step.AdapterLigationStep
        parameters: {}

    -   step_name: pcr_amplification_step.PCRAmplificationStep
        parameters:
            number_cycles: 10
            retention_percentage: 0.08
            gc_bias_constant: 1.0
            gc_bias_linear: 0.0
            gc_bias_quadratic: -100
            deletion_rate: 0.0001
            insertion_rate: 0.0001
            substitution_rate: 0.001

sequence_pipeline:
    steps:
    -   step_name: bridge_amplification_step.BridgeAmplificationStep
        parameters:
            cycles: 10
            substitution_rate: 0.01

    -   step_name:  sequence_by_synthesis_step.SequenceBySynthesisStep
        parameters:
            forward_is_5_prime: true
            paired_ends: true
            read_length: 150
            skip_rate: 0.002
            drop_rate: 0.002

    flowcell:
        coordinate_strategy: random
        flowcell_geometry:
            min_lane: 1
            max_lane: 8
            min_tile: 1101
            max_tile: 2228
            min_x: 1012
            max_x: 32816
            min_y: 998
            max_y: 49247
        lanes_to_use: [1, 2]
EOF

# ------------------------------
# Step 4: Run BEERS2 simulation, simulate 20 million reads from 50000 rna molecues
# ------------------------------
echo "Step 4: Running BEERS2 under results_simulated folder..."
conda activate /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/miniconda3/envs/beers2_env
# beers --configfile ${OUT_DIR}/simulated_TEreads/beers2_simulationConfig.yaml \
#     --profile /home/qwan/slurm \
#     --unlock

beers --configfile ${OUT_DIR}/simulated_TEreads/beers2_simulationConfig.yaml \
    --jobs 8 \
    --profile /home/qwan/slurm

echo "Simulation complete! Results saved in: ${OUT_DIR}/results_simulated"
conda deactivate





