
# for fig2D and other DNAm, normalize DNAm values for B_rest and PMM
# R3: normalize Bismark output + use cpm before filter

# pkg needed
library(tidyverse)
library(magrittr)
library(plyranges)
library(data.table)
library(methylKit)
library(parallel)
library(readr)
library(GenomicRanges)

outDir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/methylKit_output"


bam_file_list <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/*/03.BismarkPosMethylation/nr_sorted_*_pe.bam")

DNAmBAM_files_sampleInfo <- bam_file_list %>%
  as.data.frame() %>% `colnames<-`("path") %>%
  dplyr::mutate(samplesLabel = str_remove(path, "..+/")) %>%
  dplyr::mutate(samplesLabel = str_remove(samplesLabel, "nr_sorted_")) %>%
  dplyr::mutate(samplesLabel = str_remove(samplesLabel, "_pe.bam")) %>%
  dplyr::mutate(group=str_replace_all(samplesLabel, "[:digit:]", "")) %>%
  dplyr::filter(group != "B_EBV") %>%
  dplyr::mutate(group_integer = ifelse(group=="B_rest", 0, 1))


# Define a wrapper function for a single sample
# run_sample <- function(i) {
#   processBismarkAln(
#     location = DNAmBAM_files_sampleInfo$path[i],
#     sample.id = DNAmBAM_files_sampleInfo$samplesLabel[i],
#     assembly = "hg38",
#     save.folder = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/methylKit_output/methylRaw_output_cov1",
#     save.context = c("CpG"),
#     read.context = "CpG",
#     nolap = FALSE,
#     mincov = 1,
#     minqual = 20,
#     phred64 = FALSE,
#     treatment = DNAmBAM_files_sampleInfo$group_integer[i]
#   )
# }

# Use mclapply to parallelize
# bismark_methylKit_objs <- mclapply(
#   1:nrow(DNAmBAM_files_sampleInfo),
#   run_sample,
#   mc.cores = 6  # Adjust to how many cores you have
# )


file_list <- list.files("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/methylKit_output/methylRaw_output_cov1", 
                        full.names = TRUE, pattern = "CpG.txt$")
sample_ids <- DNAmBAM_files_sampleInfo$samplesLabel

# bismark_methylKit_objs <- methRead(
#   location = as.list(file_list),
#   sample.id = as.list(sample_ids),
#   assembly = "hg38",
#   # treatment = as.numeric(DNAmBAM_files_sampleInfo$group_integer),
#   treatment = rep(1, length(file_list))
# )


# saveRDS(bismark_methylKit_objs, file = paste0(outDir, "/bismark_methylKit_objs_R2_cov1.rds"))
bismark_methylKit_objs <- readRDS(paste0(outDir, "/bismark_methylKit_objs_R2_cov1.rds"))

# tests
getMethylationStats(bismark_methylKit_objs[[1]], plot=TRUE, both.strands=FALSE)
getMethylationStats(bismark_methylKit_objs[[10]], plot=TRUE, both.strands=FALSE)

getCoverageStats(bismark_methylKit_objs[[1]], plot=TRUE, both.strands=FALSE)
getCoverageStats(bismark_methylKit_objs[[10]], plot=TRUE, both.strands=FALSE)

# filterByCoverage is not working well for different batches
# bismark_methylKit_objs_filt <- filterByCoverage(bismark_methylKit_objs,
#                                lo.count=5,
#                                lo.perc=NULL,
#                                hi.count=NULL,
#                                hi.perc=99.9)

# get methylation mat data
# meth <- unite(bismark_methylKit_objs,
#                  destrand = TRUE,
#                  min.per.group = 3L) #13953837 CpG
# saveRDS(meth, file = paste0(outDir, "/bismark_methylKit_cov1_3L_meth.rds"))
meth <- readRDS(paste0(outDir, "/bismark_methylKit_cov1_3L_meth.rds"))

# get coverage
cov_mat <- getData(meth)[, grep("coverage", colnames(getData(meth)))]

# use cpm norm reads cov 
total_cov <- colSums(cov_mat, na.rm = TRUE)
norm_cov <- sweep(cov_mat, 2, total_cov / 1e6, FUN = "/")

# add filter, select CpG with normalized coverage > 1 
# keep <- apply(norm_cov, 1, function(x) all(x > 1))
keep <- apply(norm_cov, 1, function(x) any(x > 1, na.rm = TRUE))

meth_filtered <- meth[keep, ]

# or warp into a R function
#' Normalize coverage and filter CpG sites in methylKit object
#'
#' @param myobj A methylRawList or methylBase object (e.g., from methylKit)
#' @param min_norm_cov Minimum normalized coverage threshold (default = 1)
#' @param min_sample_frac Fraction of samples that must exceed threshold (default = 1, i.e., all samples)
#' @param destrand Whether to destrand the CpG data when uniting (default = FALSE)
#' @return A filtered methylBase object
#' @export
# filterByNormalizedCoverage <- function(myobj, 
#                                        min_norm_cov = 1, 
#                                        min_sample_frac = 1,
#                                        destrand = FALSE) {
#   # Unite methylation data across samples
#   meth <- unite(myobj, destrand = destrand)
#   
#   # Extract raw coverage matrix
#   cov.mat <- getData(meth)[, grep("coverage", colnames(getData(meth)))]
#   
#   # Normalize coverage by total CpGs per sample (in millions)
#   total_cov <- colSums(cov.mat)
#   norm_cov <- sweep(cov.mat, 2, total_cov / 1e6, FUN = "/")
#   
#   # Determine which CpG sites pass the threshold
#   keep <- rowMeans(norm_cov > min_norm_cov) >= min_sample_frac
#   
#   # Return filtered methylBase object
#   meth.filtered <- meth[keep, ]
#   
#   return(meth.filtered)
# }


bismark_methylKit_objs_norm <- 
  normalizeCoverage(bismark_methylKit_objs, method = "median")

meth_norm <- unite(bismark_methylKit_objs_norm,
              destrand = TRUE,
              min.per.group = 3L) #13953837 CpG
saveRDS(meth_norm, file = paste0(outDir, "/bismark_methylKit_cov1_3L_meth_norm.rds"))

# get coverage
cov_mat <- getData(meth_norm)[, grep("coverage", colnames(getData(meth_norm)))]

# use cpm norm reads cov 
total_cov <- colSums(cov_mat, na.rm = TRUE)
norm_cov <- sweep(cov_mat, 2, total_cov / 1e6, FUN = "/")

# add filter, select CpG with normalized coverage > 1 
# keep <- apply(norm_cov, 1, function(x) all(x > 1))
keep <- apply(norm_cov, 1, function(x) any(x > 1, na.rm = TRUE))

meth_norm_filtered <- meth_norm[keep, ]


# saveRDS(bismark_methylKit_objs_filt_norm, file = paste0(outDir, "/bismark_methylKit_objs_filt_norm_R3.rds"))

meth_fil <- unite(bismark_methylKit_objs_filt_norm, destrand = TRUE)
meth_fil
getCorrelation(meth_fil, plot=TRUE)
clusterSamples(meth_fil, dist="correlation", method="ward.D2", plot=TRUE)
PCASamples(meth_fil)


# Read in your chrom sizes (must match your genome assembly)
chrom_sizes <- read.table("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_hg38_p14.chrom.sizes", col.names = c("chr", "size"))
seqinfo_obj <- Seqinfo(seqnames = chrom_sizes$chr, seqlengths = chrom_sizes$size)

# Write each sample to a BigWig file
for (i in seq_along(bismark_methylKit_objs_filt_norm)) {
  
  sample_obj <- bismark_methylKit_objs_filt_norm[[i]]
  sample_df <- getData(sample_obj)
  
  # Calculate beta values
  sample_df$beta <- with(sample_df, numCs / (numCs + numTs))
  
  # Create GRanges
  gr <- GRanges(
    seqnames = sample_df$chr,
    ranges = IRanges(start = sample_df$start, end = sample_df$end),
    strand = "*",
    score = sample_df$beta
  )
  
  # Add Seqinfo
  seqinfo(gr) <- seqinfo_obj[seqlevels(gr)]
  
  # Write to BigWig
  out_path <- paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/methylKit_output/norm_bw_cpm1/",
                     sample_obj@sample.id, ".bw")
  write_bigwig(gr, out_path)
  
  message("✔ Written: ", out_path)
}



