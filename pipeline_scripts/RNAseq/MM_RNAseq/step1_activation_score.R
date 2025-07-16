
# filter L1 with reads num. ratio > random cutoff for each sample
# modified based on step1_R2; we use s5_activationScore_R4_gemini as input (feature count used)
# output dir: s6_activationScore_R3


# data paths
# dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s5_activationScore_R4"
dir_gemini = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s5_activationScore_R4_gemini"

args <- commandArgs(trailingOnly = TRUE)
# sample <- "SRR1606199"
sample <- args[1]
output_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_redo"

# .libPaths("/home/qwan/miniconda3/envs/coh/lib/R/library")
library(tidyverse)
library(data.table)
library(dplyr)
library(readr)
library(rlang)
options(scipen=999)


# read in txt files ############################################################
## Define your initial directory
message("Data input directory: ", dir_gemini)

## Search for files in the directory
score_path <- list.files(dir_gemini, recursive = TRUE,
                         pattern = paste0(sample, ".*\\.txt$"),
                         full.names = TRUE)

## Final check
if (length(score_path) < 4) {
    stop("No .txt files found in directory.")
} else {
    message("Found files: ", paste(score_path, collapse = ", "))
}

## for loop to read in all txt files
sample_txt_list <- list()
for (i in 1:length(score_path)) {
  txt_name_i <- str_remove_all(basename(score_path[i]), 
                               pattern = paste0(sample, "_"))
  txt_name_i <- str_remove_all(txt_name_i, ".txt")
  position <- str_remove_all(txt_name_i, "..+_hg38_")
  position <- str_remove_all(position, "_counts")
  dt <- fread(file = score_path[i])
  colnames(dt) <- c("Geneid", "Chr", "Start", "End", "Strand", 
                    "Length", "counts")
  dt <- dt[, c("Geneid", "Strand", "counts")]
  dt$position <- rep(position, nrow(dt))
  sample_txt_list[[txt_name_i]] <- dt
}

sample_combined_df <- sample_txt_list %>% 
  dplyr::bind_rows(.id = "category") 

# get random region results
sample_combined_random <- sample_combined_df %>% 
  dplyr::filter(str_detect(category, "random")) %>% 
  pivot_wider(id_cols = c(Geneid,Strand) , names_from = position, 
              values_from = counts) %>% 
  # dplyr::filter(!(start100bp==0 & end100bp ==0)) %>%
  dplyr::mutate(Ratio_Sense = 
                  ifelse(Strand == "+", end100bp/pmax(start100bp, 1), 
                         start100bp/pmax(end100bp, 1))) %>%
  dplyr::mutate(Ratio_Antisense_predict = 
                  ifelse(Strand == "+", start100bp/pmax(end100bp, 1), 
                         end100bp/pmax(start100bp, 1)))

# get TE region results
sample_combined_TE <- sample_combined_df %>% 
  dplyr::filter(str_detect(category, "rmsk")) %>% 
  pivot_wider(id_cols = c(Geneid,Strand) , names_from = position, 
              values_from = counts) %>% 
  # dplyr::filter(!(start100bp==0 & end100bp ==0)) %>% 
  dplyr::mutate(Ratio_Sense = 
                  ifelse(Strand == "+", end100bp/pmax(start100bp, 1), 
                         start100bp/pmax(end100bp, 1))) %>%
  dplyr::mutate(Ratio_Antisense_predict = 
                  ifelse(Strand == "+", start100bp/pmax(end100bp, 1), 
                         end100bp/pmax(start100bp, 1)))


# match TE id with rmsk annotation
## rmsk only have been filtered out hg38 patches; keep chr1:22 and chrXYM
rmsk_hg38 <- fread("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/repeatMasker_hg38")
rmsk_hg38_te_ids <- fread("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_hg38_rmsk_teCoord_ids.txt",
                          col.names = c("genoName", "genoStart", "genoEnd", "strand", "te_id"))
## -1 since we did genoEnd +1 in python (because python is [) structure)
rmsk_hg38_te_ids <- rmsk_hg38_te_ids %>% mutate(genoEnd = genoEnd-1)

## get L1 annotations
L1_info <- rmsk_hg38_te_ids %>% 
  dplyr::left_join(rmsk_hg38, 
                   by = c("genoName" ="genoName", "genoStart"="genoStart",
                          "genoEnd"="genoEnd", "strand" = "strand")) %>% 
  dplyr::select(genoName, genoStart, genoEnd, strand, 
                te_id, repName, repClass, repFamily) %>% 
  dplyr::filter(repFamily == "L1") 


# only consider one (i.e. sense) situation #####################################
# random top5 percent cutoff
random_cutoff <- sample_combined_random %>% 
  summarise(top1perc = quantile(Ratio_Sense, probs = 0.99))
message("Random region cutoff (sense): ", as.integer(random_cutoff))

# get L1 region pass the random cutoff
sample_combined_L1_pass <- sample_combined_TE %>% 
  dplyr::filter(Geneid %in% L1_info$te_id) %>% 
  dplyr::left_join(L1_info, by = c("Geneid" = "te_id")) %>% 
  dplyr::filter(Ratio_Sense > random_cutoff$top1perc) %>% 
  dplyr::filter(genoEnd - genoStart > 5900)

# objects to save: save output tsv for a sample
write_tsv(sample_combined_L1_pass, 
          file = paste0(output_dir, "/", sample, "/",
                        sample, "_predict_sense_L1score.tsv"))


# predicted anti-sense situation ###############################################
# random top5 percent cutoff
random_cutoff_predAnti <- sample_combined_random %>% 
  summarise(top1perc = quantile(Ratio_Antisense_predict, probs = 0.99))
message("Random region cutoff (antisense): ", as.integer(random_cutoff_predAnti))

# get L1 region pass the random cutoff
sample_combined_L1_predictAnti_pass <- sample_combined_TE %>% 
  dplyr::filter(Geneid %in% L1_info$te_id) %>% 
  dplyr::left_join(L1_info, by = c("Geneid" = "te_id")) %>% 
  dplyr::filter(Ratio_Antisense_predict > random_cutoff_predAnti$top1perc) %>% 
  dplyr::filter(genoEnd - genoStart > 5900)

# objects to save: save output tsv for a sample
write_tsv(sample_combined_L1_predictAnti_pass, 
          file = paste0(output_dir, "/", sample, "/",
                        sample, "_predict_antisense_L1score.tsv"))











