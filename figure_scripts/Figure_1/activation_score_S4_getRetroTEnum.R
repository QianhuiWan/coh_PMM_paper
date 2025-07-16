
# heatmap
# plot for all L1s

input_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/activationScore_R2/filtered_retroTE"

# .libPaths("/home/qwan/miniconda3/envs/coh/lib/R/library")
library(tidyverse)
library(data.table)
library(dplyr)
library(readr)
library(rlang)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(patchwork)

options(scipen=999)

# read in tsv files

# read in tsv files
## get path of tsv files
activeRetroTE_path <- list.files(input_dir, recursive = TRUE, 
                                 pattern = paste0(".*\\.tsv$"), 
                                 full.names = TRUE)

activeRetroTE_list <- list()
for (i in 1:length(activeRetroTE_path)) {
  sample <-  str_remove_all(basename(activeRetroTE_path[i]), "_.+")
  # tsv_name_i <- str_replace_all(basename(activeL1_path[i]), "^[^_]+_", "")
  tsv_name_i <- str_remove_all(basename(activeRetroTE_path[i]), ".tsv")
  tsv_file <- read_tsv(file = activeRetroTE_path[i])
  tsv_file$sampleID <- rep(sample, nrow(tsv_file))
  activeRetroTE_list[[tsv_name_i]] <- tsv_file
}

# Function to match classes to the first data frame
match_classes <- function(df, ref_df) {
  for (col in names(ref_df)) {
    if (col %in% names(df)) {
      df[[col]] <- as(df[[col]], class(ref_df[[col]]))
    }
  }
  return(df)
}

# Apply the function to the list
ref_df <- activeRetroTE_list[[1]] # Reference data frame (the first in the list)
activeRetroTE_list <- lapply(activeRetroTE_list, match_classes, ref_df = ref_df)

activeRetroTE_combined_df <- activeRetroTE_list %>% 
  dplyr::bind_rows(.id = "group") 

# get default situation (named sense) L1 activation
activeRetroTE_combined_sense <- activeRetroTE_combined_df %>% 
  dplyr::filter(str_detect(group, "_sense_retroTEscore"))

write_tsv(activeRetroTE_combined_sense, 
          file = paste0(input_dir, "/sense_retroTEs.tsv"))

# get antisense L1 activation
activeRetroTE_combined_antisense <- activeRetroTE_combined_df %>% 
  dplyr::filter(str_detect(group, "_antisense_retroTEscore"))

write_tsv(activeRetroTE_combined_antisense, 
          file = paste0(input_dir, "/antisense_retroTEs.tsv"))




# activeRetroTE_combined_sense <- read_tsv(paste0(input_dir, "/sense_retroTEs.tsv"))
# activeRetroTE_combined_antisense <- read_tsv(paste0(input_dir, "/antisense_retroTEs.tsv"))

activeRetro_all <- rbind(activeRetroTE_combined_sense, activeRetroTE_combined_antisense) %>% 
  dplyr::mutate(direction = str_remove_all(group, "_retroTEscore")) %>% 
  dplyr::mutate(direction = str_remove_all(direction, "..+_"))

# Reorder sampleID by total count before plotting
activeRetro_all_forPlot <- activeRetro_all %>%
  count(sampleID) %>%
  right_join(activeRetro_all, by = "sampleID") %>%
  mutate(sampleID = fct_reorder(sampleID, n)) 

# the n above is L1 number for each sample
# we also want to add L1PA2 and L1PA3 number and save the output
activeRetro_all_forSheet <- activeRetro_all_forPlot %>% 
  dplyr::mutate(active_retroTE = n) %>% 
  dplyr::group_by(repName) %>% 
  count(sampleID, .drop = FALSE)
  # summarise(count = n(), .groups = "drop")  # Summarise counts
  


# rmsk annotation
## rmsk only have been filtered out hg38 patches; keep chr1:22 and chrXYM
rmsk_hg38 <- fread("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/repeatMasker_hg38")
rmsk_hg38_te_ids <- fread("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_hg38_rmsk_teCoord_ids.txt",
                          col.names = c("genoName", "genoStart", "genoEnd", "strand", "te_id"))
## -1 since we did genoEnd +1 in python (because python is [) structure)
rmsk_hg38_te_ids <- rmsk_hg38_te_ids %>% mutate(genoEnd = genoEnd-1)

## get L1 annotations
retroTE_info <- rmsk_hg38_te_ids %>% 
  dplyr::left_join(rmsk_hg38, 
                   by = c("genoName" ="genoName", "genoStart"="genoStart",
                          "genoEnd"="genoEnd", "strand" = "strand")) %>% 
  dplyr::select(genoName, genoStart, genoEnd, strand, 
                te_id, repName, repClass, repFamily) %>% 
  # get LTR, LINE, SINE and Retroposon
  dplyr::filter(repClass %in% c("LTR", "LINE", "SINE", "Retroposon")) %>% 
  dplyr::filter(genoEnd - genoStart > 2000)



# Fisher enrichment test function
## background counts
bg_counts <- retroTE_info %>% 
  count(repName, name = "background")


active_counts <- activeRetro_all_forPlot %>% 
  count(repName, name = "activated")


# merge counts
merged <- full_join(active_counts, bg_counts, by = "repName") %>%
  mutate(across(c(activated, background), ~replace_na(., 0)))

# total te numbers
total_activated <- length(unique(activeRetro_all_forPlot$Geneid))
total_background <- nrow(retroTE_info)

# Fisher test: for each repeat_name
enrichment_results <- merged %>%
  rowwise() %>%
  mutate(
    other_activated = total_activated - activated,
    other_background = total_background - background,
    p_value = fisher.test(matrix(
      c(activated, background, other_activated, other_background),
      nrow = 2
    ), alternative = "greater")$p.value
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    fold_enrichment = (activated / total_activated) / (background / total_background)
  ) %>%
  arrange(p_adj)

# check sig. repeat names
head(enrichment_results)


# plot enrichment for activated retro TEs

p <- enrichment_results %>% 
  dplyr::filter(p_adj<0.05) %>%
  dplyr::arrange(p_adj) %>% 
  dplyr::filter(activated > 1) %>%
  dplyr::filter(background >10) %>%
  slice(1:10) %>% 
  ggplot(aes(x=reorder(repName, p_adj), 
             y=fold_enrichment, 
             colour = p_adj, 
             size = activated))+
  geom_point()+
  geom_hline(yintercept = 1, linetype = "dashed",linewidth = 0.6, col="red")+
  # facet_wrap(~group, nrow = 3)+
  labs(x="Retrotransposon family", y = "Observed/Expected", 
       colour = "FDR", 
       size = "Number of activated TE")+
  guides(
    colour = guide_colourbar(position = "right"),
    size   = guide_legend(position = "top")
  ) +
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.9))



ggsave(paste0(input_dir, "/repeatFam_retroTEactivation_enrich_R1_fdr0.05top10.pdf"),
       plot = p, 
       device = "pdf", width = 20, height = 12, units = "cm")




