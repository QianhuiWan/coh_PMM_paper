
# check L1 length distribution
# output dir: s6_activationScore_R3

# data paths
# dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s5_activationScore_R4"
# dir_gemini = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s5_activationScore_R4_gemini"

output_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_redo_res"

# .libPaths("/home/qwan/miniconda3/envs/coh/lib/R/library")
library(tidyverse)
library(data.table)
library(dplyr)
library(readr)
library(rlang)
library(ggplot2)
library(ggpubr)
options(scipen=999)


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

## plot distribution
L1_length_plot <- 
  L1_info %>% 
  dplyr::mutate(L1_length = genoEnd - genoStart) %>% 
  dplyr::filter(L1_length>4000) %>%
  ggplot(aes(L1_length))+
  # geom_histogram(bins = 200)+
  geom_histogram(aes(y = after_stat(density)), bins =500, fill = "skyblue", alpha = 0.8)+
  geom_density(adjust = 0.5)+
  labs(x = "L1 length", y= "Density")+
  theme_pubr()


## plot distribution in 2 y axis
### Basic counts
hist_data <- L1_info %>%
  mutate(L1_length = genoEnd - genoStart) %>%
  filter(L1_length > 4000)

max_count <- max(ggplot_build(
  ggplot(hist_data, aes(x = L1_length)) +
    geom_histogram(bins = 500)
)$data[[1]]$count)

scale_factor <- max_count * 200

L1_length_plot_2y <- ggplot(hist_data, aes(x = L1_length)) +
  geom_histogram(bins = 500, fill = "skyblue", alpha = 0.8) +
  geom_density(aes(y = after_stat(density * scale_factor)), adjust = 0.5, color = "red") +
  scale_y_continuous(
    name = "Number of L1",
    sec.axis = sec_axis(~./scale_factor, name = "Density")
  ) +
  labs(x = "L1 length") +
  theme_pubr()

ggsave(paste0(output_dir, "/L1_length_distribution.pdf"),
       plot = L1_length_plot_2y,
       device = "pdf", width = 15, height = 10, units = "cm")



