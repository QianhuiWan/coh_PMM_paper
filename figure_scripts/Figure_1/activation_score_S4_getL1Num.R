
# heatmap
# plot for all L1s

input_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/activationScore_R2/filtered_L1"

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

activeL1_combined_sense <- read_tsv(paste0(input_dir, "/sense_L1s.tsv"))

# get predicted antisense L1 activation
activeL1_combined_antisense <- read_tsv(paste0(input_dir, "/antisense_L1s.tsv"))

activeL1_all <- rbind(activeL1_combined_sense, activeL1_combined_antisense) %>% 
  dplyr::mutate(direction = str_remove_all(group, "_L1score")) %>% 
  dplyr::mutate(direction = str_remove_all(direction, "..+_"))

# Reorder sampleID by total count before plotting
activeL1_all_forPlot <- activeL1_all %>%
  count(sampleID) %>%
  right_join(activeL1_all, by = "sampleID") %>%
  mutate(sampleID = fct_reorder(sampleID, n))  %>% # or reorder(sampleID, n)
  dplyr::mutate(
    repName_less = ifelse(
      repName %in% c("L1HS", "L1PA2", "L1PA3", "L1PA4", 
                     "L1PA5", "L1PA6", "L1PA7", "L1PA8"), 
      repName, "others")) 

fig1_plot5 <- activeL1_all_forPlot %>%
  ggplot(aes(x = sampleID, fill = repName_less)) +  
  geom_bar() +
  # Make facet labels UPPERCASE
  # facet_wrap(~direction, labeller = labeller(direction = tools::toTitleCase)) +  
  coord_flip() +
  scale_fill_viridis_d(option = "turbo") +
  guides(fill = guide_legend(ncol = 2)) +   #2-column legend
  theme_pubr(base_size = 9, legend = c(0.8, 0.3)) +  # Set base font size
  theme(
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.key.size = unit(0.4, "lines"),     # make the color square smaller
    legend.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9, face = "bold")  # or just size = 9
  ) +
  labs(x = "", y = "Num. Active", fill = "Repeat names")

ggsave(paste0(input_dir, "/num_active_all_L1s_pmm_barplot_fig1E.pdf"), 
       plot = fig1_plot5, width = 10, height = 5.5, units = "cm", 
       limitsize = FALSE)

# the n above is L1 number for each sample
# we also want to add L1PA2 and L1PA3 number and save the output
activeL1_all_forSheet <- activeL1_all_forPlot %>% 
  dplyr::mutate(L1 = n) %>% 
  dplyr::group_by(repName) %>% 
  count(sampleID, .drop = FALSE)
  # summarise(count = n(), .groups = "drop")  # Summarise counts
  
activeL1PA2_forSheet <- activeL1_all_forSheet[activeL1_all_forSheet$repName=="L1PA2",] %>% 
  dplyr::mutate(L1PA2 = n) %>% 
  as.data.frame() %>% 
  dplyr::select(sampleID, L1PA2)
  
activeL1PA3_forSheet <- activeL1_all_forSheet[activeL1_all_forSheet$repName=="L1PA3",] %>% 
  dplyr::mutate(L1PA3 = n) %>% 
  as.data.frame() %>% 
  dplyr::select(sampleID, L1PA3)

activeL1_all_forSheet_save <- activeL1_all_forPlot %>% 
  as.data.frame() %>% 
  dplyr::mutate(L1 = n) %>% 
  dplyr::select(-n) %>% 
  left_join(activeL1PA2_forSheet, by = "sampleID") %>% 
  left_join(activeL1PA3_forSheet, by = "sampleID")
  
write_tsv(activeL1_all_forSheet_save, 
          file = paste0(input_dir, "/num_active_all_L1s_L1PA2_L1PA3_pmm.tsv"))




