library(tidyverse)
library(plyranges)
library(rtracklayer)
library(data.table)
library(ggplot2)
library(ggpubr)

# 1. read in activated and inactivated L1s (PMM) ========
activeL1_sense <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/activationScore/filtered_L1/sense_L1s.tsv")
activeL1_antisense <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/activationScore/filtered_L1/antisense_L1s.tsv")
all_activeL1 <- rbind(activeL1_sense, activeL1_antisense) %>% 
  dplyr::mutate(direction = str_remove_all(group, "_L1score")) %>% 
  dplyr::mutate(direction = str_remove_all(direction, "..+_"))
all_activeL1$active <- rep("Yes", nrow(all_activeL1))

map <- c(
  "CD19-GSM2171212" = "BMPC4",
  "CD19-GSM2171213" = "BMPC5",
  "CD19-GSM2171214" = "BMPC6",
  "CD19-GSM2171215" = "BMPC7")
all_activeL1$sampleID <- str_replace_all(all_activeL1$sampleID, map)

# Reorder sampleID by total count before plotting
activeL1_all_forPlot <- all_activeL1 %>%
  count(sampleID) %>%
  right_join(all_activeL1, by = "sampleID") %>%
  mutate(sampleID = fct_reorder(sampleID, n))  %>% # or reorder(sampleID, n)
  dplyr::mutate(
    repName_less = ifelse(
      repName %in% c("L1HS", "L1PA2", "L1PA3", "L1PA4", 
                     "L1PA5", "L1PA6", "L1PA7", "L1PA8"), 
      repName, "others"))

# load MM data
clinical_data_addIDs <- read_tsv(paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess", "/clinical_data_BMprimary_addID_addL1_uniqueMap.tsv"))

mm_activeL1 <- clinical_data_addIDs %>% 
  dplyr::select(sampleID = SRR_id, n = L1_num, group = L1_quantile_group) %>% 
  dplyr::mutate(group = paste0("MM (", group, ")"))

# supFig
my_comparisons <- list(
  c("BMPC", "PMM"),
  c("BMPC", "MM (Q1)"),
  c("PMM",  "MM (Q1)"),
  c("MM (Q1)", "MM (Q2)"),
  c("MM (Q1)", "MM (Q3)"),
  c("MM (Q1)", "MM (Q4)")
)

plot_df <- activeL1_all_forPlot %>%
  dplyr::mutate(group = ifelse(str_detect(sampleID, "BMPC"), "BMPC", "PMM")) %>% 
  dplyr::select(sampleID, n, group) %>% 
  unique() %>% 
  rbind(mm_activeL1) %>% 
  dplyr::mutate(
    group = factor(
      group,
      levels = c("BMPC", "MM (Q1)", "MM (Q2)", "MM (Q3)", "MM (Q4)", "PMM")))

write_tsv(plot_df, file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/cibersort/L1_num_cor/mm_pmm_L1_num.tsv")


my_comparisons <- combn(
  levels(plot_df$group),
  2,
  simplify = FALSE
)

activatL1_MMPMM_plot <- activeL1_all_forPlot %>%
  dplyr::mutate(group = ifelse(str_detect(sampleID, "BMPC"), "BMPC", "PMM")) %>% 
  dplyr::select(sampleID, n, group) %>% 
  unique() %>% 
  rbind(mm_activeL1) %>% 
  dplyr::mutate(
    group = factor(
      group,
      levels = c("BMPC", "MM (Q1)", "MM (Q2)", "MM (Q3)", "MM (Q4)", "PMM"))) %>% 
  ggplot(aes(x = group , y = n, fill = group)) +  
  geom_boxplot(outlier.shape = NA)+
  # geom_violin(
  #   trim = FALSE,
  #   width = 1.2,
  #   alpha = 0.45,
  #   color = "grey40"
  # ) +
  # geom_boxplot(
  #   width = 0.08,
  #   outlier.shape = NA,
  #   fill = "white",
  #   color = "black"
  # ) +
  geom_jitter(
    width = 0.12,
    size = 1,
    alpha = 0.65
  ) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     method.args = list(exact = FALSE),
                     # method = "t.test", 
                     p.adjust.method = "BH",
                     # step.increase = 0.07,
                     size = 2.8,
                     label = "p.signif")+
  # Make facet labels UPPERCASE
  guides(fill = "none") +   #2-column legend
  theme_pubr(base_size = 9.6, legend = c(0.8, 0.3)) +  # Set base font size
  theme(
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.key.size = unit(0.4, "lines"),     # make the color square smaller
    legend.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9, face = "bold")  # or just size = 9
  ) +
  labs(x = "", y = "Num. Active")

ggsave(paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/hg38_mappability", "/activatL1_MMPMM_plot_forSUPP.pdf"), 
       plot = activatL1_MMPMM_plot, width = 10, height = 9, units = "cm", 
       limitsize = FALSE)




