
# R2: use cov filtered input + k=101 default; use this output in paper
# set the library path:
# .libPaths("/home/qwan/miniconda3/envs/coh/lib/R/library")
# .libPaths("/home/qwan/R/Rstudio-WithModules-WithJupyter.4.3.3.sif/libs")
library(tidyverse)
library(magrittr)
library(bsseq)
library(methrix)
library(plyranges)
library(BiocParallel)
library(mixOmics)
library(gridExtra)


outDir <- '/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts'
KZFP_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs"
# get KZFP ages
KZFP_age <- read_tsv(paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/KZFP_evolutional_age.tsv"))

# get PMM and MM L1 binding KZFPs expression

## PMM
### load sample sheet with L1 subfam numbers
PMM_ss <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/PMM_sampleSheet_R2.tsv")
PMM_ss_L1PAnum <- PMM_ss %>% 
  dplyr::select(sample = sampleName, L1PA2, L1PA3)

### load lcpm PMM
# lcpm_df <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_star_allo_round4/counts_tx_round4/lcpm_TMM.tsv")
lcpm_df <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/lcpm_TMM_PMM_STARcounts.tsv")

L1binding_KZFPs <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/all_overlaps_L1_binding_KZFPs_PMM_R4.tsv")
# L1binding_KZFPs <- L1binding_KZFPs %>%
#   group_by(KZFP_name) %>%       # Group by KZFP_name
#   filter(n() >= 10) %>%
#   ungroup()
L1PA2binding_KZFPs <- L1binding_KZFPs %>% dplyr::filter(repName == "L1PA2")
L1PA2binding_KZFPs_fil <- L1PA2binding_KZFPs %>%
  dplyr::filter(KZFP_name %in% c("ZNF425", "ZNF141", "ZNF382","ZNF680", "ZNF28")) %>% 
  group_by(KZFP_name) %>%  
  summarise(num_bindingSites = n()) %>% 
  dplyr::mutate(KZFP_name =str_replace_all(KZFP_name, "ZNF28", "ZFP28")) %>% 
  dplyr::select(KZFP = KZFP_name, num_bindingSites)
  

KZFP_lcpm_PMM <- lcpm_df %>% 
  dplyr::mutate(gene_name_original = gene_name) %>% 
  dplyr::mutate(gene_name = str_replace_all(gene_name_original, "_.+", "")) %>% 
  dplyr::filter(gene_name %in% c("ZNF425", "ZNF141", "ZNF382","ZNF680", "ZFP28")) %>% 
  dplyr::select(sample, KZFP = gene_name, KZFP_lcpm =lcpm) %>% 
  dplyr::left_join(KZFP_age, by = "KZFP") %>% 
  dplyr::left_join(L1PA2binding_KZFPs_fil, by = "KZFP") %>% 
  mutate(across(c(KZFP_age, num_bindingSites), ~ (. - min(.)) / (max(.) - min(.)), 
                .names = "scaled_{.col}")) %>% 
  dplyr::mutate(repressionScore = (KZFP_lcpm*scaled_num_bindingSites)) %>% 
  group_by(sample) %>% 
  summarise(repressionScore_mean = sum(repressionScore)/5) %>% 
  dplyr::left_join(PMM_ss_L1PAnum, by = "sample") %>% 
  dplyr::mutate(disease = rep("PMM", nrow(.)))

# plot
KZFP_lcpm_PMM %>% 
ggplot(aes(x = L1PA2, y = repressionScore_mean)) + 
  geom_point() + geom_smooth(method = "lm")+ 
  # facet_wrap(~KZFP, scales = "free_x")+
  # ylim(c(-20,150))+
  ggpubr::stat_cor(label.y=c(0.9,0.9), method = "pearson",size=3) + 
  ggrepel::geom_text_repel(aes(label = sample), max.overlaps = 100,
                           box.padding = 0.5, size = 3)+
  labs(y= "Repression score")+
  ggpubr::theme_pubr(base_size = 9.6)


## MM
MM_lcpm <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR1/lcpm_TMM_not_fil.tsv")

MM_ss <- read_tsv(paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess", "/clinical_data_BMprimary_addID.tsv")) #764
MM_ss_L1PAnum <- MM_ss %>% 
  dplyr::select(sample = sampleID, L1PA2)

L1binding_KZFPs <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/all_overlaps_L1_binding_KZFPs_MM_R4.tsv")
L1binding_KZFPs <- L1binding_KZFPs %>%
  group_by(KZFP_name) %>%       # Group by KZFP_name
  filter(n() > 10) %>%
  ungroup()
L1PA2binding_KZFPs <- L1binding_KZFPs %>% dplyr::filter(repName == "L1PA2")

L1PA2binding_KZFPs_fil <- L1PA2binding_KZFPs %>%
  dplyr::filter(KZFP_name %in% c("ZNF324", "ZNF317", "ZNF136","ZNF425", 
                                 "ZNF765", "ZNF671", "ZNF17")) %>% 
  group_by(KZFP_name) %>%  
  summarise(num_bindingSites = n()) %>% 
  dplyr::select(KZFP = KZFP_name, num_bindingSites)

# plot
KZFP_lcpm_MM <- MM_lcpm %>% 
  dplyr::mutate(gene_name_original = gene_name) %>% 
  dplyr::mutate(gene_name = str_replace_all(gene_name_original, "_.+", "")) %>% 
  dplyr::filter(gene_name %in% c("ZNF324", "ZNF317", "ZNF136","ZNF425", 
                                 "ZNF765", "ZNF671", "ZNF17"))
KZFP_lcpm_MM <- KZFP_lcpm_MM %>% 
  dplyr::mutate(gene_name_original = gene_name) %>% 
  dplyr::mutate(gene_name = str_replace_all(gene_name_original, "_.+", "")) %>% 
  dplyr::filter(gene_name %in% c("ZNF324", "ZNF317", "ZNF136","ZNF425", 
                                 "ZNF765", "ZNF671", "ZNF17")) %>% 
  dplyr::select(sample, KZFP = gene_name, KZFP_lcpm =lcpm) %>% 
  dplyr::left_join(KZFP_age, by = "KZFP") %>% 
  dplyr::left_join(L1PA2binding_KZFPs_fil, by = "KZFP") %>% 
  mutate(across(c(KZFP_age, num_bindingSites), ~ (. - min(.)) / (max(.) - min(.)), 
                .names = "scaled_{.col}")) %>% 
  dplyr::mutate(repressionScore = (KZFP_lcpm*scaled_num_bindingSites)/(scaled_KZFP_age+1)) %>% 
  group_by(sample) %>% 
  summarise(repressionScore_mean = sum(repressionScore)/7) %>% 
  dplyr::left_join(MM_ss_L1PAnum, by = "sample") %>% 
  dplyr::mutate(disease = rep("MM", nrow(.)))

p <- KZFP_lcpm_PMM %>% 
  dplyr::select(sample, repressionScore_mean, L1PA2, disease) %>% 
  rbind(KZFP_lcpm_MM) %>% 
  ggplot(aes(x = L1PA2, y = repressionScore_mean)) + 
  geom_point() + geom_smooth(method = "lm")+ 
  # facet_wrap(~KZFP, scales = "free_x")+
  # ylim(c(-20,150))+
  ggpubr::stat_cor(label.y=c(0.9,0.9), method = "pearson",size=3) + 
  labs(y= "Repression score")+
  facet_wrap(~factor(disease,levels = c("PMM", "MM")), 
             scales = "free_x",ncol = 1)+
  ggpubr::theme_pubr(base_size = 9.6)

ggsave(paste0(outDir, "/repression_score_MM_PMM.pdf"), 
       plot = p,
       width = 8, height = 15, units = "cm", limitsize = FALSE)












