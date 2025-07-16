

library(tidyverse)
library(magrittr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# PMM
outdir_PMM <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_star_allo_round4/counts_tx_round4/znfTx_5UTR"
PMM_metadata <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/PMM_sampleSheet.tsv")


lcpm_df <- readr::read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_star_allo_round4/counts_tx_round4/lcpm_TMM.tsv")


# long to wide 
lcpm_df_wide <- lcpm_df %>% 
  mutate(lcpm = ifelse(lcpm < 0, 0, lcpm)) %>% 
  pivot_wider(names_from = sample, values_from = lcpm,
              id_cols = gene_name)
# filter gene no 0 (lcpm>=1) >=3
lcpm_df_fil <- lcpm_df_wide %>% 
  rowwise() %>%
  filter(sum(c_across(-gene_name) > 1) >= 3) %>% 
  pivot_longer(cols = c(starts_with("BMPC"), starts_with("PMM")), 
               names_to = "sample", 
               values_to = "lcpm")

lcpm_df_fil %>%
  # dplyr::mutate(gene_id = str_replace_all(gene_name, "_.+", ""),
  #               gene_name = str_replace_all(gene_name, ".+_", "")) %>%
  # dplyr::filter(str_detect(gene_name, "KLF")) %>%
  # dplyr::filter(gene_name %in% c("ZNF687", "E2F4", "E2F8", "MAX", "KLF6")) %>%
  # dplyr::filter(gene_name %in% c("CGAS", "AIM2", "IFI16", "DDX41",
  #                                "STING1", "ERIS", 
  #                                "MAVS", "PYCARD", "IFIH1", "DDX58",
  #                                "TBK1", "IKBKE", "IRF3", "IRF7",
  #                                "NFKB1", "RELA", "IFNB1", "CXCL10", "TNF", 
  #                                "ISG15", "MX1",
  #                                "IL6", "IL1B")) %>%
  # dplyr::filter(gene_name %in% c("HSP60", "DDIT3", "ATF4", "XBP1",
  #                                "TRAF3","TRAF6","IRF1",
  #                                "TGFB1","FOXP3", "ARG1", "VEGFA", "VEGFB",
  #                                "CSF1R", "CD163", "IL10",
  #                                "HDAC1","HDAC2", "DNMT1","DNMT3B")) %>% 
  # dplyr::filter(gene_name %in% c("AIM2", "ZBP1", "TREX1",
  #                                "ISG15", "MX1", "OAS1", "IRF7",
  #                                "TP53BP1","H2AX", "RAD51", "BRIP1", "NBN", "XRCC2",
  #                                "XRCC6", "XRCC5", "LIG4", "SLFN11", "PARG",
  #                                "ATM", "ATR", "CHEK1", "CHEK2", "BRCA1", "BRCA2",
  #                                "CGAS", "TBK1", "IRF3", "STAT1", "STAT2", "STAT3",
  #                                "SOCS1", "USP18" , "FANCD2", "PALB2")) %>% 
  # dplyr::filter(gene_name %in% c("TASOR", "TASOR2")) %>% 
  dplyr::filter(gene_name %in% c("KDM4A", "KDM1A",
                                 "SAP18", "SAP30", "ARID4A","ARID4B", "PHF8",
                                 "ASH2L", "HDAC2", "EZH2", "HDAC1", "MAX",
                                 "CHD3", "CHD4", "MBD2", "MBD3",
                                 "GATAD2A", "GATAD2B", "RBBP4", "RBBP7",
                                 "E2F4", "RBL2", "E2F1", "E2F8", "EGR1", "E2F5",
                                 "MYBL2", "FOXM1", "FOXO3A", "ZNF564", "JUND", "SOX6",
                                 "LIN54", "KLF6")) %>%
  distinct() %>%
  ggplot(aes(x=sample, y=lcpm, fill=gene_name))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(title = "Histone genes", x="")+
  facet_wrap(~gene_name, scales = "free_y", 
             # nrow = 1
  )+
  # facet_wrap(~gene_name+transcript_id, scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "none")

ggsave(paste0(outdir_PMM, "/test_", cell_type, ".pdf"),
       device = "pdf", 
       # width = 12, height = 10, 
       width = 22, height = 9, 
       units = "cm")


# PMM: lcpm from same pipe as MM
lcpm_PMM <- readr::read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/lcpm_TMM_PMM_STARcounts.tsv")

# long to wide 
lcpm_PMM_wide <- lcpm_PMM %>% 
  dplyr::select(gene_name = gene_name, sample, lcpm) %>% 
  pivot_wider(names_from = sample, values_from = lcpm,
              id_cols = gene_name)
# filter gene no 0 (lcpm>=1) >=3
lcpm_PMM_fil <- lcpm_PMM_wide %>% 
  rowwise() %>%
  filter(sum(c_across(-gene_name) > 1) >= 3) %>% 
  pivot_longer(cols = c(starts_with("BMPC"), starts_with("PMM")), 
               names_to = "sample", 
               values_to = "lcpm")

lcpm_PMM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>% 
  dplyr::left_join(PMM_metadata, by = c("sample" = "sampleName")) %>% 
  # dplyr::filter(gene_name %in% "NRF1") %>%
  # dplyr::filter(gene_name %in% c("IGSF9", "DMAP1", "CKLF", "AP4B1", "CTBP2",
  #                                "POLE4", "GPR155", "RPS3AP7", "RPSAP11",
  #                                "SH3BP5", "GPR155", "RPL17P18", "ENSG00000279070",
  #                                "R3HDM1", "INPP5B")) %>%
  dplyr::filter(gene_name %in% c("ZNF248", "ZNF33A", "ZNF425", "ZNF519", "ZNF93", "ZNF649")) %>%
  # dplyr::filter(str_detect(gene_name, "KLF")) %>%
  # dplyr::filter(gene_name %in% c("ZNF687", "E2F4", "E2F8", "MAX", "KLF6","CEBPA")) %>%
  # dplyr::filter(gene_name %in% c("KLF6","KLF4")) %>%
  # dplyr::filter(gene_name %in% c("SOCS1","CD274", "CD163","TP53", "MDM2", "MDM4", "ZEB2",
  #                                "BBC3","PMAIP1", "CDKN1A", "CDKN2A", "CDK6",
  #                                "E2F1", "MYC", "RB1", "ARF1",
  #                                "HRAS", "NRAS", "KRAS", "STING1", "IFI16",
  #                                "ENPP1", "TREX1", "USP18","SOCS1",
  #                                "MAP1LC3B")) %>%
  # cell cycle
  # dplyr::filter(gene_name %in% c("CDK4", "CDK6", "CCND1", "CCNE1")) %>%
  # dplyr::filter(gene_name %in% c("CD274", "CD163", "MRC1", "MYC", "MDM2",
  #                                "H2AX", "CDCA7", "CDCA7L", "NFKB1","NFKBIE",
  #                                "RIGI", "IFIH1", "DHX58", "SP140", "USP18",
  #                                "ISG15", "IL10", "TGFB1", "SOCS1", "SOCS3",
  #                                "NFKBIA", "PDCD1", "LAG3", "HAVCR2", "CCL2",
  #                                "LGALS9", "CEACAM1", "HMGB1","ANXA5",
  #                                "TNFAIP3", "ARG1", "CBX5", "SUV39H1",
  #                                "UHRF1", "OAS1", "POU5F1")) %>%
  # dplyr::filter(gene_name %in% c("TLR3", "TLR7", "TLR9",
  #                                "DDX58", "IFIH1", 
  #                                "STING1", "CGAS", 
  #                                "MYD88","TRIF", "MAVS", "TBK1", "JAK1", 
  #                                "IRF3", "TRF7", 
  #                                "ATM", "ATR", "CHEK2", "STAT1", "NFKB1", "SP1",
  #                                "PML"
  #                                )) %>%
  # dplyr::filter(str_detect(gene_name, "TLR")) %>%
  # dplyr::filter(str_detect(gene_name, "TCF")) %>%
  # dplyr::filter(str_detect(gene_name, "SOCS")) %>%
  # dplyr::filter(gene_name %in% c("CGAS", "AIM2", "IFI16", "DDX41",
  #                                "STING1", "ERIS",
  #                                "MAVS", "PYCARD", "IFIH1", "DDX58",
  #                                "TBK1", "IKBKE", "IRF3", "IRF7",
  #                                "NFKB1", "RELA", "IFNB1", "CXCL10", "TNF",
  #                                "ISG15", "MX1",
  #                                "IL6", "IL1B")) %>%
  # dplyr::filter(gene_name %in% c("HSP60", "DDIT3", "ATF4", "XBP1",
  #                                "TRAF3","TRAF6","IRF1",
  #                                "TGFB1","FOXP3", "ARG1", "VEGFA", "VEGFB",
  #                                "CSF1R", "CD163", "IL10",
  #                                "HDAC1","HDAC2", "DNMT1","DNMT3B")) %>%
  # dplyr::filter(gene_name %in% c("ZBP1", "TREX1",
  #                                "ISG15", "MX1", "OAS1", "IRF7",
  #                                "TP53BP1","H2AX", "RAD51", "BRIP1", "NBN", "XRCC2",
  #                                "XRCC6", "XRCC5", "LIG4", "SLFN11", "PARG",
  #                                "ATM", "ATR", "CHEK1", "CHEK2", "BRCA1", "BRCA2",
  #                                "CGAS", "TBK1", "IRF3", "STAT1", "STAT2", "STAT3",
  #                                "SOCS1", "USP18" , "FANCD2", "PALB2")) %>%
  # dplyr::filter(gene_name %in% c("TASOR", "TASOR2")) %>%
  # dplyr::filter(gene_name %in% c("KDM4A", "KDM1A", "TBP", "DNMT1","DNMT3B",
  #                                "SAP18", "SAP30", "ARID4A","ARID4B", "PHF8",
  #                                "ASH2L", "HDAC2", "EZH2", "HDAC1", "MAX",
  #                                "CHD3", "CHD4", "MBD2", "MBD3",
  #                                "GATAD2A", "GATAD2B", "RBBP4", "RBBP7",
  #                                "E2F4", "RBL2", "E2F1", "E2F8", "EGR1", "E2F5",
  #                                "MYBL2", "FOXM1", "ZNF564", "JUND", "SOX6",
  #                                "LIN54", "KLF6")) %>%
  distinct() %>%
  ggplot(aes(x=reorder(sample, L1PA2), y=lcpm, fill=gene_name))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(title = "Genes in PMM", x="")+
  facet_wrap(~gene_name, scales = "free_y", 
             # nrow = 1
  )+
  # facet_wrap(~gene_name+transcript_id, scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "none")

ggsave(paste0(outdir_PMM, "/test_PMM_KLFs_expr_v2.pdf"),
       device = "pdf",
       width = 12, height = 10,
       # width = 20, height = 20,
       units = "cm")


# sig HR KZFPS
sig_KZFPS <- read_tsv(paste0(out_dir, "/gene_HR_MMRF_sig_KZFPs.tsv")) %>% 
  dplyr::select(-group)

sig_L1_KZFPS <- read_tsv(paste0(out_dir, "/gene_HR_MMRF_sig_L1_KZFPS.tsv"))


lcpm_PMM_fil_KZFP <- 
  lcpm_PMM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>% 
  dplyr::left_join(PMM_metadata, by = c("sample" = "sampleName")) %>% 
  dplyr::filter(gene_name %in% sig_L1_KZFPS$gene_symbol) %>% 
  distinct() %>%
  dplyr::left_join(sig_L1_KZFPS, by = c("gene_name"="gene_symbol"))
  
lcpm_PMM_fil_KZFP_mean <- lcpm_PMM_fil_KZFP %>% 
  group_by(sample, Prognosis, primate_specific_KZFPs) %>% 
  summarise(mean_lcpm = mean(lcpm), .groups = 'drop') %>% 
  left_join(PMM_metadata, by = c("sample"="sampleName"))

p3_PMM <- lcpm_PMM_fil_KZFP %>% 
  group_by(sample, Prognosis) %>% 
  summarise(mean_lcpm = mean(lcpm), .groups = 'drop') %>% 
  left_join(PMM_metadata, by = c("sample"="sampleName")) %>% 
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm))+
  ggplot(aes(x=group, y=mean_lcpm))+
  geom_boxplot(outlier.shape = NA, alpha = 0.9, position = position_dodge(0.8)) +
  geom_jitter(position = position_jitter(width = 0.15), alpha = 0.5) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif"
  ) +
  labs(x ="", y="Mean expression (lcpm)")+
  facet_wrap(~Prognosis) +
  # geom_text_repel(aes(label = sample), 
  #                 position = position_jitter(width = 0.15))+
  geom_text_repel(aes(label = sample, color = Prognosis),
                  size = 3, max.overlaps = 20, segment.color = "gray25")+
  scale_color_manual(
    values = c(
      "Poor prognosis" = "darkred",  # warm soft red
      "Good prognosis" = "#457B9D"  # dusty blue
    )) +
  guides(color = "none")+
  theme_pubr()

ggsave(paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/allgene_hazard", 
              "/gene_HR_sig_L1binding_KZFPS_expre_PMMpatients.pdf"), 
       plot = p3_PMM,
       width = 12, height = 10, units = "cm", limitsize = FALSE)


p3_2_PMM <- lcpm_PMM_fil_KZFP_mean %>% 
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm))+
  ggplot(aes(x=group, y=mean_lcpm, fill = primate_specific_KZFPs))+
  geom_boxplot(outlier.shape = NA, alpha = 0.9, position = position_dodge(0.8)) +
  geom_jitter(position = position_jitter(width = 0.15), alpha = 0.5) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif"
  ) +
  labs(x ="", y="Mean expression (lcpm)")+
  facet_wrap(~Prognosis) +
  # geom_text_repel(aes(label = sample), 
  #                 position = position_jitter(width = 0.15))+
  geom_text_repel(aes(label = sample, color = primate_specific_KZFPs),
                  size = 3, max.overlaps = 15, segment.color = "gray25")+
  guides()+
  theme_pubr()

ggsave(paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/allgene_hazard", 
              "/gene_HR_sig_L1_KZFPS_expre_OP_PMMpatients.pdf"), 
       plot = p3_2_PMM,
       width = 12, height = 10, units = "cm", limitsize = FALSE)



lcpm_PMM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>% 
  dplyr::filter(gene_name %in% c("KLF6", "KLF4")) %>%
  dplyr::left_join(PMM_metadata, by = c("sample" = "sampleName")) %>% 
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm, fill=subtype, group = 1))+
  ggplot(aes(x=group, y=lcpm))+
  # geom_bar(stat = "identity", position = "dodge")+
  # geom_col() +
  # geom_line(color = "red", size = 1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, position = position_dodge(0.8)) +
  geom_jitter(position = position_jitter(width = 0.15), alpha = 0.5) +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    label = "p.signif"
  ) +
  labs(x ="", y="Gene expression (lcpm)")+
  facet_wrap(~gene_name) +
  ggpubr::theme_pubr(base_size = 9.6)

ggsave(paste0(outdir_PMM, "/PMM_KLF4KLF6_expr.pdf"),
       device = "pdf", 
       # width = 12, height = 10, 
       width = 10, height = 8, 
       units = "cm")


lcpm_PMM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>% 
  dplyr::filter(gene_name %in% c("NRF1")) %>%
  dplyr::left_join(PMM_metadata, by = c("sample" = "sampleName")) %>% 
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm, fill=subtype, group = 1))+
  ggplot(aes(x=group, y=lcpm))+
  # geom_bar(stat = "identity", position = "dodge")+
  # geom_col() +
  # geom_line(color = "red", size = 1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, position = position_dodge(0.8)) +
  geom_jitter(position = position_jitter(width = 0.15), alpha = 0.5) +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    label = "p.signif"
  ) +
  labs(x ="", y="Gene expression (lcpm)")+
  facet_wrap(~gene_name) +
  geom_text_repel(aes(label = sample),
                  size = 3, max.overlaps = 20, segment.color = "gray25")+
  ggpubr::theme_pubr(base_size = 9.6)

ggsave(paste0(outdir_PMM, "/PMM_NRF1_expr.pdf"),
       device = "pdf", 
       # width = 12, height = 10, 
       width = 10, height = 8, 
       units = "cm")


# IGSF9
# calculate repeat transcript number
PMM_metadata <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/PMM_sampleSheet.tsv")
L1_tetx_num <- read_tsv(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/activationScore_R2/filtered_L1/num_active_all_L1s_L1PA2_L1PA3_pmm.tsv")

L1_TEtx_num <- L1_tetx_num %>%
  dplyr::select(sample = sampleID, L1_fc = L1,
                L1PA2_fc = L1PA2, L1PA3_fc=L1PA3) %>% unique()
L1_TEtx_num
PMM_metadata <- PMM_metadata %>%
  left_join(L1_TEtx_num, by = c("sampleName"="sample"))


igsf9_cor_PMM <- lcpm_PMM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>%
  dplyr::left_join(PMM_metadata, by = c("sample" = "sampleName")) %>%
  dplyr::filter(gene_name %in% c("IGSF9")) %>%
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm))+
  ggplot(aes(x=lcpm, y=L1PA2_fc))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson", label.y=c(120,1), size=3)+
  labs(x ="Gene expression (lcpm)", y="L1PA2")+
  geom_text_repel(aes(label = sample),
                  size = 3, max.overlaps = 20, segment.color = "gray25")+
  guides(color = "none")+
  theme_pubr()
igsf9_cor_PMM
igsf9_cor_PMM <- lcpm_PMM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>%
  dplyr::left_join(PMM_metadata, by = c("sample" = "sampleName")) %>%
  dplyr::filter(gene_name %in% c("IGSF9")) %>%
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm))+
  ggplot(aes(x=lcpm, y=L1PA2_fc))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson", label.y=c(120,1), size=3)+
  labs(x ="Gene expression (lcpm)", y="L1PA2")+
  geom_text_repel(aes(label = sample),
                  size = 3, max.overlaps = 20, segment.color = "gray25")+
  guides(color = "none")+
  theme_pubr()
ggsave(paste0(output_dir,
              "/IGSF9_L1PA2_cor.pdf"),
       plot = igsf9_cor_PMM,
       width = 12, height = 10, units = "cm", limitsize = FALSE)
igsf9_cor_PMM <- lcpm_PMM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>%
  dplyr::left_join(PMM_metadata, by = c("sample" = "sampleName")) %>%
  dplyr::filter(gene_name %in% c("IGSF9")) %>%
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm))+
  ggplot(aes(x=lcpm, y=L1PA2_fc))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson", label.y=c(120,1), size=3)+
  labs(x ="IGSF9 expression (lcpm)", y="L1PA2")+
  geom_text_repel(aes(label = sample),
                  size = 3, max.overlaps = 20, segment.color = "gray25")+
  guides(color = "none")+
  theme_pubr()
ggsave(paste0(output_dir,
              "/IGSF9_L1PA2_cor.pdf"),
       plot = igsf9_cor_PMM,
       width = 12, height = 10, units = "cm", limitsize = FALSE)

igsf9_ZNF_cor_PMM <- lcpm_PMM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>%
  dplyr::left_join(PMM_metadata, by = c("sample" = "sampleName")) %>%
  dplyr::filter(gene_name %in% c("IGSF9", "ZNF141", "ZNF382")) %>%
  pivot_wider(id_cols = c(sample, group, L1PA2, L1PA3, L1PA2_fc),
              names_from = gene_name,
              values_from = lcpm ) %>%
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm))+
  ggplot(aes(x=IGSF9, y=ZNF382))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson", label.y=c(0.9,1), size=3)+
  labs(x ="Gene expression (lcpm)", y="")+
  geom_text_repel(aes(label = sample),
                  size = 3, max.overlaps = 20, segment.color = "gray25")+
  guides(color = "none")+
  theme_pubr()
igsf9_ZNF_cor_PMM
igsf9_ZNF_cor_PMM
igsf9_ZNF_cor_PMM <- lcpm_PMM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>%
  dplyr::left_join(PMM_metadata, by = c("sample" = "sampleName")) %>%
  dplyr::filter(gene_name %in% c("IGSF9", "ZNF141", "ZNF382")) %>%
  pivot_wider(id_cols = c(sample, group, L1PA2, L1PA3, L1PA2_fc),
              names_from = gene_name,
              values_from = lcpm ) %>%
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm))+
  ggplot(aes(x=IGSF9, y=ZNF382))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson", label.y=c(3.2,1), size=3)+
  labs(x ="IGSF9 expression (lcpm)", y="ZNF382 expression (lcpm)")+
  geom_text_repel(aes(label = sample),
                  size = 3, max.overlaps = 20, segment.color = "gray25")+
  guides(color = "none")+
  theme_pubr()
igsf9_ZNF_cor_PMM
igsf9_ZNF141_cor_PMM <- lcpm_PMM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>%
  dplyr::left_join(PMM_metadata, by = c("sample" = "sampleName")) %>%
  dplyr::filter(gene_name %in% c("IGSF9", "ZNF141", "ZNF382")) %>%
  pivot_wider(id_cols = c(sample, group, L1PA2, L1PA3, L1PA2_fc),
              names_from = gene_name,
              values_from = lcpm ) %>%
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm))+
  ggplot(aes(x=IGSF9, y=ZNF141))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson", label.y=c(3.2,1), size=3)+
  labs(x ="IGSF9 expression (lcpm)", y="ZNF382 expression (lcpm)")+
  geom_text_repel(aes(label = sample),
                  size = 3, max.overlaps = 20, segment.color = "gray25")+
  guides(color = "none")+
  theme_pubr()
igsf9_ZNF141_cor_PMM
igsf9_ZNF141_cor_PMM <- lcpm_PMM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>%
  dplyr::left_join(PMM_metadata, by = c("sample" = "sampleName")) %>%
  dplyr::filter(gene_name %in% c("IGSF9", "ZNF141", "ZNF382")) %>%
  pivot_wider(id_cols = c(sample, group, L1PA2, L1PA3, L1PA2_fc),
              names_from = gene_name,
              values_from = lcpm ) %>%
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm))+
  ggplot(aes(x=IGSF9, y=ZNF141))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson", label.y=c(6.2,1), size=3)+
  labs(x ="IGSF9 expression (lcpm)", y="ZNF141 expression (lcpm)")+
  geom_text_repel(aes(label = sample),
                  size = 3, max.overlaps = 20, segment.color = "gray25")+
  guides(color = "none")+
  theme_pubr()

igsf9_ZNF382_cor_PMM <- lcpm_PMM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>%
  dplyr::left_join(PMM_metadata, by = c("sample" = "sampleName")) %>%
  dplyr::filter(gene_name %in% c("IGSF9", "ZNF141", "ZNF382")) %>%
  pivot_wider(id_cols = c(sample, group, L1PA2, L1PA3, L1PA2_fc),
              names_from = gene_name,
              values_from = lcpm ) %>%
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm))+
  ggplot(aes(x=IGSF9, y=ZNF382))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson", label.y=c(3.2,1), size=3)+
  labs(x ="IGSF9 expression (lcpm)", y="ZNF382 expression (lcpm)")+
  geom_text_repel(aes(label = sample),
                  size = 3, max.overlaps = 20, segment.color = "gray25")+
  guides(color = "none")+
  theme_pubr()
ggsave(paste0(output_dir,
              "/IGSF9_ZNF382_cor.pdf"),
       plot = igsf9_ZNF382_cor_PMM,
       width = 12, height = 10, units = "cm", limitsize = FALSE)

igsf9_ZNF141_cor_PMM <- lcpm_PMM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>%
  dplyr::left_join(PMM_metadata, by = c("sample" = "sampleName")) %>%
  dplyr::filter(gene_name %in% c("IGSF9", "ZNF141", "ZNF382")) %>%
  pivot_wider(id_cols = c(sample, group, L1PA2, L1PA3, L1PA2_fc),
              names_from = gene_name,
              values_from = lcpm ) %>%
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm))+
  ggplot(aes(x=IGSF9, y=ZNF141))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson", label.y=c(6.2,1), size=3)+
  labs(x ="IGSF9 expression (lcpm)", y="ZNF141 expression (lcpm)")+
  geom_text_repel(aes(label = sample),
                  size = 3, max.overlaps = 20, segment.color = "gray25")+
  guides(color = "none")+
  theme_pubr()

ggsave(paste0(output_dir,
              "/IGSF9_ZNF141_cor.pdf"),
       plot = igsf9_ZNF141_cor_PMM,
       width = 12, height = 10, units = "cm", limitsize = FALSE)

Figure5_B <- ggarrange(igsf9_ZNF141_cor_PMM, igsf9_ZNF382_cor_PMM,
                       igsf9_cor_PMM,
                       labels = c("", "", ""),
                       font.label = list(size = 9.6, face = "bold",
                                         color ="black"),
                       nrow = 3, ncol=1,
                       vjust = c(0, 0),hjust = c(0, 1),
                       common.legend = FALSE)

Figure5_B <- ggarrange(igsf9_ZNF141_cor_PMM, igsf9_ZNF382_cor_PMM,
                       igsf9_cor_PMM,
                       labels = c("", "", ""),
                       font.label = list(size = 9.6, face = "bold",
                                         color ="black"),
                       nrow = 3, ncol=1,
                       vjust = c(0, 0),hjust = c(0, 1),
                       common.legend = FALSE)
ggsave(paste0(output_dir,
              "/Figure5_B.pdf"),
       plot = Figure5_B,
       width = 8, height = 10, units = "cm", limitsize = FALSE)
8*0.2
20*0.2
10*0.2
ggsave(paste0(output_dir,
              "/Figure5_B.pdf"),
       plot = Figure5_B,
       width = 9.6, height = 12, units = "cm", limitsize = FALSE)
Figure5_B <- ggarrange(igsf9_ZNF141_cor_PMM, igsf9_ZNF382_cor_PMM,
                       igsf9_cor_PMM,
                       labels = c("", "", ""),
                       font.label = list(size = 9.6, face = "bold",
                                         color ="black"),
                       nrow = 3, ncol=1,
                       common.legend = FALSE)
Figure5_B
Figure5_B <- ggarrange(igsf9_ZNF141_cor_PMM, igsf9_ZNF382_cor_PMM,
                       igsf9_cor_PMM,
                       labels = c("", "", ""),
                       font.label = list(size = 9.6, face = "bold",
                                         color ="black"),
                       nrow = 3, ncol=1, align = "hv",
                       common.legend = FALSE)
Figure5_B
ggsave(paste0(output_dir,
              "/Figure5_B.pdf"),
       plot = Figure5_B,
       width = 16, height = 20, units = "cm", limitsize = FALSE)
Figure5_B <- ggarrange(igsf9_ZNF141_cor_PMM, igsf9_ZNF382_cor_PMM,
                       igsf9_cor_PMM,
                       labels = c("", "", ""),
                       font.label = list(size = 9.6, face = "bold",
                                         color ="black"),
                       nrow = 3, ncol=1, align = "hv",
                       common.legend = FALSE)
ggsave(paste0(output_dir,
              "/Figure5_B.pdf"),
       plot = Figure5_B,
       width = 16, height = 20, units = "cm", limitsize = FALSE)




Figure5_B <- ggarrange(igsf9_ZNF141_cor_PMM, igsf9_ZNF382_cor_PMM,
                       igsf9_cor_PMM,
                       labels = c("", "", ""),
                       font.label = list(size = 9.6, face = "bold",
                                         color ="black"),
                       nrow = 3, ncol=1,
                       vjust = c(0, 0),hjust = c(0, 1),
                       common.legend = FALSE)


# MM
outdir_MM = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR1"
MM_meta <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/clinical_data_BMprimary_addID.tsv")
MM_meta_update <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/mmrf_subtypes/clinical_data_BMprimary_addID_addSubtypes.tsv")
PR_subtypeIDs <- MM_meta_update[MM_meta_update$RNA_Subtype_Name=="PR",]$sampleID

lcpm <- read.table(paste0(outdir_MM, "/lcpm_TMM_fil.tsv"), 
                   sep = "\t", header = TRUE, row.names = 1) %>%
  as.data.frame() 

lcpm_MM <- lcpm %>% 
  rownames_to_column(var="gene_name") %>% 
  dplyr::mutate(gene_id_version = str_replace_all(gene_name, "_.+", "")) %>% 
  dplyr::mutate(gene_id = str_remove(gene_id_version, "\\..*")) %>% 
  dplyr::mutate(gene_symbol = str_replace_all(gene_name, ".+_", "")) %>% 
  dplyr::filter(!(gene_id == "ENSG00000228314")) %>%  # redundant Pseudogene id
  pivot_longer(cols = -c(gene_name, gene_id_version, gene_id, gene_symbol), 
               names_to = "sampleID", values_to = "lcpm") %>% 
  dplyr::mutate(lcpm = ifelse(lcpm < 0, 0, lcpm))

# long to wide 
lcpm_MM_wide <- lcpm_MM %>% 
  dplyr::select(gene_name = gene_name, sample = sampleID, lcpm) %>% 
  pivot_wider(names_from = sample, values_from = lcpm,
              id_cols = gene_name) 
# filter gene no 0 (lcpm>=1) >=3
lcpm_MM_fil <- lcpm_MM_wide %>% 
  # rowwise() %>%
  # filter(sum(c_across(-gene_name) > 1) >= 3) %>% 
  pivot_longer(cols = c(starts_with("SRR")), 
               names_to = "sample", 
               values_to = "lcpm")


lcpm_MM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, "_.+", ""),
                gene_name = str_replace_all(gene_name, ".+_", "")) %>%
  dplyr::left_join(MM_meta,by = c("sample" = "SRR_id")) %>% 
  # dplyr::filter(gene_name %in% c("NRF1")) %>%
  dplyr::filter(gene_name %in% c("IGSF9")) %>%
  # dplyr::filter(str_detect(gene_name, "KLF")) %>%
  # dplyr::filter(gene_name %in% c("ZNF687", "E2F4", "E2F8", "MAX", "KLF6")) %>%
  # dplyr::filter(gene_name %in% c("CGAS", "AIM2", "IFI16", "DDX41",
  #                                "STING1", "ERIS",
  #                                "MAVS", "PYCARD", "IFIH1", "DDX58",
  #                                "TBK1", "IKBKE", "IRF3", "IRF7",
  #                                "NFKB1", "RELA", "IFNB1", "CXCL10", "TNF",
  #                                "ISG15", "MX1",
  #                                "IL6", "IL1B")) %>%
  # dplyr::filter(gene_name %in% c("HSP60", "DDIT3", "ATF4", "XBP1",
  #                                "TRAF3","TRAF6","IRF1",
  #                                "TGFB1","FOXP3", "ARG1", "VEGFA", "VEGFB",
  #                                "CSF1R", "CD163", "IL10",
  #                                "HDAC1","HDAC2", "DNMT1","DNMT3B")) %>%
  # dplyr::filter(gene_name %in% c("ZBP1", "TREX1",
  #                                "ISG15", "MX1", "OAS1", "IRF7",
  #                                "TP53BP1","H2AX", "RAD51", "BRIP1", "NBN", "XRCC2",
  #                                "XRCC6", "XRCC5", "LIG4", "SLFN11", "PARG",
  #                                "ATM", "ATR", "CHEK1", "CHEK2", "BRCA1", "BRCA2",
  #                                "CGAS", "TBK1", "IRF3", "STAT1", "STAT2", "STAT3",
  #                                "SOCS1", "USP18" , "FANCD2", "PALB2")) %>%
  # dplyr::filter(gene_name %in% c("TASOR", "TASOR2")) %>%
  # dplyr::filter(gene_name %in% c("KDM4A", "KDM1A", "TBP",
  #                                "SAP18", "SAP30", "ARID4A","ARID4B", "PHF8",
  #                                "ASH2L", "HDAC2", "EZH2", "HDAC1", "MAX",
  #                                "CHD3", "CHD4", "MBD2", "MBD3",
  #                                "GATAD2A", "GATAD2B", "RBBP4", "RBBP7",
  #                                "E2F4", "RBL2", "E2F1", "E2F8", "EGR1", "E2F5",
  #                                "MYBL2", "FOXM1", "ZNF564", "JUND", "SOX6",
  #                                "LIN54", "KLF6")) %>%
  distinct() %>%
  ggplot(aes(x=reorder(sample, L1PA2), y=lcpm, fill=gene_name))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(title = "Genes in MM", x="")+
  facet_wrap(~gene_name, scales = "free_y", 
             # nrow = 1
  )+
  # facet_wrap(~gene_name+transcript_id, scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "none")

lcpm_MM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, "_.+", ""),
                gene_name = str_replace_all(gene_name, ".+_", "")) %>%
  # dplyr::filter(gene_name %in% c("KLF6", "KLF4")) %>%
  dplyr::filter(gene_name %in% c("IGSF9")) %>%
  dplyr::left_join(MM_meta,by = c("sample" = "sampleID")) %>% 
  dplyr::mutate(subtype = ifelse(sample %in% na.omit(PR_subtypeIDs), "PR", 
                                 "Others")) %>% 
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm, fill=subtype, group = 1))+
  ggplot(aes(x=subtype, y=lcpm))+
  # geom_bar(stat = "identity", position = "dodge")+
  # geom_col() +
  # geom_line(color = "red", size = 1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, position = position_dodge(0.8)) +
  geom_jitter(position = position_jitter(width = 0.15), alpha = 0.5) +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    label = "p.signif"
  ) +
  labs(x ="", y="Gene expression (lcpm)")+
  facet_wrap(~gene_name) +
  ggpubr::theme_pubr(base_size = 9.6)
  
ggsave(paste0(outdir_MM, "/MM_KLF4KLF6_expr.pdf"),
       device = "pdf", 
       # width = 12, height = 10, 
       width = 10, height = 8, 
       units = "cm")


lcpm_MM_fil %>%
  dplyr::mutate(gene_id = str_replace_all(gene_name, "_.+", ""),
                gene_name = str_replace_all(gene_name, ".+_", "")) %>%
  dplyr::filter(gene_name %in% c("NRF1")) %>%
  dplyr::left_join(MM_meta,by = c("sample" = "sampleID")) %>% 
  dplyr::mutate(subtype = ifelse(sample %in% na.omit(PR_subtypeIDs), "PR", 
                                 "Others")) %>% 
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm, fill=subtype, group = 1))+
  ggplot(aes(x=subtype, y=lcpm))+
  # geom_bar(stat = "identity", position = "dodge")+
  # geom_col() +
  # geom_line(color = "red", size = 1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, position = position_dodge(0.8)) +
  geom_jitter(position = position_jitter(width = 0.15), alpha = 0.5) +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    label = "p.signif"
  ) +
  labs(x ="", y="Gene expression (lcpm)")+
  facet_wrap(~gene_name) +
  ggpubr::theme_pubr(base_size = 9.6)

ggsave(paste0(outdir_MM, "/MM_NRF1_expr.pdf"),
       device = "pdf", 
       # width = 12, height = 10, 
       width = 10, height = 8, 
       units = "cm")





