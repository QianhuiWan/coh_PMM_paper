
# detect KZFP bindiing TFs, step 1, get ZNF141 ZNF382 TSS regions
# TSS regions start with -1000bp upstream and end with the end of 5UTR regions
# S3 : get some plots

args <- commandArgs(trailingOnly = TRUE)
# set cell type variable
cell_type <- "H1"
# cell_type <- "HepG2"
# cell_type <- args[1]

outdir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_star_allo_round4/counts_tx_round4/znfTx_5UTR"

# load pkgs
library(tidyverse)
library(data.table)
library(rlang)
library(rtracklayer)
library(biomaRt)
library(plyranges)
library(ggplot2)
library(ComplexHeatmap)
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg38)
huGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38


# ====== Step 1: Prepare ZNF141 and ZNF382 Regions and other data needed ======

## load gtf annotation files
data <- import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/gencode_v46_chr_patch_hapl_scaff_annotation.gtf")
data <- as.data.frame(data, stringsAsFactors = FALSE)
data <- data[data$type == "transcript",]

deGenes <- readr::read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_star_allo_round4/counts_tx_round4/DE_res.tsv")
dat <- readr::read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_star_allo_round4/counts_tx_round4/topTable_res.tsv")
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
# get SD and mean
bmpc_lcpm_SD_mean <- lcpm_df %>% 
  mutate(lcpm = ifelse(lcpm < 0, 0, lcpm)) %>% 
  dplyr::filter(str_detect(sample, "BMPC")) %>% 
  dplyr::group_by(gene_name) %>% 
  dplyr::summarise(SD_lcpm = sd(lcpm, na.rm =TRUE),
                   mean_lcpm = mean(lcpm, na.rm =TRUE),
                   SD_mean_lcpm = SD_lcpm/mean_lcpm)

pmm_lcpm_SD_mean <- lcpm_df %>% 
  mutate(lcpm = ifelse(lcpm < 0, 0, lcpm)) %>% 
  dplyr::filter(str_detect(sample, "PMM")) %>% 
  dplyr::group_by(gene_name) %>% 
  dplyr::summarise(SD_lcpm = sd(lcpm, na.rm =TRUE),
                   mean_lcpm = mean(lcpm, na.rm =TRUE),
                   SD_mean_lcpm = SD_lcpm/mean_lcpm)


tss1000_ZNFtx_GR <- read_tsv(
  file = paste0(outdir, "/tss1000bp_ZNF141ZNF382tx_rangeReduced.tsv")) %>% 
  as_granges()

# get binding TFs
all_overlaps_ZNF141ZNF382 <- 
  read_tsv(paste0(outdir, "/all_overlaps_ZNF141ZNF382_bindingTFs_", 
                 cell_type,".tsv"))


# filter TFs: keep expressed TFs in PMM and BMPC
all_overlaps_ZNF141ZNF382_fil <- all_overlaps_ZNF141ZNF382 %>%
  dplyr::filter(log10_qVale > 1.3) %>% # FDR<0.05
  dplyr::filter(TF_name %in% lcpm_df_fil$gene_name)
  # dplyr::filter(TF_name %in% dat[dat$FDR < 0.05, ]$gene_name) %>%
  # dplyr::filter(TF_name %in% dat[dat$logFC < -0.5, ]$gene_name) %>%
  # dplyr::filter(TF_name %in% bmpc_lcpm_SD_mean[bmpc_lcpm_SD_mean$mean_lcpm>1,]$gene_name) %>%
  # dplyr::filter(TF_name %in% pmm_lcpm_SD_mean[pmm_lcpm_SD_mean$mean_lcpm>0,]$gene_name) %>%
  # dplyr::filter(TF_name %in% pmm_lcpm_SD_mean[pmm_lcpm_SD_mean$SD_mean_lcpm >=0.1, ]$gene_name)

# ====== Step 4: Summarize Overlaps ======
# Count overlaps per KZFP
ZNF141ZNF382_bind_counts <- all_overlaps_ZNF141ZNF382_fil %>%
  group_by(TF_name) %>%       # Group by KZFP_name
  # filter(n() > 10) %>%
  # ungroup()  %>%
  group_by(gene_name) %>%
  summarise(overlap_count = n_distinct(TF_name), 
            singnalValue = mean(singnalValue)) %>%
  arrange(desc(overlap_count))

# ====== Step 5: Visualize Results ==========

# Plot the most frequently overlapping repeat regions
p1 <- 
  ZNF141ZNF382_bind_counts %>% 
  # ggplot(aes(x = reorder(repName, overlap_count), y = overlap_count)) +
  ggplot(aes(x = reorder(gene_name, singnalValue), y = overlap_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "",
       x = "KZFPs promoters",
       y = "Number of binding TFs") +
  theme_minimal()

ggsave(paste0(outdir, "/ZNF141ZNF382_bind_counts_", cell_type, ".pdf"),
       plot = p1, 
       device = "pdf", width = 10, height = 7, units = "cm")



# plot signalValue heatmap ######==========
mat <- all_overlaps_ZNF141ZNF382_fil %>% 
  # dplyr::mutate(TF_id = paste0(TF_id, "_", TF_name)) %>% 
  dplyr::select(gene_name, singnalValue, TF_id, TF_name) %>% 
  pivot_wider(
    names_from = TF_name, 
    values_from = singnalValue,
    id_cols = gene_name, 
    values_fill = list(singnalValue = 0),  # Fill missing values with 0
    values_fn = list(singnalValue = mean)  # Handle duplicates by averaging
  ) %>% 
  column_to_rownames(var = "gene_name") %>% 
  as.matrix()

scaled_mat = t(scale(t(mat)))

library(circlize)
col_fun = colorRamp2(c(0, 50, max(mat)), c("white", "lightblue",  "darkblue"))

pdf(file = paste0(outdir, "/all_overlaps_ZNF141ZNF382_fil_SVheatmap_", 
                  cell_type, ".pdf"), 
    # width = 8, height = 4)
    width = 20, height = 4)
ComplexHeatmap::Heatmap(mat,name = "Singal value", col = col_fun, # mean SV
                        show_row_dend = FALSE, row_names_side = "left", 
                        column_title_side = "bottom",
                        row_title = "KZFP promoter", 
                        column_title = paste0("Binding TFs (", cell_type, ")"))
dev.off()

# add KZFP Gviz plots for PHF8, KDM4A, KDM1A, HDAC2, ASH2L and EZH2 chipSeq data


# add cor plot between L1PA2 numbers and these gene expression
PMM_metadata <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/PMM_sampleSheet.tsv")
# PMM_tetx <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_star_allo_round4/counts_tx_round3/tetx.tsv")

lcpm_df_fil_addL1num <- lcpm_df_fil %>% 
  dplyr::filter(gene_name %in% c(colnames(mat), rownames(mat), "DNMT1", "DNMT3B")) %>% 
  dplyr::filter(str_detect(sample, "PMM")) %>% 
  dplyr::left_join(PMM_metadata, by = c("sample" = "sampleName"))

p3 <- lcpm_df_fil_addL1num %>% 
  dplyr::filter(gene_name %in% c(
    # "TBP", "KDM4A", 
    "KDM1A", "HDAC2", "EZH2", "DNMT1", "DNMT3B"
    # "ZNF141", "ZNF382"
    )) %>%
  dplyr::mutate(gene_name =
                  factor(gene_name,
                         levels = c("KDM1A", "HDAC2", "EZH2", "DNMT1", "DNMT3B"))) %>%
  # dplyr::filter(gene_name %in% c(c("ZNF687", "E2F4", "E2F8", "KLF6", "KLF4", "MAX")
  #                                # "ZNF141", "ZNF382"
  #                                )) %>%
  # dplyr::mutate(gene_name =
  #                 factor(gene_name,
  #                        levels = c("ZNF687", "E2F4", "E2F8", "KLF6", "MAX"))) %>%
  ggplot(aes(x = lcpm, y = L1PA2)) + 
  geom_point() + geom_smooth(method = "lm")+ 
  facet_wrap(~gene_name, scales = "free_x", nrow = 2)+
  ggpubr::stat_cor(method = "pearson", size=3, label.y = 150) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3)+
  theme_set(ggpubr::theme_pubr(base_size = 12))

ggsave(paste0(outdir, "/key_TFs_L1PA2_cor_", cell_type, "_v2.pdf"),
       plot = p3, 
       device = "pdf", 
       width = 23, height = 10,
       # width = 22, height = 35, 
       units = "cm")

# check gene expression
# lcpm_df %>%
lcpm_df_fil %>%
  # dplyr::mutate(gene_id = str_replace_all(gene_name, "_.+", ""),
  #               gene_name = str_replace_all(gene_name, ".+_", "")) %>%
  # dplyr::filter(str_detect(gene_name, "KLF")) %>%
  # dplyr::filter(gene_name %in% c("ZNF687", "E2F4", "E2F8", "MAX", "KLF6")) %>%
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
  dplyr::filter(gene_name %in% c("KDM4A", "KDM1A",
                                 "SAP18", "SAP30", "ARID4A","ARID4B", "PHF8",
                                 "ASH2L", "HDAC2", "EZH2", "HDAC1", "MAX",
                                 "CHD3", "CHD4", "MBD2", "MBD3",
                                 "GATAD2A", "GATAD2B", "RBBP4", "RBBP7",
                                 "E2F4", "RBL2", "E2F1", "E2F8", "EGR1", "E2F5",
                                 "MYBL2", "FOXM1", "ZNF564", "JUND", "SOX6",
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

ggsave(paste0(outdir, "/test_", cell_type, ".pdf"),
       device = "pdf", 
       # width = 12, height = 10, 
       width = 22, height = 9, 
       units = "cm")


