
# in this round3, change to methrix R packge
library(tidyverse)
library(magrittr)
library(methrix)
library(stringr)
library(plyranges)
library(rtracklayer)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MafDb.1Kgenomes.phase3.GRCh38)

# read bed.graph files
# files <- list.files("/Users/qianhuiwan/gitHubRepo/coh_PMM/data/testRes_DNAm/reformatedCpGmeth", pattern = "formatted_perc.txt", full.names = TRUE) 
files <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/*/03.DNMtoolsPosMethylation/*_CpG.meth.gz")

print(basename(files))

outDir <- '/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1'

# get sample id and group
sampleID <- str_replace_all(basename(files), "_CpG.meth.gz", "")
group <- str_replace_all(sampleID, "[:digit:]", "")
# group <- ifelse(group == "B_rest", 0, ifelse(group == "PMM", 1, 2))

# meth_snps_fil_mat <- readRDS(paste0(outDir, "/meth_snps_fil_matGR_withXY.rds"))
meth_snps_fil_mat <- readRDS("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1/meth_snps_cov_fil_matGR_withXY.rds")

PMM_sampleSheet <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/gene_counts/PMM_sampleSheet_R3_addBMPC.tsv")
PMM_sampleSheet <- PMM_sampleSheet[PMM_sampleSheet$group != "CD19+BMPC",] %>%
  dplyr::filter(str_detect(sampleName, "PMM"))
PMM_sampleSheet_DNAm <- read_tsv("/home/qwan/githubRepo/coh_PMM/sampleSheet/RNA/PMM_metadata_supFile.tsv") %>% 
  dplyr::select(sampleName = ID, Infilteration) %>% 
  dplyr::left_join(PMM_sampleSheet[,c("sampleName", "L1", "L1PA2", "L1PA3")], by = "sampleName")

# check infilteration cor ######################################################
# v_obj <- voom(dge, design)
# E <- v_obj$E  # voom 的 logCPM matrix (genes x samples)
# pc <- prcomp(t(E), scale. = TRUE)
library(bsseq)
bsseq_obj <- readRDS("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/methylKit_output/bismark_methylKit_bsseq_obj.rds")
meth_snps_fil_mat2 <- getMeth(bsseq_obj, type = "raw")

meth_snps_fil_mat_t <- meth_snps_fil_mat2 %>% as.data.frame() %>% 
  dplyr::mutate(cpg_id = paste0("cpg", 1:nrow(.))) %>% 
  dplyr::select(cpg_id, starts_with("PMM")) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames(var = "cpg_id")

# meth_snps_fil_mat_t <- meth_snps_fil_mat %>% as.data.frame() %>% 
#   dplyr::mutate(cpg_id = paste0(seqnames, "_", start, "_", end)) %>% 
#   dplyr::select(cpg_id, starts_with("PMM")) %>% 
#   tibble::column_to_rownames(var = "cpg_id")

mat <- as.matrix(meth_snps_fil_mat_t)  # df = 你 column_to_rownames 后的矩阵
row_sds <- apply(mat, 1, sd, na.rm = TRUE)    # 每行标准差
keep <- is.finite(row_sds) & row_sds > 0
mat_filt <- mat[keep, , drop = FALSE]

min_frac <- 0.8
keep2 <- rowMeans(!is.na(mat_filt)) >= min_frac
mat_filt <- mat_filt[keep2, , drop = FALSE]

mat_filt[is.nan(mat_filt)] <- NA
mat_filt[is.infinite(mat_filt)] <- NA

mat_noNA <- mat_filt[ rowSums(is.na(mat_filt)) == 0, , drop = FALSE ]
dim(mat_noNA)

pc <- prcomp(t(mat_noNA), scale. = TRUE)

saveRDS(pc, file = paste0(outDir, "/meth_snps_fil_mat_PCs_v2.rds"))

# pc <- readRDS(file = paste0(outDir, "/meth_snps_fil_mat_PCs_v2.rds"))

df_pc <- data.frame(PC1=pc$x[,1], PC2=pc$x[,2]) %>% 
  tibble::rownames_to_column(var = "sampleName") %>% 
  dplyr::left_join(PMM_sampleSheet_DNAm, by = "sampleName") 

# cor.test(PMM_sampleSheet$Infilteration, normalized_lcpm[str_detect(rownames(normalized_lcpm), pattern = "ZNF141"),])

ggplot(df_pc, aes(PC1, PC2)) +
  geom_point(size=3) +
  theme_classic() +
  labs(title="PCA (DNAm)")

df_pmm <- df_pc
df_pmm$sample <- df_pmm$sampleName
df_pmm$infiltration <- df_pmm$Infilteration

# cor.test(df_pmm$PC1, df_pmm$infiltration, method = "spearman")
# cor.test(df_pmm$PC2, df_pmm$infiltration, method = "spearman")
# 
# ct1 <- cor.test(df_pmm$PC1, df_pmm$infiltration, method="spearman", exact=FALSE)
# ct2 <- cor.test(df_pmm$PC2, df_pmm$infiltration, method="spearman", exact=FALSE)
# 
# lab1 <- sprintf("Spearman ρ = %.3f\np = %.3g", ct1$estimate, ct1$p.value)
# lab2 <- sprintf("Spearman ρ = %.3f\np = %.3g", ct2$estimate, ct2$p.value)

# Pearson correlation tests
cor.test(df_pmm$PC1, df_pmm$Infilteration, method = "pearson")
cor.test(df_pmm$PC2, df_pmm$Infilteration, method = "pearson")

ct1 <- cor.test(df_pmm$PC1, df_pmm$infiltration, method = "pearson")
ct2 <- cor.test(df_pmm$PC2, df_pmm$infiltration, method = "pearson")

lab1 <- sprintf("Pearson r = %.3f\np = %.3g", ct1$estimate, ct1$p.value)
lab2 <- sprintf("Pearson r = %.3f\np = %.3g", ct2$estimate, ct2$p.value)


p1 <- ggplot(df_pmm, aes(infiltration, PC1)) +
  geom_point(size=3) +
  geom_smooth(method="lm", se=FALSE) +
  geom_text_repel(aes(label=sample), size=3, max.overlaps=Inf) +
  annotate("text",
           x = Inf, y = Inf, label = lab1,
           hjust = 1.1, vjust = 1.1, size = 3.5) +
  theme_classic() +
  labs(title="PMM DNAm: PC1 vs Infiltration")

p2 <- ggplot(df_pmm, aes(infiltration, PC2)) +
  geom_point(size=3) +
  geom_smooth(method="lm", se=FALSE) +
  geom_text_repel(aes(label=sample), size=3, max.overlaps=Inf) +
  annotate("text",
           x = Inf, y = Inf, label = lab2,
           hjust = 1.1, vjust = 2, size = 3.5) +
  theme_classic() +
  labs(title="PMM DNAm: PC2 vs Infiltration")

# pw <- p1 / p2 / p3
pw <- p1 / p2 
pw

## 4) save patchwork to pdf (vector)
ggsave(
  filename = paste0(outDir, "/pmm_DNAm_pc_infiltrationCor_v2.pdf"),
  plot = pw,
  width = 5,
  height = 6)







