# R2: featureCount L1 number and epigenetic regulators used in this round

# .libPaths()
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# use count matrix from feature counts
# get gene counts stats, PMM
dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap"
dir_addBMPC <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/"
output_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/gene_counts"

# load packages
library(tidyverse)
library(magrittr)
library(stringr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(rtracklayer)
library(ggpubr)
library(edgeR)
library(limma)
library(missMethyl)
library(msigdbr)
library(patchwork)

# read in Epigenetic regulators (ER genes) and save the gene list back
ER_genes <- read_excel("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/pub_data/epigenetic_regulators/KEPI_A_2139067_SM6382.xlsx",
  sheet = 2, skip = 2
) %>%
  `colnames<-`(c("gene_id", "gene_symbol", "class"))

# write_tsv(ER_genes,
#           file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/pub_data/epigenetic_regulators/ER_genes.tsv")
KZFPs <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/kzfps_hg38_all.tsv")

# load the PMM metadata
PMM_metadata_supFile <- read_tsv(file = "/home/qwan/githubRepo/coh_PMM/sampleSheet/RNA/PMM_metadata_supFile.tsv")


# getNum function
getNum <- function(path, pattern) {
  # get path and sample names
  counts_path <- list.files(path,
    recursive = TRUE, pattern = pattern,
    full.names = TRUE
  )
  # sample_names <- str_remove_all(basename(counts_path), pattern = "_..+")
  sample_names <- str_remove_all(basename(counts_path),
    pattern = "_ReadsPerGene.out.tab"
  )
  # get ref. genome version
  refGenome <- "hg38"

  # load data
  df_list <- list()
  for (i in 1:length(sample_names)) {
    if (str_detect(counts_path[i], "TEs.bed")) {
      # load TEs.bed including FPKM and TPM
      data <- read.table(counts_path[i],
        header = FALSE, sep = "\t",
        stringsAsFactors = FALSE
      )
      colnames(data) <- c(
        "seqnames", "start", "end",
        "gene_id", "transcript_id", "strand", "repeat",
        "cov", "FPKM", "TPM"
      )
      data$cov <- as.numeric(data$cov)
      data$FPKM <- as.numeric(data$FPKM)
      data$TPM <- as.numeric(data$TPM)
    } else if (str_detect(counts_path[i], "ballgown.gtf")) {
      # load tx.gtf including FPKM and TPM
      data <- import(counts_path[i])
      data <- as.data.frame(data, stringsAsFactors = FALSE)
      data <- data[data$type == "transcript", ]
      data$cov <- as.numeric(data$cov)
      data$FPKM <- as.numeric(data$FPKM)
      data$TPM <- as.numeric(data$TPM)
    } else if (str_detect(counts_path[i], "alignedReads.txt")) {
      # load aligned reads number
      data <- read.table(file = counts_path[i], stringsAsFactors = FALSE)
      colnames(data) <- "num_aligned_reads"
    } else if (str_detect(counts_path[i], "gene_counts.txt")) {
      # load aligned reads number
      data <- read.table(
        file = counts_path[i], stringsAsFactors = FALSE,
        header = TRUE
      )
      # colnames(data) <- c("Geneid", "Chr","Start", "End", "Strand",
      #                     "Length", "crick", "watson")
      # data <- data %>% dplyr::mutate(totalCount = crick + watson)
      colnames(data) <- c(
        "Geneid", "Chr", "Start", "End", "Strand",
        "Length", "totalCount"
      )
    } else if (str_detect(counts_path[i], "_ReadsPerGene.out.tab")) {
      # load aligned reads number
      data <- read.table(
        file = counts_path[i], stringsAsFactors = FALSE,
        header = FALSE, skip = 4
      )
      colnames(data) <- c("Geneid", "totalCount", "watson", "crick")
    } else if (str_detect(counts_path[i], "transcript_count_matrix.csv")) {
      # load csv files
      data <- read_csv(counts_path[i])
    } else if (str_detect(counts_path[i], "gene_count_matrix.csv")) {
      data <- read_csv(counts_path[i])
    } else if (str_detect(counts_path[i], "TEtxs_count_matrix.csv")) {
      data <- read.table(file = counts_path[i], stringsAsFactors = FALSE)
      colnames(data) <- c("transcript_id", "counts", "repeat")
    }

    # add data df to df_list
    df_list[[i]] <- as.data.frame(data)
  }

  # name df_list
  names(df_list) <- sample_names

  # retun the df_list
  return(df_list)
}


## use getNum FUN to get all transcripts #######################################

df_list_11 <- getNum(
  path = paste0(dir),
  pattern = "_ReadsPerGene.out.tab"
)

df_list_addBMPC <- getNum(
  path = paste0(dir_addBMPC),
  pattern = "_ReadsPerGene.out.tab"
)

df_list <- c(df_list_11, df_list_addBMPC)

# saveRDS(df_list, file =  paste0(output_dir, "/allSample_genes_list_STARcounts.rds"))
# df_list <- readRDS(file =  paste0(output_dir, "/allSample_genes_list_STARcounts.rds"))

# bind list for gene numbers
## only keep total count > 0 and plot the library size for rnaseq samples
g_count_list <- lapply(df_list, function(x) {
  nrow(x[x[, 2] > 0, ])
})
num_totalG <- do.call("rbind", g_count_list)
num_totalG <- as.data.frame(num_totalG, stringsAsFactors = FALSE)
colnames(num_totalG) <- "num_G"

num_totalG %>%
  rownames_to_column(var = "samples") %>%
  ggplot(aes(x = samples, y = num_G)) +
  geom_bar(stat = "identity") +
  labs(title = "detected genes", x = "", y = "number of total genes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# bind list for gene matrix
## gene list to df and add id column
gene_df <- df_list %>%
  bind_rows(.id = "sampleID")
## save all gene counts to tsv
readr::write_tsv(gene_df,
  file = paste0(output_dir, "/genes_starCounts_PMM_R3_addBMPC.tsv")
)

# load gtf annotation files for human genes
data <- rtracklayer::import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/teAnno_round3/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf")
data <- as.data.frame(data, stringsAsFactors = FALSE)
data <- data[data$type == "transcript", ]

# get unique gene_id and gene_name
genes_anno <- data %>%
  dplyr::filter(gene_id %in% unique(gene_df$Geneid)) %>%
  dplyr::select(gene_id, gene_name) %>%
  unique()


# gene counts analysis
count_matrix <- gene_df %>%
  dplyr::left_join(genes_anno, by = c("Geneid" = "gene_id")) %>%
  # dplyr::mutate(Geneid = paste0(gene_name, "_", Geneid)) %>%
  dplyr::mutate(Geneid = paste0(Geneid, "_", gene_name)) %>%
  # dplyr::select(-watson, -crick) %>% #unique() %>%
  pivot_wider(
    names_from = sampleID, values_from = totalCount,
    id_cols = Geneid, values_fill = 0
  ) %>%
  column_to_rownames(var = "Geneid") %>%
  as.matrix()
# count_matrix[is.na(count_matrix)] <- 0

# Save as a .tsv file
write.table(count_matrix,
  file = paste0(output_dir, "/gene_count_matrix_STARcounts_PMM_R3_addBMPC.tsv"),
  sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE
)

# get sample meta data
## load L1 number info
L1_tetx_num <- read_tsv(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/activationScore/filtered_L1/num_active_all_L1s_L1PA2_L1PA3_pmm.tsv")
L1_tetx_num <- L1_tetx_num %>%
  dplyr::select(sampleID, L1, L1PA2, L1PA3) %>%
  unique()

L1_tetx_num_addBMPC <- read_tsv(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_add_bmpc_R1/filtered_L1/num_active_all_L1s_L1PA2_L1PA3_pmm.tsv")
L1_tetx_num_addBMPC <- L1_tetx_num_addBMPC %>%
  dplyr::select(sampleID, L1, L1PA2, L1PA3) %>%
  unique()

L1_tetx_num <- rbind(L1_tetx_num, L1_tetx_num_addBMPC)

L1_tetx_num %>% write_tsv(paste0(output_dir, "/L1_tetx_num_addBMPC.tsv"))

PMM_sampleSheet <-
  colnames(count_matrix) %>%
  as.data.frame() %>%
  `colnames<-`(c("sampleName")) %>%
  dplyr::mutate(group = ifelse(str_detect(sampleName, "PMM|BMPC"),
    str_replace_all(sampleName, "[:digit:]", ""),
    ifelse(str_detect(sampleName, "CD19-"), "CD19-BMPC", "CD19+BMPC")
  )) %>%
  dplyr::mutate(group2 = ifelse(group == "PMM", "PMM", "BMPC")) %>%
  left_join(L1_tetx_num, by = c("sampleName" = "sampleID")) %>%
  dplyr::mutate(batch = ifelse(group == "CD19-BMPC", "single_end", "pair_end")) %>%
  left_join(PMM_metadata_supFile[, c("ID", "Infilteration")], by = c("sampleName" = "ID")) %>%
  dplyr::mutate(Infilteration = ifelse(is.na(Infilteration), 1, Infilteration))

PMM_sampleSheet %>% write_tsv(paste0(output_dir, "/PMM_sampleSheet_R3_addBMPC.tsv"))


# normalize counts with TMM method
# count_matrix <- read.table(paste0(output_dir, "/gene_count_matrix_STARcounts_PMM_R3_addBMPC.tsv"), check.names = FALSE)
# PMM_sampleSheet <- read_tsv(paste0(output_dir, "/PMM_sampleSheet_R3_addBMPC.tsv"))
table(colnames(count_matrix) == PMM_sampleSheet$sampleName)
## only use CD19- BMPC
PMM_sampleSheet <- PMM_sampleSheet[PMM_sampleSheet$group != "CD19+BMPC", ]
count_matrix <- count_matrix[, colnames(count_matrix) %in% PMM_sampleSheet$sampleName]
table(colnames(count_matrix) == PMM_sampleSheet$sampleName)

dge <- DGEList(counts = count_matrix, group = PMM_sampleSheet$group2)

# Define a threshold for filtering
keep <- rowSums(cpm(dge, normalized.lib.sizes = TRUE) >= 1) >= 8
# keep <- filterByExpr(dge, group=dge$samples$group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
# normalized_lcpm <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)
normalized_lcpm <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)

# cor.test(PMM_sampleSheet$Infilteration, normalized_lcpm[str_detect(rownames(normalized_lcpm), pattern = "ZNF141"),])
# correct for batch after normalization
normalized_lcpm_bc <- removeBatchEffect(
  normalized_lcpm,
  batch  = PMM_sampleSheet$batch, # single/paired end lib batch
  design = model.matrix(~ 0 + PMM_sampleSheet$group2) # keep PMM vs BMPC variations
)

# cor.test(PMM_sampleSheet$Infilteration, normalized_lcpm_bc[str_detect(rownames(normalized_lcpm_bc), pattern = "ZNF141"),])

# stopifnot(all(PMM_sampleSheet$sampleName == colnames(normalized_lcpm)))
# normalized_lcpm_bc <- removeBatchEffect(
#   normalized_lcpm,
#   batch = PMM_sampleSheet$batch,
#   covariates = PMM_sampleSheet$Infilteration,
#   design = model.matrix(~0 + PMM_sampleSheet$group2)
# )


lcpm_df <- normalized_lcpm_bc %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_name") %>%
  pivot_longer(
    cols = c(
      starts_with("BMPC"), starts_with("PMM"),
      starts_with("CD19")
    ),
    names_to = "sample",
    values_to = "lcpm"
  ) # %>%
# mutate(lcpm = ifelse(lcpm < 0, 0, lcpm))

readr::write_tsv(lcpm_df, file = paste0(output_dir, "/lcpm_TMM_PMM_STARcounts_R3_addBMPC_bc.tsv"))
# readr::write_tsv(lcpm_df, file = paste0(output_dir,"/lcpm_TMM_PMM_STARcounts_R3_addBMPC_bc_infilterAdj.tsv"))

# DE
## design matrix, use group as batch
# design <- model.matrix(~0 + group2 + batch + Infilteration, data=PMM_sampleSheet)
design <- model.matrix(~ 0 + group2 + batch, data = PMM_sampleSheet)
colnames(design) <- make.names(colnames(design))
colnames(design)
## make contrast
con <- makeContrasts(
  group = group2PMM - group2BMPC,
  levels = design
)
v <-
  voom(dge, design) %>%
  lmFit(design) %>%
  contrasts.fit(con) %>%
  eBayes()
dt <- decideTests(v)
print(summary(dt))
plotSA(v, main = "Final model: Mean-variance trend")

tfit <- treat(v, lfc = 1)
dt_tfit <- decideTests(tfit)
summary(dt_tfit)

## get topTable output
dat <- topTable(v, coef = "group", number = Inf) %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  dplyr::mutate(gene_id_version = str_replace_all(gene_name, "_.+", "")) %>%
  dplyr::mutate(gene_id = str_remove(gene_id_version, "\\..*")) %>%
  dplyr::mutate(gene_symbol = str_replace_all(gene_name, ".+_", "")) %>%
  # distinct(gene_name, .keep_all = TRUE) %>%
  dplyr::select(-t, -B, FDR = `adj.P.Val`) %>%
  mutate(Sig = ifelse(FDR < 0.05 & logFC > 2, "Up",
    ifelse(FDR < 0.05 & logFC < -2, "Down", "NS")
  ))
### save DE results to tsv file:
write_tsv(dat, file = paste0(output_dir, "/topTable_res_addBMPC_bc.tsv"))
# write_tsv(dat, file = paste0(output_dir,"/topTable_res_addBMPC_bc_infilterAdj.tsv"))

dat_de <- dat %>%
  dplyr::select(
    gene_symbol, gene_id, logFC, AveExpr, P.Value,
    FDR, Sig
  ) %>%
  dplyr::filter(Sig != "NS")
write_tsv(dat_de, file = paste0(output_dir, "/deTable_res_addBMPC_bc.tsv"))

# ER in DE?
dat_DE_ERs <- dat %>%
  dplyr::filter(gene_id %in% ER_genes$gene_id) %>%
  dplyr::filter(Sig != "NS") %>%
  dplyr::left_join(ER_genes, by = "gene_id")

dat_DE_ERs %>%
  ggplot(aes(x = class, fill = Sig)) +
  geom_bar() +
  # facet_wrap(~Sig)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

dat_DE_ERs %>%
  dplyr::filter(gene_symbol.x %in% c(
    "ZNF141", "MPHOSPH8",
    "ATF7IP", "SIRT1",
    "TASOR", "MORC2", "PPHLN1", "SETDB1",
    "TASOR2", "IRF2",
    "DNMT1", "DNMT3A"
  ))
# dat_DE_ERs filter output: 0 rows

# check RIGI/cGAS genes in DE
dat_DE <- dat %>%
  dplyr::filter(Sig != "NS")
## check RIGI genes
dat_DE %>%
  dplyr::filter(gene_symbol %in% c(
    "RIGI", "IFIH1", "DHX58",
    "MAVS", "TBK1",
    "CHUK", "IKBKB",
    "IRF3", "NFKB1", "RELA"
  ))
## check RIGI genes
dat_DE %>%
  dplyr::filter(gene_symbol %in% c(
    "CGAS", "STING1",
    "TBK1",
    "CHUK", "IKBKB",
    "IRF3", "NFKB1", "RELA"
  ))
# dat_DE filter output: only NFKB1

## diff var analysis
# DiffVar
fitvar <- varFit(dge, design = design, coef = c(1, 2)) # coef 1 and 2
fitvar2 <- contrasts.varFit(fitvar, contrasts = con)
summary(decideTests(fitvar2))

topDV <- topVar(fitvar2, coef = 1L, number = nrow(dge$counts)) %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  dplyr::mutate(gene_id_version = str_replace_all(gene_name, "_.+", "")) %>%
  dplyr::mutate(gene_id = str_remove(gene_id_version, "\\..*")) %>%
  dplyr::mutate(gene_symbol = str_replace_all(gene_name, ".+_", "")) %>%
  # distinct(gene_name, .keep_all = TRUE) %>%
  dplyr::select(-t, FDR = `Adj.P.Value`) %>%
  mutate(Sig = ifelse(FDR < 0.05 & LogVarRatio > 2, "Up",
    ifelse(FDR < 0.05 & LogVarRatio < -2, "Down", "NS")
  ))

head(topDV)

### save DE results to tsv file:
write_tsv(topDV, file = paste0(output_dir, "/topTable_DVres_addBMPC_bc.tsv"))
# write_tsv(topDV, file = paste0(output_dir,"/topTable_DVres_addBMPC_bc_infilterAdj.tsv"))

topDV_de <- topDV %>%
  dplyr::select(
    gene_symbol, gene_id, LogVarRatio, SampleVar, P.Value,
    FDR, Sig
  ) %>%
  dplyr::filter(Sig != "NS")
write_tsv(topDV_de, file = paste0(output_dir, "/dvTable_DVres_addBMPC_bc.tsv"))

# ER in DV?
topDV_ERs <- topDV %>%
  dplyr::filter(gene_id %in% ER_genes$gene_id) %>%
  dplyr::filter(Sig != "NS") %>%
  dplyr::left_join(ER_genes, by = "gene_id")

topDV_ERs %>%
  ggplot(aes(x = class, fill = Sig)) +
  geom_bar() +
  # facet_wrap(~Sig)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

topDV_ERs %>%
  dplyr::filter(gene_symbol.x %in% c(
    "ZNF141", "MPHOSPH8",
    "ATF7IP", "SIRT1",
    "TASOR", "MORC2", "PPHLN1", "SETDB1",
    "TASOR2", "IRF2",
    "DNMT1", "DNMT3A"
  ))
# topDV_ERs filter output 0 rows

# check RIGI/cGAS genes in DV
topDV_significant <- topDV %>%
  dplyr::filter(Sig != "NS")
## check RIGI genes
topDV_significant %>%
  dplyr::filter(gene_symbol %in% c(
    "RIGI", "IFIH1", "DHX58",
    "MAVS", "TBK1",
    "CHUK", "IKBKB",
    "IRF3", "NFKB1", "RELA"
  ))
## check cGAS Sting genes
topDV_significant %>%
  dplyr::filter(gene_symbol %in% c(
    "CGAS", "STING1",
    "TBK1",
    "CHUK", "IKBKB",
    "IRF3", "NFKB1", "RELA"
  ))
# topDV_significant filter output: only NFKB1 and IKBKB

# plot NFKB1 and IKBKB genes
IRGI_DE <- lcpm_df %>%
  dplyr::mutate(
    gene_id = str_replace_all(gene_name, "_.+", ""),
    gene_name = str_replace_all(gene_name, ".+_", "")
  ) %>%
  dplyr::left_join(PMM_sampleSheet, by = c("sample" = "sampleName")) %>%
  # dplyr::filter(gene_name %in% c("NFKB1", "IKBKB")) %>%
  dplyr::filter(gene_name %in% c("NFKB1")) %>%
  mutate(sample2 = case_when(
    str_detect(sample, "^BMPC\\d+$") ~ "BMPC",
    str_detect(sample, "^CD19") ~ "BMPC",
    TRUE ~ "PMM"
  )) %>%
  ggplot(aes(x = sample2, y = lcpm)) +
  labs(x = "") +
  # geom_boxplot(outlier.shape = NA, alpha = 0.9, position = position_dodge(0.8)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.15), alpha = 0.5) +
  facet_wrap(~gene_name, scales = "free_x", nrow = 2) +
  # ylim(c(-20,150))+
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", vjust = 1) +
  ggpubr::theme_pubr(base_size = 9.6)
ggsave(paste0(output_dir, "/DE_DV_IRGI_DE.pdf"),
  plot = IRGI_DE,
  width = 10, height = 8, units = "cm", limitsize = FALSE
)

IRGI_DV <- lcpm_df %>%
  dplyr::mutate(
    gene_id = str_replace_all(gene_name, "_.+", ""),
    gene_name = str_replace_all(gene_name, ".+_", "")
  ) %>%
  dplyr::left_join(PMM_sampleSheet, by = c("sample" = "sampleName")) %>%
  dplyr::filter(gene_name %in% c("NFKB1", "IKBKB")) %>%
  dplyr::mutate(gene_name = factor(gene_name, levels = c("NFKB1", "IKBKB"))) %>%
  mutate(sample2 = case_when(
    str_detect(sample, "^BMPC\\d+$") ~ "BMPC",
    str_detect(sample, "^CD19") ~ "BMPC",
    TRUE ~ sample
  )) %>%
  ggplot(aes(x = lcpm, y = L1PA2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~gene_name, scales = "free_x") +
  # ylim(c(-20,150))+
  ggpubr::stat_cor(label.y = c(135, 135), method = "pearson", size = 3) +
  # ggrepel::geom_text_repel(aes(label = sample2), max.overlaps = 2,
  #                          box.padding = 0.5, size = 5)+
  ggrepel::geom_text_repel(
    aes(label = ifelse(str_detect(sample2, "PMM"), sample2, NA_character_)),
    max.overlaps = Inf, box.padding = 0.5, size = 5, na.rm = TRUE
  ) +
  ggpubr::theme_pubr(base_size = 9.6)
ggsave(paste0(output_dir, "/DE_DV_IRGI_DV.pdf"),
  plot = IRGI_DV,
  width = 10, height = 8, units = "cm", limitsize = FALSE
)


# volcano plot to show top var KZFPs #################
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
})

# =========================
# 0) 你需要准备的数据
# =========================
# topDV: 至少包含这三列：
#   gene_symbol (chr), LogVarRatio (num), FDR (num)
# kzfp_list: KZFP 基因 symbol 向量（chr）
# epi_list : epigenetic regulator 基因 symbol 向量（chr）

# ---- 例子（你要替换成你自己的） ----
# kzfp_list <- readLines("kzfp_symbols.txt")
# epi_list  <- readLines("epi_regulators_symbols.txt")

# =========================
# 1) 参数（可改）
# =========================
fdr_cut <- 0.05
lfc_cut <- 2

# 每组希望总共标多少个（默认 10 = Up5 + Down5）
N_total_per_set <- 10
N_each_dir_max <- 5 # Up/Down 每边最多取 5

# top 的定义：按 “FDR最小” 或 “|log2FC|最大”
rank_by <- "FDR" # 可选： "FDR" 或者 "abs_log2FC"
# =========================

# =========================
# 2) 数据清洗 + volcano 坐标
# =========================
df <- topDV %>%
  mutate(
    gene_symbol = as.character(gene_symbol),
    log2FC = as.numeric(LogVarRatio),
    FDR = as.numeric(FDR)
  ) %>%
  filter(!is.na(gene_symbol), !is.na(log2FC), !is.na(FDR), FDR > 0) %>%
  mutate(
    neg_log10_FDR = -log10(FDR),
    sig_dir = case_when(
      FDR < fdr_cut & log2FC >= lfc_cut ~ "Up",
      FDR < fdr_cut & log2FC <= -lfc_cut ~ "Down",
      TRUE ~ "NS"
    ),
    is_kzfp = gene_id %in% KZFPs$ensembl_gene_id,
    is_epi = gene_id %in% ER_genes$gene_id
  )

# =========================
# 3) 通用函数：尽量凑够 N_total（Up/Down 不足会从另一边补）
# =========================
pick_top_ud_fill <- function(dat,
                             n_total = 10,
                             n_each_dir_max = 5,
                             rank_by = c("FDR", "abs_log2FC")) {
  rank_by <- match.arg(rank_by)

  dat2 <- dat %>% filter(sig_dir %in% c("Up", "Down"))
  if (nrow(dat2) == 0) {
    return(dat2)
  } # 直接返回空表（列结构保留）

  # 定义排序 key：越靠前越“top”
  if (rank_by == "FDR") {
    dat2 <- dat2 %>% arrange(FDR)
  } else {
    dat2 <- dat2 %>% arrange(desc(abs(log2FC)))
  }

  up_all <- dat2 %>% filter(sig_dir == "Up")
  down_all <- dat2 %>% filter(sig_dir == "Down")

  # 先各取最多 n_each_dir_max
  up_pick <- up_all %>% slice_head(n = min(n_each_dir_max, nrow(up_all)))
  down_pick <- down_all %>% slice_head(n = min(n_each_dir_max, nrow(down_all)))

  picked <- bind_rows(up_pick, down_pick) %>%
    distinct(gene_symbol, .keep_all = TRUE)

  # 如果已经够了（或总数本来就不足 n_total），直接返回
  if (nrow(picked) >= n_total || nrow(dat2) <= n_total) {
    # 如果 picked 超过 n_total（理论上不会，因为最多 5+5=10），还是保险裁一下
    return(picked %>% slice_head(n = min(n_total, nrow(picked))))
  }

  # 还不够：从剩下的里补齐（按同样排序）
  remain <- dat2 %>%
    filter(!gene_symbol %in% picked$gene_symbol)

  need <- n_total - nrow(picked)
  fill <- remain %>% slice_head(n = min(need, nrow(remain)))

  bind_rows(picked, fill) %>%
    distinct(gene_symbol, .keep_all = TRUE)
}

# =========================
# 4) 选三套要标注的基因（每套尽量 10 个）
# =========================
top_global <- pick_top_ud_fill(df,
  n_total = N_total_per_set,
  n_each_dir_max = N_each_dir_max,
  rank_by = rank_by
) %>%
  mutate(label_group = "Global top")

top_kzfp <- pick_top_ud_fill(df %>% filter(is_kzfp),
  n_total = N_total_per_set,
  n_each_dir_max = N_each_dir_max,
  rank_by = rank_by
) %>%
  mutate(label_group = "KZFP top")

top_epi <- pick_top_ud_fill(df %>% filter(is_epi),
  n_total = N_total_per_set,
  n_each_dir_max = N_each_dir_max,
  rank_by = rank_by
) %>%
  mutate(label_group = "EpiReg top")

# 合并所有要标注的基因（去重）
top_genes <- bind_rows(top_global, top_kzfp, top_epi) %>%
  distinct(gene_symbol, .keep_all = TRUE)

# 可选：看看每组实际取到了几个
cat("Labeled Global:", nrow(top_global), "\n")
cat("Labeled KZFP  :", nrow(top_kzfp), "\n")
cat("Labeled EpiReg:", nrow(top_epi), "\n")

# =========================
# 5) 给点上色分类（KZFP / Epi / Other）
#    如果一个基因同时在 KZFP 和 Epi（很少），这里优先按 KZFP 归类
# =========================
df <- df %>%
  mutate(
    highlight = case_when(
      is_kzfp & sig_dir == "Up" ~ "KZFP Up",
      is_kzfp & sig_dir == "Down" ~ "KZFP Down",
      !is_kzfp & is_epi & sig_dir == "Up" ~ "EpiReg Up",
      !is_kzfp & is_epi & sig_dir == "Down" ~ "EpiReg Down",
      TRUE ~ "Other genes"
    )
  )

# =========================
# 6) 画 volcano（先灰点，再彩点）
# =========================
top_genes <- top_genes %>%
  dplyr::filter(sig_dir == "Up" & label_group %in% c("KZFP top", "EpiReg top"))
p <- ggplot(df, aes(x = log2FC, y = neg_log10_FDR, color = highlight)) +
  geom_point(
    data = subset(df, highlight == "Other genes"),
    size = 1, alpha = 0.5
  ) +
  geom_point(
    data = subset(df, highlight != "Other genes"),
    size = 2.8, alpha = 0.9
  ) +
  scale_color_manual(
    values = c(
      "KZFP Up"     = "red",
      "KZFP Down"   = "blue",
      "EpiReg Up"   = "orange",
      "EpiReg Down" = "deepskyblue3",
      "Other genes" = "grey80"
    ),
    labels = c(
      "KZFP Up"     = paste0("KZFP Up (n=", sum(df$highlight == "KZFP Up"), ")"),
      "KZFP Down"   = paste0("KZFP Down (n=", sum(df$highlight == "KZFP Down"), ")"),
      "EpiReg Up"   = paste0("EpiReg Up (n=", sum(df$highlight == "EpiReg Up"), ")"),
      "EpiReg Down" = paste0("EpiReg Down (n=", sum(df$highlight == "EpiReg Down"), ")"),
      "Other genes" = "Other genes"
    )
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  theme_pubr(base_size = 9.6) +
  labs(
    x = "PMM vs BMPC: log2(Variation change)",
    y = "-log10(FDR)",
    color = "Highlighted genes"
  ) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = top_genes,
    inherit.aes = FALSE,
    aes(x = log2FC, y = neg_log10_FDR, label = gene_symbol),
    size = 3.6,
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf,
    min.segment.length = 0
  )


print(p)

ggsave(paste0(output_dir, "/volcano_KZFP_Epi_GlobalTopDV_markUp.pdf"),
  plot = p,
  width = 8, height = 8
)


## GSEA ########################################################################
library(msigdbr)
as.data.frame(msigdbr_collections())
m_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- m_df[, c("gs_name", "ensembl_gene")] # "human_gene_symbol",

### Create a ranked list
ranked_list <- dat %>%
  dplyr::select(gene_id, logFC) %>%
  arrange(desc(logFC)) %>%
  deframe() # Convert to named vector (gene names as names, logFC as values)

## rank result by p values
dum <- dat$FDR
## Avoid Inf when taking log of 0
# dum[dum == 0] <- 1e-300
score <- -log10(dum)
# It is important to differentiate up and down regulation
score <- ifelse(dat$logFC >= 0, score, -score)
names(score) <- dat$gene_id
score <- sort(score, decreasing = TRUE)
## Add jitter to avoid ties
# set.seed(52)
# score <- score + rnorm(length(score), mean = 0, sd = 1e-6)
# score <- sort(score, decreasing = TRUE)

## perform GSEA
library(clusterProfiler)
em_1 <- GSEA(ranked_list,
  TERM2GENE = hallmark, eps = 0,
  pvalueCutoff = 0.05, pAdjustMethod = "BH"
)

em_2 <- GSEA(score,
  TERM2GENE = hallmark, eps = 0,
  pvalueCutoff = 0.05, pAdjustMethod = "BH"
)

# Create a dotplot of top 30 pathways
dotplot(em_1, orderBy = "NES", x = "NES", showCategory = 20)
dotplot(em_2, orderBy = "NES", x = "NES", showCategory = 20)

dotplot(em_1, orderBy = "NES", x = "NES", showCategory = 20) +
  # coord_flip() +
  scale_color_gradient(low = "red", high = "blue") +
  # scale_size(range = c(2, 10)) +
  labs(
    title = "",
    x = "Normalized Enrichment Score (NES)",
    y = "",
    # color = "Adjusted p-value",
    size = "Gene Count"
  ) +
  theme_pubr(base_size = 12, legend = "right") +
  theme(axis.text.y = element_text(size = 9.6))


ggsave(paste0(output_dir, "/lmVoomLfcRanked_gsea_NES_top20_adjustDataSetBatches_v.pdf"),
  width = 8, height = 12
)

library(org.Hs.eg.db)
# all_gene <- data$gene_id %>% unique()
# all_gene <- str_remove_all(all_gene, "\\..*$")
ego <- enrichGO(
  # gene          = dat[dat$Sig == "Up",]$gene_id,
  gene          = dat[dat$Sig != "NS", ]$gene_id,
  # gene          = dat[dat$Sig != "Down",]$gene_id,
  universe      = dat$gene_id, # optional but recommended
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL", # or "ENSEMBL" etc, SYMBOL match your IDs
  ont           = "ALL", # "BP", "MF", or "CC"
  # pvalueCutoff  = 0.25,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  readable      = TRUE
)

ego_s <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
# dotplot(ego_s, showCategory = 10, split = "ONTOLOGY") +
#   facet_grid(ONTOLOGY ~ ., scales = "free_y")


p_ego_s <- dotplot(ego_s, showCategory = 10, split = "ONTOLOGY") +
  # facet_grid(ONTOLOGY ~ ., scales = "free_y", cols = 3) +
  facet_wrap(vars(ONTOLOGY), ncol = 3, scales = "free_y") +
  # scale_fill_continuous()
  scale_color_gradient(low = "red", high = "blue") +
  # scale_size(range = c(2, 10)) +
  labs(
    title = "",
    x = "GeneRatio",
    y = "",
    # color = "Adjusted p-value",
    size = "Gene Count"
  ) +
  theme_pubr(base_size = 9.6, legend = "right") +
  theme(axis.text.y = element_text(size = 9.6))
# p <- dotplot(ego, showCategory = 10)
p_ego_s
ggsave(paste0(output_dir, "/lmVoom_allDE_GO_top10_adjustDataSetBatches.pdf"),
  plot = p_ego_s,
  width = 15, height = 8
)


# plot lcpm of some genes:
p1 <-
  lcpm_df %>%
  dplyr::mutate(
    gene_id = str_replace_all(gene_name, "_.+", ""),
    gene_name = str_replace_all(gene_name, ".+_", "")
  ) %>%
  dplyr::filter(gene_name %in% c(
    "ZNF141", "MPHOSPH8",
    "ATF7IP", "SIRT1",
    "TASOR", "MORC2", "PPHLN1", "SETDB1",
    "TASOR2", "IRF2",
    "DNMT1", "DNMT3A"
  )) %>%
  # distinct() %>%
  dplyr::left_join(PMM_sampleSheet, by = c("sample" = "sampleName")) %>%
  # ggplot(aes(x=sample, y=lcpm, fill=gene_name))+
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm))+
  ggplot(aes(x = L1PA2, y = lcpm)) +
  geom_point() +
  # geom_bar(stat = "identity")+
  ggrepel::geom_text_repel(aes(label = sample), size = 3) +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(method = "pearson") +
  # geom_bar(stat = "identity", position = "dodge")+
  labs(title = "HUSH complex genes", x = "") +
  facet_wrap(~gene_name) +
  # facet_wrap(~gene_name+transcript_id, scales = "free_y")+
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(paste0(output_dir, "/lcpm_HUSHgenes_R3_STARcounts.pdf"),
  plot = p1,
  width = 25, height = 20, units = "cm", limitsize = FALSE
)


# Create a data frame with genes and their annotations
key_gene_df <- data.frame(
  Gene = c(
    "ZNF141", "ZNF382", "ZNF429", "ZNF483", "ZNF860",
    "NFKB1", "JAK1", "IFIH1", "IRF1",
    "SNAI2", "TWIST1",
    "E2F1", "E2F2", "E2F5", "E2F7", "E2F8", "RB1", "CDK1",
    "SUZ12", "EZH2", "PRC1", "EED",
    "TP53", "MDM2"
  ),
  Category = c(
    "Young L1 ZNF", "Young L1 ZNF", "Young L1 ZNF",
    "Stem Cell ZNF", "Sig.Down ZNF in PMM",
    "NFKB", "NFKB", "IFN", "IFN",
    "EMT", "EMT",
    "E2F", "E2F", "E2F", "E2F", "E2F", "E2F", "E2F",
    "PRC2", "PRC2", "PRC2", "PRC2",
    "TP53 Pathway", "TP53 Pathway"
  )
)

# NFKB(NFKB -> proinflammatory cytokines), IFN(IRF3 -> type I IFNs)
# PRC2, EZH2 bind to young L1s
p2 <-
  lcpm_df %>%
  dplyr::mutate(
    gene_id = str_replace_all(gene_name, "_.+", ""),
    gene_name = str_replace_all(gene_name, ".+_", "")
  ) %>%
  dplyr::filter(gene_name %in% key_gene_df$Gene) %>%
  dplyr::left_join(key_gene_df, by = c("gene_name" = "Gene")) %>%
  dplyr::left_join(PMM_sampleSheet, by = c("sample" = "sampleName")) %>%
  # ggplot(aes(x=sample_index, y=lcpm, fill=gene_name))+
  ggplot(aes(x = L1PA2, y = lcpm)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  labs(title = "key genes (ZNF, NFKB, EMT, E2F, PRC and TP53)", x = "") +
  facet_wrap(~ Category + gene_name, scales = "free_y", ncol = 5) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(paste0(output_dir, "/lcpm_keygenes_R2_STARcounts.pdf"),
  plot = p2,
  width = 25, height = 35, units = "cm", limitsize = FALSE
)


# check infilteration cor ######################################################
# v_obj <- voom(dge, design)
# E <- v_obj$E  # voom 的 logCPM matrix (genes x samples)
# pc <- prcomp(t(E), scale. = TRUE)
pc <- prcomp(t(normalized_lcpm_bc), scale. = TRUE)

df_pc <- data.frame(
  PC1 = pc$x[, 1], PC2 = pc$x[, 2],
  infiltration = PMM_sampleSheet$Infilteration,
  group = PMM_sampleSheet$group,
  L1PA2 = PMM_sampleSheet$L1PA2,
  ZNF141 = normalized_lcpm_bc[str_detect(rownames(normalized_lcpm_bc), pattern = "_ZNF141"), ]
)

# cor.test(PMM_sampleSheet$Infilteration, normalized_lcpm[str_detect(rownames(normalized_lcpm), pattern = "ZNF141"),])


ggplot(df_pc, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(title = "PCA (voom logCPM)", color = "Infiltration")

# cor.test(df_pc$PC1, df_pc$infiltration, method="pearson")
# cor.test(df_pc$PC2, df_pc$infiltration, method="pearson")
# cor.test(df_pc$L1PA2, df_pc$infiltration, method="pearson")
cor.test(df_pc$L1PA2, df_pc$ZNF141, method = "pearson")

df_pmm <- subset(df_pc, group == "PMM")
df_pmm$sample <- rownames(df_pmm)


# cor.test(df_pmm$PC1, df_pmm$infiltration, method = "spearman")
# cor.test(df_pmm$PC2, df_pmm$infiltration, method = "spearman")
#
# ct1 <- cor.test(df_pmm$PC1, df_pmm$infiltration, method="spearman", exact=FALSE)
# ct2 <- cor.test(df_pmm$PC2, df_pmm$infiltration, method="spearman", exact=FALSE)
#
# lab1 <- sprintf("Spearman ρ = %.3f\np = %.3g", ct1$estimate, ct1$p.value)
# lab2 <- sprintf("Spearman ρ = %.3f\np = %.3g", ct2$estimate, ct2$p.value)

# Pearson correlation tests
cor.test(df_pmm$PC1, df_pmm$infiltration, method = "pearson")
cor.test(df_pmm$PC2, df_pmm$infiltration, method = "pearson")
cor.test(df_pmm$L1PA2, df_pmm$infiltration, method = "pearson")
cor.test(df_pmm$ZNF141, df_pmm$infiltration, method = "pearson")
cor.test(df_pmm$L1PA2, df_pmm$ZNF141, method = "pearson")

ct1 <- cor.test(df_pmm$PC1, df_pmm$infiltration, method = "pearson")
ct2 <- cor.test(df_pmm$PC2, df_pmm$infiltration, method = "pearson")
ct3 <- cor.test(df_pmm$L1PA2, df_pmm$infiltration, method = "pearson")

ct4 <- cor.test(df_pmm$ZNF141, df_pmm$infiltration, method = "pearson")
ct5 <- cor.test(df_pmm$L1PA2, df_pmm$ZNF141, method = "pearson")

lab1 <- sprintf("Pearson r = %.3f\np = %.3g", ct1$estimate, ct1$p.value)
lab2 <- sprintf("Pearson r = %.3f\np = %.3g", ct2$estimate, ct2$p.value)
lab3 <- sprintf("Pearson r = %.3f\np = %.3g", ct3$estimate, ct3$p.value)

lab4 <- sprintf("Pearson r = %.3f\np = %.3g", ct4$estimate, ct4$p.value)
lab5 <- sprintf("Pearson r = %.3f\np = %.3g", ct5$estimate, ct5$p.value)


p1 <- ggplot(df_pmm, aes(infiltration, PC1)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf) +
  annotate("text",
    x = Inf, y = Inf, label = lab1,
    hjust = 1.1, vjust = 1.1, size = 3.5
  ) +
  theme_classic() +
  labs(title = "PMM RNAseq: PC1 vs Infiltration")

p2 <- ggplot(df_pmm, aes(infiltration, PC2)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf) +
  annotate("text",
    x = Inf, y = Inf, label = lab2,
    hjust = 1.1, vjust = 2, size = 3.5
  ) +
  theme_classic() +
  labs(title = "PMM RNAseq: PC2 vs Infiltration")

p3 <- ggplot(df_pmm, aes(infiltration, L1PA2)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf) +
  annotate("text",
    x = Inf, y = Inf, label = lab3,
    hjust = 1.1, vjust = 2.5, size = 3.5
  ) +
  theme_classic() +
  labs(title = "PMM RNAseq: L1PA2 vs Infiltration")

p4 <- ggplot(df_pmm, aes(infiltration, ZNF141)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf) +
  annotate("text",
    x = Inf, y = Inf, label = lab4,
    hjust = 1.1, vjust = 2.5, size = 3.5
  ) +
  theme_classic() +
  labs(title = "PMM RNAseq: ZNF141 vs Infiltration")

p5 <- ggplot(df_pmm, aes(L1PA2, ZNF141)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf) +
  annotate("text",
    x = Inf, y = Inf, label = lab5,
    hjust = 1.1, vjust = 2.5, size = 3.5
  ) +
  theme_classic() +
  labs(title = "PMM RNAseq: L1PA2 vs ZNF141")

# pw <- p1 / p2 / p3
pw <- p1 / p2
pw

## 4) save patchwork to pdf (vector)
ggsave(
  filename = paste0(output_dir, "/pmm_RNAseq_pcL1_infiltraCor_v2.pdf"),
  plot = pw,
  width = 5,
  height = 6
)

# pw2 <- p4/p5
pw2 <- p3 / p5

ggsave(
  filename = paste0(output_dir, "/pmm_RNAseq_pcL1_infiltraCor_ZNF141_v2.pdf"),
  plot = pw2,
  width = 5,
  height = 6
)


# ------------------------------------------------------------------------------
# Add L1 number vs ZHX1 and SMARCA2 gene cor plot
# ------------------------------------------------------------------------------
p_zhx1_smarca2 <- lcpm_df %>%
  dplyr::mutate(
    gene_id = str_replace_all(gene_name, "_.+", ""),
    gene_name = str_replace_all(gene_name, ".+_", "")
  ) %>%
  dplyr::left_join(PMM_sampleSheet, by = c("sample" = "sampleName")) %>%
  dplyr::filter(gene_name %in% c("ZHX1", "SMARCA2")) %>%
  # dplyr::mutate(gene_name = factor(gene_name, levels = c("ZHX1", "SMARCA2"))) %>%
  mutate(sample2 = case_when(
    str_detect(sample, "^BMPC\\d+$") ~ "BMPC",
    str_detect(sample, "^CD19") ~ "BMPC",
    TRUE ~ sample
  )) %>%
  ggplot(aes(x = L1PA2, y = lcpm)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~gene_name, scales = "free_y") +
  ggpubr::stat_cor(method = "pearson", size = 3) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(str_detect(sample2, "PMM"), sample2, NA_character_)),
    max.overlaps = Inf, box.padding = 0.5, size = 3, na.rm = TRUE
  ) +
  ggpubr::theme_pubr(base_size = 9.6) +
  labs(title = "PMM RNAseq: L1PA2 vs ZHX1/SMARCA2")

ggsave(paste0(output_dir, "/L1PA2_vs_ZHX1_SMARCA2_cor.pdf"),
  plot = p_zhx1_smarca2,
  width = 15, height = 8, units = "cm", limitsize = FALSE
)
