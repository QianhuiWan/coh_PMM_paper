# use count matrix from feature counts
# filterByExpr, TMM normalize
# get gene counts stats, MMRF, use limma voom method
# and then get DE topTable

dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s1_checkStrandness_R4_gemini_uniqueMap"
# output_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR2"
output_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3_uniqueMap/filtered_L1_redo_res"
# input_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s1_checkStrandness_R4"
# dir_gemini = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s1_checkStrandness_R4_gemini"

# load packages
library(tidyverse)
library(magrittr)
library(stringr)
library(rtracklayer)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(DESeq2)
library(edgeR)
library(biomaRt)
library(readxl)
library(clusterProfiler)
library(msigdbr)
library(enrichR)
library(enrichplot)
library(ggplot2)

# load sample sheet: sample in 4 groups by activated L1
# patientsIn4Groups <-
#   read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3_uniqueMap/filtered_L1_redo_res/all_activeL1PA2_patientsIn4Groups_764.tsv")

patientsIn4Groups <-
  read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3_uniqueMap/filtered_L1_redo_res/all_activeL1s_patientsIn4Groups_764.tsv")

# MM_metadata <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/clinical_data_BMprimary_addID.tsv")
# MM_metadata_update <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/clinical_data_BMprimary_addID_addSubtypes_uniqueMap.tsv")
MM_metadata_update <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/clinical_data_BMprimary_addID_addL1_addSubtypes_uniqueMap.tsv")

MM_metadata_update <- MM_metadata_update[match(patientsIn4Groups$sampleID, MM_metadata_update$SRR_id), ]
table(patientsIn4Groups$sampleID == MM_metadata_update$SRR_id)

patientsIn4Groups$low_purity <- MM_metadata_update$`Low purity`
patientsIn4Groups$sex <- MM_metadata_update$gender
patientsIn4Groups$age <- MM_metadata_update$age_at_index

ER_genes <- read_excel("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/pub_data/epigenetic_regulators/KEPI_A_2139067_SM6382.xlsx",
  sheet = 2, skip = 2
) %>%
  `colnames<-`(c("gene_id", "gene_symbol", "class"))

## use getNum FUN to get all transcripts #######################################
# Save as a .tsv file
# write.table(count_matrix, file =  paste0(dir, "/gene_count_matrix.tsv"),
#             sep = "\t", row.names = TRUE, col.names = TRUE)
# load count matrix
count_matrix <-
  read.table(paste0(dir, "/gene_count_matrix_STARcounts_MM_764.tsv"),
    sep = "\t", header = TRUE, row.names = 1
  ) %>%
  as.matrix()

# write.table(count_matrix, file =  paste0(dir, "/gene_count_matrix_STARcounts.tsv"),
#             sep = "\t", row.names = TRUE, col.names = TRUE)

# DE
## match meta with matrix column names (i.e. sample IDs)
# count_matrix <-
#   count_matrix[,colnames(count_matrix)%in%patientsIn4Groups$sampleID]
# count_matrix <-
#   count_matrix[, match(patientsIn4Groups$sampleID, colnames(count_matrix))]
table(colnames(count_matrix) == patientsIn4Groups$sampleID)
table(colnames(count_matrix) == MM_metadata_update$SRR_id)

dge <- DGEList(
  counts = count_matrix,
  group = patientsIn4Groups$L1_quantile_group
)

# Define a threshold for filtering
# Keep genes with CPM > 1 in at least 1 group
keep <- rowSums(cpm(dge, normalized.lib.sizes = TRUE) >= 1) >= min(table(patientsIn4Groups$L1_quantile_group))
# keep <- rowSums(cpm(dge) > 2) >= 2
# keep <- rowSums(cpm(dge) >= 2) >= floor(870*0.75)
# keep <- filterByExpr(dge, group=dge$samples$group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
# normalize and filter
normalized_lcpm <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)
# save filtered normalized lcpm for later plotting
# write.table(normalized_lcpm, file =  paste0(output_dir, "/lcpm_TMM_fil.tsv"),
#             sep = "\t", row.names = TRUE, col.names = TRUE)


# cor.test(MM_metadata_update$`Low purity`, normalized_lcpm[str_detect(rownames(normalized_lcpm), pattern = "ZNF425"),])
# correct for batch after normalization
idx <- match(colnames(normalized_lcpm), MM_metadata_update$SRR_id)
lp <- MM_metadata_update$`Low purity`[idx]
age <- MM_metadata_update$age_at_index[idx]


# 2) keep samples with complete covariates
keep <- !is.na(lp) & !is.na(age)

expr2 <- normalized_lcpm[, keep, drop = FALSE]

design2 <- model.matrix(~ 0 + patientsIn4Groups$L1_quantile_group[idx][keep])

# 3) build covariate matrix (must be numeric, no NA)
Xcov <- cbind(low_purity = lp[keep])

normalized_lcpm_bc <- removeBatchEffect(
  expr2,
  covariates = Xcov,
  design = design2
)

lcpm_df <- normalized_lcpm_bc %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_name") %>%
  pivot_longer(
    cols = c(starts_with("SRR")),
    names_to = "sample",
    values_to = "lcpm"
  ) # %>%
# mutate(lcpm = ifelse(lcpm < 0, 0, lcpm))

readr::write_tsv(lcpm_df, file = paste0(output_dir, "/lcpm_TMM_MM_STARcounts_l1groups_bc.tsv"))


lcpm_df <- expr2 %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_name") %>%
  pivot_longer(
    cols = c(starts_with("SRR")),
    names_to = "sample",
    values_to = "lcpm"
  ) # %>%
# mutate(lcpm = ifelse(lcpm < 0, 0, lcpm))

readr::write_tsv(lcpm_df, file = paste0(output_dir, "/lcpm_TMM_MM_STARcounts.tsv"))

## design matrix
dge <- dge[, !is.na(patientsIn4Groups$low_purity)]

design <- model.matrix(~ 0 + L1_quantile_group + low_purity, data = patientsIn4Groups[!is.na(patientsIn4Groups$low_purity), ])
## make contrast
con <- makeContrasts(
  group = L1_quantile_groupQ4 - L1_quantile_groupQ1,
  levels = design
) # TumorS - TumorM

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

plotMD(tfit,
  column = 1, status = dt[, 1], main = colnames(tfit)[1],
  xlim = c(-8, 13)
)

## get topTable output
dat <- topTable(v, coef = "group", number = Inf) %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  dplyr::mutate(gene_id_version = str_replace_all(gene_name, "_.+", "")) %>%
  dplyr::mutate(gene_id = str_remove(gene_id_version, "\\..*")) %>%
  dplyr::mutate(gene_symbol = str_replace_all(gene_name, ".+_", "")) %>%
  # distinct(gene_name, .keep_all = TRUE) %>%
  dplyr::select(-t, -B, FDR = `adj.P.Val`) %>%
  mutate(Sig = ifelse(FDR < 0.05 & logFC > 1, "Up",
    ifelse(FDR < 0.05 & logFC < -1, "Down", "NS")
  ))
table(dat$Sig)

# get geneID from geneID version
library(biomaRt)
## load mart: #i.e. Human genes (GRCh38.p14)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
## get gene ids
valid_genes <- getBM(
  attributes = c(
    "hgnc_symbol", "ensembl_gene_id",
    "ensembl_gene_id_version"
  ),
  filters = "ensembl_gene_id_version",
  values = dat$gene_id_version,
  mart = ensembl
)
# update dat
dat <-
  dat %>%
  dplyr::left_join(valid_genes,
    by = c("gene_id_version" = "ensembl_gene_id_version")
  )

write_tsv(dat, file = paste0(output_dir, "/lmVoom_topTable_L1groupsQ4Q1_dat_purityAdj.tsv"))


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
# dat_DE filter output: 0


## diff var analysis
# DiffVar
library(missMethyl)
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
  mutate(Sig = ifelse(FDR < 0.05 & LogVarRatio > 1, "Up",
    ifelse(FDR < 0.05 & LogVarRatio < -1, "Down", "NS")
  ))

head(topDV)
table(topDV$Sig) # 0


## GSEA
library(msigdbr)
library(clusterProfiler)
as.data.frame(msigdbr_collections())
m_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- m_df[, c("gs_name", "db_ensembl_gene")] # "human_gene_symbol",

### Create a ranked list
ranked_list <- dat %>%
  dplyr::select(gene_id, logFC) %>%
  arrange(desc(logFC)) %>%
  deframe() # Convert to named vector (gene names as names, logFC as values)

## rank result by p values
dum <- dat$FDR
# Avoid Inf when taking log of 0
dum[which(dum == 0)] <- 1e-300
score <- -log10(dum)
# It is important to differentiate up and down regulation
score <- ifelse(dat$logFC >= 0, score, -1 * score)
names(score) <- dat$gene_id
score <- sort(score, decreasing = T)
set.seed(51)
score <- score + rnorm(length(score), mean = 0, sd = 1e-6)
score <- sort(score, decreasing = TRUE)

em_1 <- GSEA(ranked_list,
  TERM2GENE = hallmark, eps = 0,
  pvalueCutoff = 0.05, pAdjustMethod = "BH"
)

em_2 <- GSEA(score,
  TERM2GENE = hallmark,
  pvalueCutoff = 0.05, pAdjustMethod = "BH"
)


# Create a dotplot of top 30 pathways
dotplot(em_1, orderBy = "NES", x = "NES", showCategory = 20)
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_NES_top20_purityAdj.pdf"), width = 9, height = 12)

dotplot(em_2, orderBy = "NES", x = "NES", showCategory = 20) + ggpubr::theme_pubr()
ggsave(paste0(output_dir, "/lmVoomPval_gsea_NES_top20.pdf"), width = 9, height = 12)

em_plot <- dotplot(em_1, orderBy = "NES", x = "NES", showCategory = 20) +
  # coord_flip() +
  scale_color_gradient(low = "red", high = "blue") +
  scale_size(range = c(2, 10)) +
  labs(
    title = "",
    x = "Normalized Enrichment Score (NES)",
    y = "",
    fill = "Adjusted p-value",
    size = "Gene count"
  ) +
  theme_pubr(base_size = 12, legend = "right") +
  theme(axis.text.y = element_text(size = 9.6))

ggsave(paste0(output_dir, "/lmVoomLfc_gsea_NES_top20_L1groups_purityAdj.pdf"),
  plot = em_plot,
  width = 8, height = 12
)


# dotplot(em_1, orderBy="NES", x="NES", showCategory = 25)
# ggsave(paste0(output_dir, "/lmVoomLfc_gsea_NES_top25.pdf"), width = 9, height = 12)
#
# dotplot(em_2, orderBy="NES", x="NES", showCategory = 25)
# ggsave(paste0(output_dir, "/lmVoomPval_gsea_NES_top25.pdf"), width = 9, height = 12)

# plot GSEA plot for key pathways
gseaplot(em_1, geneSetID = "HALLMARK_P53_PATHWAY", title = "HALLMARK_P53_PATHWAY")
gseaplot2(em_1, geneSetID = "HALLMARK_P53_PATHWAY", title = "HALLMARK_P53_PATHWAY")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_P53_PATHWAY.pdf"), width = 8, height = 6)

gseaplot(em_1, geneSetID = "HALLMARK_E2F_TARGETS", title = "HALLMARK_E2F_TARGETS")
gseaplot2(em_1, geneSetID = "HALLMARK_E2F_TARGETS", title = "HALLMARK_E2F_TARGETS")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_E2F_TARGETS.pdf"), width = 8, height = 6)

gseaplot(em_1, geneSetID = "HALLMARK_G2M_CHECKPOINT", title = "HALLMARK_G2M_CHECKPOINT")
gseaplot2(em_1, geneSetID = "HALLMARK_G2M_CHECKPOINT", title = "HALLMARK_G2M_CHECKPOINT")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_G2M_CHECKPOINT.pdf"), width = 8, height = 6)


gseaplot(em_1, geneSetID = "HALLMARK_INTERFERON_ALPHA_RESPONSE", title = "HALLMARK_INTERFERON_ALPHA_RESPONSE")
gseaplot2(em_1, geneSetID = "HALLMARK_INTERFERON_ALPHA_RESPONSE", title = "HALLMARK_INTERFERON_ALPHA_RESPONSE")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_INTERFERON_ALPHA_RESPONSE.pdf"), width = 8, height = 6)

gseaplot(em_1, geneSetID = "HALLMARK_INTERFERON_GAMMA_RESPONSE", title = "HALLMARK_INTERFERON_GAMMA_RESPONSE")
gseaplot2(em_1, geneSetID = "HALLMARK_INTERFERON_GAMMA_RESPONSE", title = "HALLMARK_INTERFERON_GAMMA_RESPONSE")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_INTERFERON_GAMMA_RESPONSE.pdf"), width = 8, height = 6)

gseaplot(em_1, geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB", title = "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
gseaplot2(em_1, geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB", title = "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_TNFA_SIGNALING_VIA_NFKB.pdf"), width = 8, height = 6)


# plot gene expression ##############################
MM_meta_update <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/mmrf_subtypes/clinical_data_BMprimary_addID_addSubtypes_uniqueMap.tsv")

patientsIn4Groups <-
  read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3_uniqueMap/filtered_L1_redo_res/all_activeL1PA2_patientsIn4Groups_764.tsv")

MM_meta_final <- MM_meta_update %>%
  dplyr::select(-(sample_index_antisenseL1:L1_quantile_group)) %>%
  dplyr::left_join(patientsIn4Groups, by = c("SRR_id" = "sampleID"))

PR_subtypeIDs <- MM_meta_update[MM_meta_update$RNA_Subtype_Name == "PR", ]$sampleID


lcpm_MM <- lcpm_df %>%
  # rownames_to_column(var="gene_name") %>%
  dplyr::mutate(gene_id_version = str_replace_all(gene_name, "_.+", "")) %>%
  dplyr::mutate(gene_id = str_remove(gene_id_version, "\\..*")) %>%
  dplyr::mutate(gene_symbol = str_replace_all(gene_name, ".+_", "")) # %>%
# dplyr::filter(!(gene_id == "ENSG00000228314")) %>%  # redundant Pseudogene id
# pivot_longer(cols = -c(gene_name, gene_id_version, gene_id, gene_symbol),
#              names_to = "sampleID", values_to = "lcpm") #%>%
# dplyr::mutate(lcpm = ifelse(lcpm < 0, 0, lcpm))

# long to wide
lcpm_MM_wide <- lcpm_MM %>%
  # dplyr::select(gene_name = gene_name, sample = sampleID, lcpm) %>%
  dplyr::select(gene_name = gene_name, sample, lcpm) %>%
  pivot_wider(
    names_from = sample, values_from = lcpm,
    id_cols = gene_name
  )
# filter gene no 0 (lcpm>=1) >=3
lcpm_MM_fil <- lcpm_MM_wide %>%
  # rowwise() %>%
  # filter(sum(c_across(-gene_name) > 1) >= 3) %>%
  pivot_longer(
    cols = c(starts_with("SRR")),
    names_to = "sample",
    values_to = "lcpm"
  )


p1 <- lcpm_MM %>%
  dplyr::mutate(
    gene_id = str_replace_all(gene_name, "_.+", ""),
    gene_name = str_replace_all(gene_name, ".+_", "")
  ) %>%
  # dplyr::filter(gene_name %in% c("ZNF141", "MPHOSPH8",
  #                                "ATF7IP", "SIRT1",
  #                                "TASOR", "MORC2", "PPHLN1", "SETDB1",
  #                                "TASOR2", "IRF2",
  #                                "DNMT1", "DNMT3A")) %>%
  # dplyr::filter(gene_name %in% c("KLF6", "KLF4")) %>%
  # dplyr::filter(gene_name %in% c("ADAR")) %>%
  # dplyr::filter(gene_name %in% c("NFKB1")) %>%
  # dplyr::filter(gene_name %in% c("IGSF9", "ZNF28", "ZNF93", "DTL", "CKAP2L", "ZNF141",
  #                                "ZNF382", "ZNF793")) %>%
  dplyr::filter(gene_name %in% c("DTL", "CKAP2L", "ZNF860", "ZNF107")) %>%
  dplyr::left_join(MM_meta_final, by = c("sample" = "sampleID")) %>%
  dplyr::mutate(subtype = ifelse(sample %in% na.omit(PR_subtypeIDs), "PR",
    "Others"
  )) %>%
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm, fill=subtype, group = 1))+
  ggplot(aes(x = subtype, y = lcpm)) +
  # geom_bar(stat = "identity", position = "dodge")+
  # geom_col() +
  # geom_line(color = "red", size = 1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, position = position_dodge(0.8)) +
  geom_jitter(position = position_jitter(width = 0.15), alpha = 0.5) +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    label = "p.signif", vjust = 1
  ) +
  labs(x = "", y = "Gene expression (lcpm)") +
  facet_wrap(~gene_name, scales = "free_x", nrow = 2) +
  ggpubr::theme_pubr(base_size = 9.6)

ggsave(paste0(output_dir, "/lcpm_IRGI_NFKB1_PR.pdf"),
  plot = p1,
  width = 10, height = 8, units = "cm", limitsize = FALSE
)

ggsave(paste0(output_dir, "/lcpm_DTL_CKAP2L_PR.pdf"),
  plot = p1,
  width = 10, height = 8, units = "cm", limitsize = FALSE
)

## explore ADAR:
p2 <- lcpm_MM %>%
  dplyr::mutate(
    gene_id = stringr::str_replace_all(gene_name, "_.+", ""),
    gene_name = stringr::str_replace_all(gene_name, ".+_", "")
  ) %>%
  dplyr::filter(sample %in% MM_meta_final$SRR_id) %>%
  dplyr::left_join(MM_meta_final, by = c("sample" = "SRR_id")) %>%
  # dplyr::filter(gene_name %in% c("KLF6", "KLF4")) %>%
  # dplyr::filter(gene_name %in% c("ADAR")) %>%
  # dplyr::filter(gene_name %in% c("NFKB1")) %>%
  dplyr::filter(gene_name %in% c("DTL", "CKAP2L", "ZNF860", "ZNF107")) %>%
  # dplyr::filter(gene_name %in% c("IGSF9", "ZNF28", "ZNF93", "DTL", "CKAP2L", "ZNF141",
  #                                "ZNF382", "ZNF793")) %>%
  distinct() %>%
  ggplot(aes(x = L1PA2, y = lcpm)) +
  # ggplot(aes(x = lcpm , y = L1PA2)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") + # correlation coefficient
  labs(
    title = "",
    x = "L1PA2",
    y = "lcpm"
  ) +
  facet_wrap(~gene_name, scales = "free_y") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(output_dir, "/lcpm_DTL_CKAP2L.pdf"),
  device = "pdf", plot = p2,
  width = 20, height = 16, units = "cm", limitsize = FALSE
)
