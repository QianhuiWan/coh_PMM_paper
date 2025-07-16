
# use count matrix from feature counts
# filterByExpr, TMM normalize
# get gene counts stats, MMRF, use limma voom method
# and then get DE topTable

dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2"
dir_gemini = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2_gemini"
# output_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR2"
output_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_redo_res"
# input_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s1_checkStrandness_R4"
# dir_gemini = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s1_checkStrandness_R4_gemini"

# load packages
library(tidyverse)
library(magrittr)
library(stringr)
library(rtracklayer)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(edgeR)
library(biomaRt)

library(clusterProfiler)
library(msigdbr)
library(enrichR)
library(enrichplot)
library(ggplot2)

# load sample sheet: sample in 4 groups by activated L1
patientsIn4Groups <- 
  read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_redo_res/all_activeL1PA2_patientsIn4Groups_764.tsv")

# MM_metadata <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/clinical_data_BMprimary_addID.tsv")
MM_metadata_update <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/clinical_data_BMprimary_addID_addSubtypes.tsv")

# only get Q1 and Q4 count mat data


## use getNum FUN to get all transcripts #######################################
# Save as a .tsv file
# write.table(count_matrix, file =  paste0(dir, "/gene_count_matrix.tsv"),
#             sep = "\t", row.names = TRUE, col.names = TRUE)
# load count matrix
count_matrix <- 
  read.table(paste0(dir, "/gene_count_matrix_STARcounts_MM_764.tsv"), 
             sep = "\t", header = TRUE, row.names = 1) %>% 
  as.matrix()

# write.table(count_matrix, file =  paste0(dir, "/gene_count_matrix_STARcounts.tsv"), 
#             sep = "\t", row.names = TRUE, col.names = TRUE)

# DE
## match meta with matrix column names (i.e. sample IDs)
# count_matrix <- 
#   count_matrix[,colnames(count_matrix)%in%patientsIn4Groups$sampleID]
# count_matrix <- 
#   count_matrix[, match(patientsIn4Groups$sampleID, colnames(count_matrix))]

dge <- DGEList(counts = count_matrix, 
               group =  patientsIn4Groups$L1_quantile_group)

# Define a threshold for filtering
# Keep genes with CPM > 1 in at least 1 group
keep <- rowSums(cpm(dge, normalized.lib.sizes = TRUE) >=1 ) >= min(table(patientsIn4Groups$L1_quantile_group))
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


## design matrix
design <- model.matrix(~ 0+L1_quantile_group, data=patientsIn4Groups)
## make contrast
con <- makeContrasts(group = L1_quantile_groupQ4-L1_quantile_groupQ1,
                     levels=design) # TumorS - TumorM
v <- 
  voom(dge, design) %>%
  lmFit(design) %>%
  contrasts.fit(con) %>%
  eBayes()
dt <- decideTests(v)
print(summary(dt))
plotSA(v, main="Final model: Mean-variance trend")

tfit <- treat(v, lfc=1)
dt_tfit <- decideTests(tfit)
summary(dt_tfit)

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))

## get topTable output
dat <- topTable(v, coef = "group", number=Inf) %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  dplyr::mutate(gene_id_version = str_replace_all(gene_name, "_.+", "")) %>% 
  dplyr::mutate(gene_id = str_remove(gene_id_version, "\\..*")) %>% 
  dplyr::mutate(gene_symbol = str_replace_all(gene_name, ".+_", "")) %>% 
  # distinct(gene_name, .keep_all = TRUE) %>%
  dplyr::select(-t, -B, FDR = `adj.P.Val`) %>%
  mutate(Sig = ifelse(FDR < 0.05 & logFC > 2, "Up", 
                      ifelse(FDR < 0.05 & logFC < -2, "Down", "NS"))) 

# get geneID from geneID version
library(biomaRt)
## load mart: #i.e. Human genes (GRCh38.p14)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
## get gene ids
valid_genes <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id",
                                    "ensembl_gene_id_version"),
                     filters = "ensembl_gene_id_version",
                     values = dat$gene_id_version,
                     mart = ensembl)
# update dat
dat <- 
  dat %>% 
  dplyr::left_join(valid_genes, 
                   by = c("gene_id_version"="ensembl_gene_id_version"))

write_tsv(dat, file = paste0(output_dir, "/lmVoom_topTable_Q4Q1_dat.tsv"))


## GSEA
library(msigdbr)
library(clusterProfiler)
as.data.frame(msigdbr_collections())
m_df = msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- m_df[,c("gs_name", "db_ensembl_gene")] # "human_gene_symbol",

### Create a ranked list
ranked_list <- dat %>%
  dplyr::select(gene_id, logFC) %>% 
  arrange(desc(logFC)) %>%
  deframe()  # Convert to named vector (gene names as names, logFC as values)

## rank result by p values
dum <- dat$FDR
#Avoid Inf when taking log of 0 
dum[which(dum==0)] <- 1e-300
score <- -log10(dum)
#It is important to differentiate up and down regulation
score <- ifelse(dat$logFC>=0, score, -1*score)
names(score) <- dat$gene_id
score <- sort(score, decreasing = T)
set.seed(51)
score <- score + rnorm(length(score), mean = 0, sd = 1e-6)
score <- sort(score, decreasing = TRUE)

em_1 <- GSEA(ranked_list, TERM2GENE=hallmark, eps = 0,
           pvalueCutoff=0.25, pAdjustMethod="BH")

em_2 <- GSEA(score, TERM2GENE=hallmark,
           pvalueCutoff=0.25, pAdjustMethod="BH")


# Create a dotplot of top 30 pathways
dotplot(em_1, orderBy="NES", x="NES", showCategory = 20)
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_NES_top20.pdf"), width = 9, height = 12)

dotplot(em_2, orderBy="NES", x="NES", showCategory = 20)
ggsave(paste0(output_dir, "/lmVoomPval_gsea_NES_top20.pdf"), width = 9, height = 12)

dotplot(em_1, orderBy="NES", x="NES", showCategory = 25)
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_NES_top25.pdf"), width = 9, height = 12)

dotplot(em_2, orderBy="NES", x="NES", showCategory = 25)
ggsave(paste0(output_dir, "/lmVoomPval_gsea_NES_top25.pdf"), width = 9, height = 12)

# plot GSEA plot for key pathways
gseaplot(em, geneSetID="HALLMARK_P53_PATHWAY", title="HALLMARK_P53_PATHWAY")
gseaplot2(em, geneSetID="HALLMARK_P53_PATHWAY", title="HALLMARK_P53_PATHWAY")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_P53_PATHWAY.pdf"), width = 8, height = 6)

gseaplot(em, geneSetID="HALLMARK_E2F_TARGETS", title="HALLMARK_E2F_TARGETS")
gseaplot2(em, geneSetID="HALLMARK_E2F_TARGETS", title="HALLMARK_E2F_TARGETS")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_E2F_TARGETS.pdf"), width = 8, height = 6)

gseaplot(em, geneSetID="HALLMARK_G2M_CHECKPOINT", title="HALLMARK_G2M_CHECKPOINT")
gseaplot2(em, geneSetID="HALLMARK_G2M_CHECKPOINT", title="HALLMARK_G2M_CHECKPOINT")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_G2M_CHECKPOINT.pdf"), width = 8, height = 6)


gseaplot(em, geneSetID="HALLMARK_INTERFERON_ALPHA_RESPONSE", title="HALLMARK_INTERFERON_ALPHA_RESPONSE")
gseaplot2(em, geneSetID="HALLMARK_INTERFERON_ALPHA_RESPONSE", title="HALLMARK_INTERFERON_ALPHA_RESPONSE")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_INTERFERON_ALPHA_RESPONSE.pdf"), width = 8, height = 6)

gseaplot(em, geneSetID="HALLMARK_INTERFERON_GAMMA_RESPONSE", title="HALLMARK_INTERFERON_GAMMA_RESPONSE")
gseaplot2(em, geneSetID="HALLMARK_INTERFERON_GAMMA_RESPONSE", title="HALLMARK_INTERFERON_GAMMA_RESPONSE")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_INTERFERON_GAMMA_RESPONSE.pdf"), width = 8, height = 6)

gseaplot(em, geneSetID="HALLMARK_TNFA_SIGNALING_VIA_NFKB", title="HALLMARK_TNFA_SIGNALING_VIA_NFKB")
gseaplot2(em, geneSetID="HALLMARK_TNFA_SIGNALING_VIA_NFKB", title="HALLMARK_TNFA_SIGNALING_VIA_NFKB")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_TNFA_SIGNALING_VIA_NFKB.pdf"), width = 8, height = 6)







