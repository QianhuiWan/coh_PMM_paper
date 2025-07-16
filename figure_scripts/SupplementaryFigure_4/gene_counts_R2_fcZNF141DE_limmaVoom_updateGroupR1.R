
# use count matrix from feature counts
# filterByExpr, TMM normalize
# get gene counts stats, MMRF, use limma voom method
# and then get DE topTable

dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2"
dir_gemini = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2_gemini"
output_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR1"
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
  read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R2/all_activeL1s_patientsIn4Groups.tsv")

patientsIn4Groups <- 
  read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/all_activeL1PA2_patientsIn4Groups.tsv")

ZNF141_info <- read_tsv(paste0(output_dir, "/ZNF141_info.tsv"))

## use getNum FUN to get all transcripts #######################################
# Save as a .tsv file
# write.table(count_matrix, file =  paste0(dir, "/gene_count_matrix.tsv"),
#             sep = "\t", row.names = TRUE, col.names = TRUE)
# load count matrix
count_matrix <- 
  read.table(paste0(dir, "/gene_count_matrix.tsv"), 
             sep = "\t", header = TRUE, row.names = 1) %>% 
  as.matrix()

# write.table(count_matrix, file =  paste0(dir, "/gene_count_matrix_STARcounts.tsv"), 
#             sep = "\t", row.names = TRUE, col.names = TRUE)

# DE
## match meta with matrix column names (i.e. sample IDs)
count_matrix <- 
  count_matrix[,colnames(count_matrix)%in%ZNF141_info$sampleID]
count_matrix <- 
  count_matrix[, match(ZNF141_info$sampleID, colnames(count_matrix))]

dge <- DGEList(counts = count_matrix, 
               group = ZNF141_info$ZNF141group)
# Define a threshold for filtering
# Keep genes with CPM > 1 in at least 1 group
# keep <- rowSums(cpm(dge) >=1 ) >= floor(870/4)
# keep <- rowSums(cpm(dge) > 2) >= 2
# keep <- rowSums(cpm(dge) >= 2) >= floor(870*0.75)

keep <- filterByExpr(dge, group=dge$samples$group)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge, method = "TMM")

# normalize and filter
normalized_lcpm <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)
# save filtered normalized lcpm for later plotting
# write.table(normalized_lcpm, file =  paste0(output_dir, "/lcpm_TMM_fil_ZNFgroups.tsv"),
#             sep = "\t", row.names = TRUE, col.names = TRUE)


# plot MDS
# lcpm <- cpm(count_matrix, log=TRUE)
# plotMDS(lcpm, top = nrow(lcpm), gene.selection = 'common', 
#         labels=colnames(count_matrix))

lcpm <- cpm(count_matrix, log=TRUE)
library(RColorBrewer)
col.group <- ZNF141_info$ZNF141group %>% 
  factor(levels = unique(ZNF141_info$ZNF141group))
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
par(mfrow = c(1, 2), mar = c(2, 2, 2, 1))  # Reduce margins
plotMDS(lcpm, top = nrow(lcpm), gene.selection = 'common', 
        labels=colnames(count_matrix), col=as.character(col.group))
legend("topright", col=levels(col.group), pch=c(15,16), 
       legend=c("Q1","Q2", "Q3", "Q4"))

plotMDS(normalized_lcpm, top = nrow(normalized_lcpm), gene.selection = 'common', 
        labels=colnames(count_matrix), col=as.character(col.group))
legend("topright", col=levels(col.group), pch=c(15,16), 
       legend=c("Q1","Q2", "Q3", "Q4"))

# Set up the plotting area with 5 rows
par(mfrow = c(6, 1), mar = c(2, 2, 2, 1))  # Reduce margins
boxplot(normalized_lcpm[,1:100], las=2, main="")
boxplot(normalized_lcpm[,101:200], las=2, main="")
boxplot(normalized_lcpm[,201:300], las=2, main="")
boxplot(normalized_lcpm[,601:700], las=2, main="")
boxplot(normalized_lcpm[,701:800], las=2, main="")
boxplot(normalized_lcpm[,801:870], las=2, main="")
title(main="Normalised data",ylab="Log-cpm")

## design matrix
design <- model.matrix(~0 + ZNF141group, data=ZNF141_info)
## make contrast
con <- makeContrasts(group = ZNF141groupnoZNF141-ZNF141grouphighZNF141,
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

write_tsv(dat, file = paste0(output_dir, "/lmVoom_topTable_NOvsHighZNF141_dat.tsv"))
# write_tsv(dat, file = paste0(output_dir, "/lmVoom_topTable_HighvsNOZNF141_dat.tsv"))

## GSEA
as.data.frame(msigdbr_collections())
m_df = msigdbr(species = "Homo sapiens", category = "H")
hallmark <- m_df[,c("gs_name", "human_ensembl_gene")] # "human_gene_symbol",

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

## perform GSEA
em <- GSEA(ranked_list, TERM2GENE=hallmark, eps = 0,
           pvalueCutoff=0.25, pAdjustMethod="BH")

em2 <- GSEA(score, TERM2GENE=hallmark,
           pvalueCutoff=0.25, pAdjustMethod="BH")


# Create a dotplot of top 30 pathways
dotplot(em, orderBy="NES", x="NES", showCategory = 10)
ggsave(paste0(output_dir, "/lmVoomLfc_ZNF141group_gsea_NES_top10.pdf"), width = 12, height = 5)

dotplot(em2, orderBy="NES", x="NES", showCategory = 10)
ggsave(paste0(output_dir, "/lmVoomFDR_ZNF141group_gsea_NES_top10.pdf"), width = 10, height = 8)

# plot GSEA plot for key pathways
gseaplot(em, geneSetID="HALLMARK_P53_PATHWAY", title="HALLMARK_P53_PATHWAY")
gseaplot2(em, geneSetID="HALLMARK_P53_PATHWAY", title="HALLMARK_P53_PATHWAY")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_P53_PATHWAY.pdf"), width = 8, height = 6)

gseaplot(em, geneSetID="HALLMARK_E2F_TARGETS", title="HALLMARK_E2F_TARGETS")
gseaplot2(em, geneSetID="HALLMARK_E2F_TARGETS", title="HALLMARK_E2F_TARGETS")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_E2F_TARGETS.pdf"), width = 8, height = 6)

gseaplot(em, geneSetID="HALLMARK_OXIDATIVE_PHOSPHORYLATION", title="HALLMARK_OXIDATIVE_PHOSPHORYLATION")
gseaplot2(em, geneSetID="HALLMARK_OXIDATIVE_PHOSPHORYLATION", title="HALLMARK_OXIDATIVE_PHOSPHORYLATION")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_G2M_CHECKPOINT.pdf"), width = 8, height = 6)

gseaplot(em2, geneSetID="HALLMARK_OXIDATIVE_PHOSPHORYLATION", title="HALLMARK_OXIDATIVE_PHOSPHORYLATION")
gseaplot2(em2, geneSetID="HALLMARK_OXIDATIVE_PHOSPHORYLATION", title="HALLMARK_OXIDATIVE_PHOSPHORYLATION")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_G2M_CHECKPOINT.pdf"), width = 8, height = 6)

gseaplot(em, geneSetID="HALLMARK_INTERFERON_ALPHA_RESPONSE", title="HALLMARK_INTERFERON_ALPHA_RESPONSE")
gseaplot2(em, geneSetID="HALLMARK_INTERFERON_ALPHA_RESPONSE", title="HALLMARK_INTERFERON_ALPHA_RESPONSE")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_INTERFERON_ALPHA_RESPONSE.pdf"), width = 8, height = 6)

gseaplot(em, geneSetID="HALLMARK_INTERFERON_GAMMA_RESPONSE", title="HALLMARK_INTERFERON_GAMMA_RESPONSE")
gseaplot2(em, geneSetID="HALLMARK_INTERFERON_GAMMA_RESPONSE", title="HALLMARK_INTERFERON_GAMMA_RESPONSE")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_INTERFERON_GAMMA_RESPONSE.pdf"), width = 8, height = 6)

gseaplot(em2, geneSetID="HALLMARK_MYC_TARGETS_V1", title="HALLMARK_MYC_TARGETS_V1")
gseaplot2(em2, geneSetID="HALLMARK_MYC_TARGETS_V1", title="HALLMARK_MYC_TARGETS_V1")
ggsave(paste0(output_dir, "/lmVoomLfc_gsea_HALLMARK_TNFA_SIGNALING_VIA_NFKB.pdf"), width = 8, height = 6)







