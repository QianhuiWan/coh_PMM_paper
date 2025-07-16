
# only use BM primary samples
# R2: only to polish figures compare to R1
# + update primate-specific KZFPs
# + consider removing low purity samples for HR analysis if necessary
# + mark active L1 binding KZFPs

.libPaths()
library(tidyverse)
library(TCGAbiolinks)
library(edgeR)
library(survival)
library(survminer)
library(biomaRt)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(Hmisc)
library(corrplot)
library(cluster)
library(dynamicTreeCut)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
library(plyranges)
options(scipen=999)

out_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/allgene_hazard"

MM_meta <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/clinical_data_BMprimary_addID.tsv")

# get TCGA count mat, match barcode to clinial data
raw_counts_df_BMprimary <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/raw_counts_geneCounts_MMRF_BMprimary.tsv") %>% 
  tibble::column_to_rownames(var = "gene_id")

## create DGElist
### get group info.
if (unique(MM_meta$barcode == colnames(raw_counts_df_BMprimary))) {
  groups <- MM_meta$L1_quantile_group
}

### get DGElist
dge_rawCounts <- DGEList(counts=raw_counts_df_BMprimary, 
                         # remove.zeros = TRUE,
                         group = groups)
### filter by expr
keep <- filterByExpr(dge_rawCounts, group=dge_rawCounts$samples$group)
dge_rawCounts <- dge_rawCounts[keep, , keep.lib.sizes = FALSE]

### normalize counts with TMM method
dge_rawCounts <- calcNormFactors(dge_rawCounts, method = "TMM")
normalized_lcpm <- cpm(dge_rawCounts, normalized.lib.sizes = TRUE, log = TRUE)

# other input needed
KZFPs <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/kzfps_hg38_all.tsv")
kzfp_ranges <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/kzfps_hg38_all_geneCoordinates.tsv") %>% 
  dplyr::mutate(seqnames = paste0("chr", seqnames))
primate_specific_KZFPs <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/KZFP_evolutional_primate_specific.tsv")
KZFP_cluster_chr19 <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/KZFP_clusters_coord_hg38_chr19_noStrand.tsv")
L1binding_KZFPs <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/all_overlaps_L1_binding_KZFPs_MM_R4.tsv")

L1binding_KZFPs <- L1binding_KZFPs %>%
  group_by(KZFP_name) %>%       # Group by KZFP_name
  filter(n() > 10) %>%
  ungroup()

MM_subtypes <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/mmrf_subtypes/clinical_data_BMprimary_addID_addSubtypes.tsv")

#############################################################################
# step1: get the hazard roatio from survival analysis for all genes in MMRF dataset
#############################################################################
# or use function:
# Make sure genes are in rows, samples in columns
# expr_mat: rows = genes, cols = samples (colnames match clinical_data$sample_id)
# replace value<0 to 0 in exprMat
normalized_lcpm[normalized_lcpm < 0] <- 0
# Transpose to long format
expr_df <- as.data.frame(t(normalized_lcpm))  # now: rows = samples, columns = genes
expr_df$barcode <- rownames(expr_df)
gene_map <- data.frame(
  original = colnames(expr_df),
  safe = make.names(colnames(expr_df)))
colnames(expr_df) <- gene_map$safe
# Merge with clinical data
merged_df <- inner_join(MM_meta, expr_df, by = "barcode")

# List of genes to run Cox on
gene_list <- setdiff(colnames(expr_df), "barcode")  # exclude ID column

# Function to fit Cox model for each gene
fit_cox <- function(gene) {
  #formula <- as.formula(paste("Surv(time, status) ~", gene, "+ age + gender"))
  formula <- as.formula(paste("Surv(days_to_last_known_disease_status, OS) ~", gene))
  model <- tryCatch(
    coxph(formula, data = merged_df),
    error = function(e) return(NULL)
  )
  if (!is.null(model)) {
    s <- summary(model)
    tibble(
      gene = gene,
      HR = s$coefficients[,"exp(coef)"],
      lower95 = s$conf.int[,"lower .95"],
      upper95 = s$conf.int[,"upper .95"],
      pval = s$coefficients[,"Pr(>|z|)"]
    )
  } else {
    tibble(gene = gene, HR = NA, lower95 = NA, upper95 = NA, pval = NA)
  }
}

# Apply to all genes
results <- map_dfr(gene_list, fit_cox)

results <- results %>%
  dplyr::left_join(gene_map, by = c("gene" = "safe")) %>%
  dplyr::select(original, everything())

results$FDR <- p.adjust(results$pval, method = "BH")


write_tsv(results, file = paste0(out_dir, "/allgene_HR_MMRF.tsv"))

results <- read_tsv(paste0(out_dir, "/allgene_HR_MMRF.tsv"))

sig_genes <- results %>% filter(FDR < 0.05)

#############################################################################
# step 2: plot volcano plot: hazard roatio against -log2FDR
#############################################################################

results_forPlot <- results %>%
  dplyr::filter(!is.infinite(HR)) %>% 
  dplyr::filter(!is.na(HR)) %>%
  # dplyr::filter(original != "IFIT1P1_ENSG00000215515.2") %>% 
  dplyr::filter(HR <= 5) %>% 
  dplyr::mutate(gene_id_version = str_replace_all(original, ".+_", "")) %>% 
  dplyr::mutate(gene_id = str_remove(gene_id_version, "\\..*")) %>% 
  dplyr::mutate(gene_symbol = str_replace_all(original, "_.+", "")) %>% 
  mutate(isKZFP = ifelse(gene_symbol %in% KZFPs$hgnc_symbol, "Y", "N")) %>% 
  mutate(isL1KZFP = ifelse(gene_symbol %in% L1binding_KZFPs$KZFP_name, "Y", "N")) %>% 
  mutate(
    log2HR = log2(HR),
    neg_log2_FDR = -log2(FDR),
    sig = case_when(
      FDR < 0.05 & HR > 1 & isL1KZFP=="Y" ~ "Poor prognosis L1 KZFPs (HR>1)",
      FDR < 0.05 & HR < 1 & isL1KZFP=="Y" ~ "Good prognosis L1 KZFPs (HR<1)",
      TRUE ~ "Other genes"
    )
  ) 

top5_genes <- results_forPlot %>%
  dplyr::filter(sig != "Other genes") %>% 
  # arrange(desc(neg_log2_FDR)) %>%
  arrange(desc(HR), desc(neg_log2_FDR)) %>%
  slice_head(n = 11)

# tail5_genes <- results_forPlot %>%
#   dplyr::filter(sig != "Other genes") %>% 
#   # arrange(desc(neg_log2_FDR)) %>%
#   arrange(HR, desc(neg_log2_FDR)) %>%
#   slice_head(n = 6)
top_genes <- rbind(top5_genes)

p1 <- ggplot(results_forPlot, aes(x = HR, y = neg_log2_FDR, color = sig)) +
  # Plot gray dots first
  geom_point(data = subset(results_forPlot, sig == "Other genes"),
             size = 1, alpha = 0.6) +
  # Then plot red and blue dots on top
  geom_point(data = subset(results_forPlot, sig != "Other genes"),
             size = 3.5, alpha = 0.9) +
  scale_color_manual(
    values = c(
      "Poor prognosis L1 KZFPs (HR>1)" = "red", 
      "Good prognosis L1 KZFPs (HR<1)" = "blue", 
      "Other genes" = "grey80"),
    labels = c(
      "Poor prognosis L1 KZFPs (HR>1)" = "Poor prognosis\n L1 binding KZFPs (n=2)",
      "Good prognosis L1 KZFPs (HR<1)" = "Good prognosis\n L1 binding KZFPs (n=3)",
      "Other genes" = "Other genes")) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  theme_pubr() +
  labs(#title = "Gene-wise hazard ratios",
    x = "Hazard ratio",
    y = "-log2(FDR)",
    color = "Genes") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log2(0.05), linetype = "dashed", color = "black")+
  geom_text_repel(
    data = top_genes,
    aes(label = gene_symbol),  # change to your gene name column
    size = 4,
    box.padding = 0.3,
    max.overlaps = Inf)

ggsave(paste0(out_dir, 
              "/allgene_HR_R3.pdf"), 
       plot = p1,
       width = 12, height = 15, units = "cm", limitsize = FALSE)


########################enrichment##############################################
# step 3: for sig prognostic genes, plot the percentage of KZFP, TF and protein coding genes (PCG)
################################################################################
# biomart gene biotype annotation
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 

# for all genes
all_gene_id_version <- str_replace_all(results$original, ".+_", "")
all_genes_PCG <-  getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", 
                 "description", "gene_biotype"),
  filters = "ensembl_gene_id_version",
  values = all_gene_id_version,
  mart = ensembl) %>% 
  dplyr::filter(gene_biotype =="protein_coding")


sig_genes_id_version <- str_replace_all(sig_genes$original, ".+_", "")
sig_genes_PCG <-  getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", 
                 "description", "gene_biotype"),
  filters = "ensembl_gene_id_version",
  values = sig_genes_id_version,
  mart = ensembl) %>% 
  dplyr::filter(gene_biotype =="protein_coding")

# GO term for DNA-binding TFs
go_tf <- "GO:0003700"  # DNA-binding transcription factor activity

tf_genes <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_gene_id_version",
                 "external_gene_name", 
                 "description", "gene_biotype"),
  filters = "go",
  values = "GO:0003700",
  mart = ensembl
)

sig_genes_forPlot <- sig_genes %>% 
  dplyr::mutate(gene_id_version = str_replace_all(original, ".+_", "")) %>% 
  dplyr::mutate(gene_id = str_remove(gene_id_version, "\\..*")) %>% 
  dplyr::mutate(gene_symbol = str_replace_all(original, "_.+", "")) %>% 
  dplyr::mutate(Prognosis = ifelse(HR > 1, "Poor prognosis", 
                                   "Good prognosis")) %>% 
  dplyr::mutate(primate_specific_KZFPs = 
                  ifelse(gene_symbol %in% primate_specific_KZFPs$KZFP, "P", "O")) %>% 
  dplyr::mutate(isKZFP = ifelse(gene_symbol %in% KZFPs$hgnc_symbol, "Y", "N")) %>% 
  mutate(isL1KZFP = ifelse(gene_symbol %in% L1binding_KZFPs$KZFP_name, "Y", "N")) %>% 
  dplyr::mutate(isTF = ifelse(gene_symbol %in% tf_genes$external_gene_name, "Y", "N")) %>% 
  dplyr::mutate(isPCG = ifelse(gene_symbol %in% sig_genes_PCG$external_gene_name, "Y", "N")) %>% 
  dplyr::mutate(
    group = case_when(
      isKZFP=="Y" ~ "KZFPs",
      isKZFP=="N" & isTF=="Y" ~ "TFs",
      isKZFP=="N" & isTF=="N" & isPCG=="Y" ~ "PCGs",
      TRUE ~ "Other genes"
    )
  )

# get those genes and then plot percentage
p2 <- sig_genes_forPlot %>% 
  dplyr::filter(group != "Other genes") %>%
  ggplot(aes(x = group, fill = Prognosis))+
  geom_bar(position = "fill")+
  # scale_fill_viridis_d(option = "turbo") + #, direction = -1
  scale_fill_manual(
    values = c(
      "Poor prognosis" = "darkred",  # warm soft red
      "Good prognosis" = "#457B9D"  # dusty blue
    )) +
  labs(x = "", y = "Proportion (%)")+
  theme_pubr()
# KZFPs = 43; TFs=86; PCGs=1704
ggsave(paste0(out_dir, 
              "/gene_HR_sig_percentage.pdf"), 
       plot = p2,
       width = 12, height = 10, units = "cm", limitsize = FALSE)



#########################step 4: sig prognostic KZFPs########################
# step 4: for sig prognostic KZFPs, plot expression (primate/hominoid vs older KZFP groups)
#############################################################################
# get primate/hominoid only KZFPs: n=463
sig_KZFPS <- sig_genes_forPlot %>% 
  dplyr::filter(group == "KZFPs")

sig_L1_KZFPS <- sig_genes_forPlot %>% 
  dplyr::filter(isL1KZFP == "Y")

# Merge expr_df with clinical data + MM_subtypes
# save naming convention for expr_df
MM_subtypes_fil <- MM_subtypes %>% 
  dplyr::filter(Reason_For_Collection == "Baseline")

merged_expr_df <- inner_join(MM_subtypes_fil, expr_df, by = "barcode") 

# filter for sample with RNA_subtype info.
# and filter sig KZFPs 
merged_expr_df_fil <- merged_expr_df %>% 
  dplyr::filter(!is.na(RNA_Subtype_Name)) %>% 
  dplyr::filter(RNA_Subtype_Name != "Low purity") %>% 
  dplyr::mutate(RNA_subtypes = ifelse(RNA_Subtype_Name == "PR", 
                                      "PR", "Other subtypes")) %>% 
  dplyr::select(barcode, RNA_subtypes, sig_L1_KZFPS$gene) %>% 
  pivot_longer(cols = -c(barcode, RNA_subtypes), 
               names_to = "geneID", values_to = "lcpm") %>% 
  dplyr::left_join(gene_map, by = c("geneID" = "safe")) %>%
  dplyr::select(original, everything()) %>% 
  left_join(sig_L1_KZFPS, by = "original")

merged_expr_df_fil_mean <- merged_expr_df_fil %>% 
  group_by(barcode, Prognosis, primate_specific_KZFPs) %>% 
  summarise(mean_lcpm = mean(lcpm), .groups = 'drop') %>% 
  left_join(MM_subtypes_fil, by = "barcode") %>% 
  dplyr::filter(!is.na(RNA_Subtype_Name)) %>% 
  dplyr::filter(RNA_Subtype_Name != "Low purity") %>% 
  dplyr::mutate(RNA_subtypes = ifelse(RNA_Subtype_Name == "PR", 
                                      "PR", "Other subtypes"))
  
# plot KZFP expression between PR and other groups
p3 <- merged_expr_df_fil_mean %>% 
  # merged_expr_df_fil %>% 
  ggplot(aes(x = RNA_subtypes, 
             # y = lcpm, 
             y = mean_lcpm)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, position = position_dodge(0.8)) +
  geom_jitter(position = position_jitter(width = 0.15), alpha = 0.5) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif"
  ) +
  labs(x ="", y="Mean expression (lcpm)")+
  facet_wrap(~Prognosis) +
  theme_pubr()

# all sig HR KZFPs are older, not primate-specific
ggsave(paste0(out_dir, 
              "/gene_HR_sig_L1_KZFPS_expre.pdf"), 
       plot = p3,
       width = 12, height = 10, units = "cm", limitsize = FALSE)

# plot KZFP expression between PR and other groups
p3_2 <- merged_expr_df_fil_mean %>%
  ggplot(aes(x = RNA_subtypes, y = mean_lcpm,
             fill = primate_specific_KZFPs)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.9,
    position = position_dodge(width = 0.8)) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8),
    alpha = 0.5,
    size = 1
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  labs(x = "", y = "Mean expression (lcpm)") +
  facet_wrap(~Prognosis) +
  theme_pubr()

# all sig HR KZFPs are older, not primate-specific
ggsave(paste0(out_dir, 
              "/gene_HR_sig_L1_KZFPS_expre_OandP.pdf"), 
       plot = p3_2,
       width = 12, height = 10, units = "cm", limitsize = FALSE)

write_tsv(sig_L1_KZFPS, file = paste0(out_dir, "/gene_HR_MMRF_sig_L1_KZFPS.tsv"))

#######################step 5: cor##############################################
# step 5: plot correlation for these sig prognostic KZFPs
################################################################################
# sig_L1_KZFPS %>% head
# poor_kzfp_genes <- sig_L1_KZFPS[sig_L1_KZFPS$Prognosis=="Poor prognosis",]

expr_mat <- normalized_lcpm[sig_L1_KZFPS$original, ]
expr_for_corr <- t(expr_mat)  # genes in columns
colnames(expr_for_corr) <- str_replace_all(colnames(expr_for_corr), "_..+", "")
corr_res <- rcorr(expr_for_corr, type = "pearson")

corrplot(corr_res$r,
         method = "color",
         type = "upper",
         order = "hclust",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         p.mat = corr_res$P,
         sig.level = 0.05,       # Only show significant correlations
         insig = "blank",        # Leave blanks for non-sig
         tl.cex = 0.8)

# Extract r and p-value matrices
r_mat <- corr_res$r
p_mat <- corr_res$P
mean_expr <- colMeans(expr_for_corr, na.rm = TRUE)
chr_info <- kzfp_ranges %>% 
  dplyr::filter(gene_name %in% colnames(r_mat)) %>% 
  dplyr::arrange(match(gene_name, rownames(r_mat))) %>% 
  dplyr::filter(!seqnames%in% c("chrHSCHR19_1_CTG2", "chrHG2066_PATCH"))
  
# Mask insignificant correlations (set them to 0 or NA)
r_mat_masked <- r_mat
r_mat_masked[p_mat > 0.05] <- 0  # or NA

# add prognosis bar, add OP bar and add chr19 cluster bar for KZFP annotation
kzfp_info <- sig_L1_KZFPS %>% 
  dplyr::filter(gene_symbol %in% colnames(r_mat)) %>% 
  dplyr::arrange(match(gene_symbol, rownames(r_mat)))
# for chr19 kzfp clusters
KZFP_cluster_chr19_GR <- KZFP_cluster_chr19 %>% as_granges()
sig_kzfp_ranges_GR <- chr_info %>% as_granges()
cluster_info <- sig_kzfp_ranges_GR %>% 
  join_overlap_left(KZFP_cluster_chr19_GR) 
  # mutate(clusterID = str_replace_na(clusterID))

set.seed(166) 
ha = HeatmapAnnotation(
  Mean_lcpm = anno_barplot(mean_expr, gp = gpar(fill = "steelblue")),
  Prognosis = kzfp_info$Prognosis,
  Primate_specific = kzfp_info$primate_specific_KZFPs,
  Chromosome = chr_info$seqnames,
  Chr19_Cluster = cluster_info$clusterID,
  
  col = list(
    Prognosis = c("Good prognosis" = "#457B9D", "Poor prognosis" = "darkred"),
    Primate_specific = c("P" = "#00BFC4", "O" = "#F8766D"),
    Chromosome = structure(
      circlize::rand_color(length(unique(chr_info$seqnames))),
      names = unique(chr_info$seqnames)
    ),
    Chr19_Cluster = structure(
      circlize::rand_color(length(na.omit(unique(cluster_info$clusterID)))),
      names = na.omit(unique(cluster_info$clusterID))
    )
  ),
  annotation_name_side = "right"
  # show_legend = c(TRUE, FALSE)  # Hide Chromosome legend here
)
# Draw heatmap and place the chromosome legend on the right
# draw(ht, annotation_legend_list = list(lgd_chr), annotation_legend_side = "right")



# Plot Heatmap
# p4
pdf(paste0(out_dir, 
           "/gene_HR_sig_L1_KZFPS_correlation_heatmap.pdf"), 
    width = 10, height = 8)
Heatmap(
  r_mat_masked,
  name = "Pearson\ncor",
  top_annotation = ha,
  col = colorRamp2(
    c(-1, 0, 1),
    c("#2166AC", "white", "#B2182B")  # Blue → White → Red
  ),
  show_row_dend = FALSE,
  show_column_dend = FALSE, 
  show_row_names = TRUE, row_names_side = "left",
  clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson",
  column_title = "Gene-Gene Correlation",
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  row_title = NULL,
  heatmap_legend_param = list(title_position = "topcenter")
)
dev.off()



########################step6:##################################################
# step6: plot KZFP cluster sample group survival plot
################################################################################
## separate MM samples into 2 groups according to sig prognostic KZFPs
sig_L1_KZFPS <- sig_L1_KZFPS$original
expr_sig <- normalized_lcpm[sig_L1_KZFPS, ]  # genes × samples
pca <- prcomp(t(expr_sig), scale. = TRUE)
sample_score <- pca$x[, 2]  # use PC1

# risk_group <- ifelse(sample_score > median(sample_score), "Low", "High")

risk_group <- ifelse(sample_score > median(sample_score), "High", "Low")

library(mclust)
expr_scaled <- scale(expr_sig)
mc <- Mclust(t(expr_scaled), G = 2)
clusters <- mc$classification
plot(mc, what = "classification")
plot(mc, what = "uncertainty")

risk_group <- ifelse(clusters ==2, "Low", "High")

## plot the survival plot
# Create survival input data
surv_df <- MM_subtypes_fil %>%
  mutate(risk_group = risk_group[match(barcode, colnames(expr_sig))])

# Make the Surv object
surv_obj <- Surv(time = surv_df$days_to_last_known_disease_status, 
                 event = surv_df$OS)

# Fit the model
fit <- survfit(surv_obj ~ risk_group, data = surv_df)

# Plot it
p5 <- ggsurvplot(
  fit,
  data = surv_df,
  pval = TRUE,
  # risk.table = TRUE,
  palette = c("red", "blue"),
  title = "Survival by KZFP-based Risk Group",
  legend.title = "Risk Group",
  xlab = "Time (days)",  
  surv.median.line = "hv"
)

ggsave(paste0(out_dir, 
              "/gene_HR_sig_L1binding_KZFPS_2groups_surv_pc1pc2.pdf"), 
       plot = p5$plot,
       width = 18, height = 16, units = "cm", limitsize = FALSE)


# set.seed(52)
# km <- kmeans(t(expr_scaled), centers = 2)  # or however many clusters you want
# clusters <- km$cluster  # numeric vector, one value per sample

dist_mat <- dist(t(expr_scaled), method = "euclidean")  # or "manhattan", etc.
# Hierarchical clustering
hc <- hclust(dist_mat, method = "ward.D2")  # or "average", "complete", "single"
clusters <- cutreeDynamic(dendro = hc, distM = as.matrix(dist_mat), deepSplit = 2)

library(mclust)

mc <- Mclust(t(expr_scaled), G = 2)
clusters <- mc$classification


# Example: assign one color per cluster
cluster_colors <- c("red", "blue")
col_vector <- cluster_colors[clusters]

pdf(paste0(out_dir, "/pca_PC1andPC2_forL1bindingKZFPs.pdf"), , width = 7, height = 6)
plot(pca$x[,1:2],
     col = col_vector,
     pch = 19,
     xlab = paste0("PC1 (", round(100 * summary(pca)$importance[2,1], 1), "%)"),
     ylab = paste0("PC2 (", round(100 * summary(pca)$importance[2,2], 1), "%)"),
     main = "PCA: Samples Colored by Cluster")
legend("topright",
       legend = paste("Cluster", sort(unique(clusters))),
       col = cluster_colors,
       pch = 19)
dev.off()





