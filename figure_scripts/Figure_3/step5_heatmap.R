# use count matrix from feature counts
# get gene counts stats, MMRF
# get DE topTable
pathway_genes_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2"
# dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR1"
dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3_uniqueMap/filtered_L1_redo_res"
# laod packages
library(tidyverse)
library(magrittr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(ComplexHeatmap)
library(circlize)
library(biomaRt)
# Load required packages
library(dplyr)
library(stringr)
library(purrr)
library(tibble)


# load sample sheet: sample in 4 groups by activated L1
patientsIn4Groups <-
  read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3_uniqueMap/filtered_L1_redo_res/all_activeL1PA2_patientsIn4Groups_764.tsv")

# patientsIn4Groups <- 
#   read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3_uniqueMap/filtered_L1_redo_res/all_activeL1s_patientsIn4Groups_764.tsv")

MM_metadata_update <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/clinical_data_BMprimary_addID_addSubtypes_uniqueMap.tsv")
MM_metadata_update <- MM_metadata_update[match(patientsIn4Groups$sampleID, MM_metadata_update$SRR_id),]
table(patientsIn4Groups$sampleID == MM_metadata_update$SRR_id)

patientsIn4Groups$low_purity <- MM_metadata_update$`Low purity`
patientsIn4Groups$sex <- MM_metadata_update$gender
patientsIn4Groups$age <- MM_metadata_update$age_at_index

# patientsIn4Groups <- patientsIn4Groups %>% 
#   dplyr::filter(sampleID %in% MM_metadata$sampleID) %>% 
#   arrange(match(sampleID, MM_metadata$sampleID))

# load topTable output
# topTable_dat <- read_tsv(paste0(dir, "/lmVoom_topTable_Q4Q1_dat.tsv")) %>% 
#   dplyr::select(gene_id, logFC, FDR)
# topTable_dat <- read_tsv(paste0(dir, "/lmVoom_topTable_Q4Q1_dat_purityAdj.tsv")) %>% 
#   dplyr::select(gene_id, logFC, FDR)
topTable_dat <- read_tsv(paste0(dir, "/lmVoom_topTable_L1groupsQ4Q1_dat_purityAdj_r2.tsv")) %>% 
  dplyr::select(gene_id, logFC, FDR)

# lcpm <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s1_checkStrandness_R4_gemini_uniqueMap/lcpm_TMM_764_STARcounts_cpmFilterApplied.tsv")
# lcpm <- read_tsv(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3_uniqueMap/filtered_L1_redo_res/lcpm_TMM_MM_STARcounts_l1groups_bc.tsv")

lcpm <- read_tsv(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3_uniqueMap/filtered_L1_redo_res/lcpm_TMM_MM_STARcounts.tsv")

MM_metadata_update <- MM_metadata_update[MM_metadata_update$SRR_id %in% lcpm$sample,]
# lcpm <- lcpm[,colnames(lcpm)%in%MM_metadata$sampleID]
# lcpm <- lcpm[,match(colnames(lcpm), MM_metadata$sampleID)]
# lcpm <- lcpm %>%
#   rownames_to_column(var="gene_name")

# load TMM lcpm matrix
library(dplyr)
normalized_lcpm_fil_DEtb <- lcpm %>% 
  # read.table(paste0(dir, "/lcpm_TMM_fil.tsv"), 
  #            sep = "\t", header = TRUE, row.names = 1) |> 
  # as_tibble(rownames="gene_name") |> 
  # dplyr::select(colnames(.) %in% patientsIn4Groups$sampleID) %>% 
  # dplyr::select(match(colnames(.), patientsIn4Groups$sampleID)) %>%
  dplyr::mutate(gene_id_version = str_replace_all(gene_name, "_.+", "")) %>% 
  dplyr::mutate(gene_id = str_remove(gene_id_version, "\\..*")) %>% 
  dplyr::mutate(gene_symbol = str_replace_all(gene_name, ".+_", "")) %>% 
  # dplyr::filter(!(gene_id == "ENSG00000228314")) %>%  # redundant Pseudogene id
  dplyr::left_join(topTable_dat, by = "gene_id") %>% 
  dplyr::filter(FDR<0.05 & abs(logFC)>1) %>%
  # dplyr::select(-FDR, -logFC) %>%
  # mutate_at(vars(-gene_name, -gene_id_version, 
  #                -gene_id, -gene_symbol), scale) |>
  # pivot_longer(cols = -c(gene_name, gene_id_version, gene_id, gene_symbol), 
  #              names_to = "sample", values_to = "lcpm") %>% 
  dplyr::left_join(MM_metadata_update, by = c("sample"="sampleID")) %>% 
  dplyr::mutate(lcpm = ifelse(lcpm < 0, 0, lcpm))


DE_up_genes <- normalized_lcpm_fil_DEtb %>% dplyr::filter(logFC>0) 
DE_down_genes <- normalized_lcpm_fil_DEtb %>% dplyr::filter(logFC<0)


normalized_lcpm_fil_DE_mat <- normalized_lcpm_fil_DEtb %>% 
  pivot_wider(id_cols = gene_id, names_from = sample, 
              values_from = lcpm, values_fill = 0) %>% 
  column_to_rownames(var = "gene_id") %>% 
  as.matrix()

normalized_lcpm_fil_DE_mat <- normalized_lcpm_fil_DE_mat[,match(MM_metadata_update$sampleID, colnames(normalized_lcpm_fil_DE_mat))]

# # get protein coding DE genes
# # Ensembl mart
# ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# 
# # gene biotype
# gene_info <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"),
#                    filters = "ensembl_gene_id",
#                    values = rownames(normalized_lcpm_fil_DE_mat),
#                    mart = ensembl)
# # protein-coding
# protein_coding_DEgenes <- gene_info %>%
#   filter(gene_biotype == "protein_coding") %>%
#   pull(ensembl_gene_id)
# 
# normalized_lcpm_fil_DE_mat <- normalized_lcpm_fil_DE_mat[rownames(normalized_lcpm_fil_DE_mat) %in% protein_coding_DEgenes,
#                                                          match(MM_metadata_update$sampleID, colnames(normalized_lcpm_fil_DE_mat))]

# update gene mapping to pathways
library(msigdbr)
library(dplyr)

# Get full gene-to-pathway mapping (Reactome, KEGG, GO, Hallmark, etc.)
all_mappings <- msigdbr(species = "Homo sapiens", db_species = "HS", 
                        collection = "C2") #subcollection = "CP:KEGG_LEGACY"
# For simplicity, keep only the pathway name and gene symbol
gene_pathway_map <- all_mappings %>%
  dplyr::select(gene_symbol, gs_collection, gs_subcollection, gs_name, ensembl_gene)


# 1. Define expanded keyword match rules
keyword_patterns <- list(
  # NFKB = "NFKB|TNFA|INFLAMMATORY|IKB|REL",
  INTERFERON = "INTERFERON|IFN|ANTIVIRAL",
  P53 = "P53_PATHWAY|DNA_DAMAGE|CELL_CYCLE_ARREST",
  E2F = "E2F|CELL_CYCLE|G1_S_TRANSITION"
)

# 2. Build gene sets per keyword based on matched pathway names
keyword_gene_sets <- imap(keyword_patterns, function(pattern, keyword) {
  gene_pathway_map %>%
    filter(str_detect(gs_name, regex(pattern, ignore_case = TRUE))) %>%
    pull(ensembl_gene) %>%
    unique()
})

# 3. Get all DE gene IDs
# de_gene_ids <- rownames(normalized_lcpm_fil_DE_mat)
de_up_gene_ids <- DE_up_genes$gene_id %>% unique()
de_down_gene_ids <- DE_down_genes$gene_id %>% unique()

# 4. Assign each DE gene to top and second keyword by overlap in shared pathways
## up genes
de_up_gene_keyword_ranked <- map_dfr(de_up_gene_ids, function(gene) {
  gene_paths <- gene_pathway_map %>%
    filter(ensembl_gene == gene) %>%
    pull(gs_name)
  
  partner_genes <- gene_pathway_map %>%
    filter(gs_name %in% gene_paths) %>%
    pull(ensembl_gene)
  
  keyword_scores <- map_int(keyword_gene_sets, ~ sum(partner_genes %in% .x))
  keyword_scores_df <- tibble(keyword = names(keyword_scores), score = keyword_scores) %>%
    arrange(desc(score))
  
  tibble(
    de_gene_id = gene,
    top_keyword = ifelse(nrow(keyword_scores_df) >= 1 && keyword_scores_df$score[1] > 0, keyword_scores_df$keyword[1], "Unclassified"),
    top_score = keyword_scores_df$score[1],
    second_keyword = ifelse(nrow(keyword_scores_df) >= 2 && keyword_scores_df$score[2] > 0, keyword_scores_df$keyword[2], "Unclassified"),
    second_score = ifelse(nrow(keyword_scores_df) >= 2, keyword_scores_df$score[2], "Unclassified")
  )
})


## down genes
de_down_gene_keyword_ranked <- map_dfr(de_down_gene_ids, function(gene) {
  gene_paths <- gene_pathway_map %>%
    filter(ensembl_gene == gene) %>%
    pull(gs_name)
  
  partner_genes <- gene_pathway_map %>%
    filter(gs_name %in% gene_paths) %>%
    pull(ensembl_gene)
  
  keyword_scores <- map_int(keyword_gene_sets, ~ sum(partner_genes %in% .x))
  keyword_scores_df <- tibble(keyword = names(keyword_scores), score = keyword_scores) %>%
    arrange(desc(score))
  
  tibble(
    de_gene_id = gene,
    top_keyword = ifelse(nrow(keyword_scores_df) >= 1 && keyword_scores_df$score[1] > 0, keyword_scores_df$keyword[1], "Unclassified"),
    top_score = keyword_scores_df$score[1],
    second_keyword = ifelse(nrow(keyword_scores_df) >= 2 && keyword_scores_df$score[2] > 0, keyword_scores_df$keyword[2], "Unclassified"),
    second_score = ifelse(nrow(keyword_scores_df) >= 2, keyword_scores_df$score[2], "Unclassified")
  )
})

print(head(de_up_gene_keyword_ranked, 10))
table(de_up_gene_keyword_ranked$top_keyword)
table(de_up_gene_keyword_ranked$second_keyword)

print(head(de_down_gene_keyword_ranked, 10))
table(de_down_gene_keyword_ranked$top_keyword)
table(de_down_gene_keyword_ranked$second_keyword)

de_down_gene_keyword_ranked$second_keyword <- ifelse(
  de_down_gene_keyword_ranked$second_keyword == "E2F",
  de_down_gene_keyword_ranked$top_keyword,
  de_down_gene_keyword_ranked$second_keyword
)

# de_gene_keyword_final <- 
#   tibble(de_gene_id = c(de_up_gene_keyword_ranked$de_gene_id, 
#                         de_down_gene_keyword_ranked$de_gene_id),
#          Pathway = c(de_up_gene_keyword_ranked$second_keyword, 
#                      de_down_gene_keyword_ranked$top_keyword)) %>% 
#   arrange(match(de_gene_id, rownames(normalized_lcpm_fil_DE_mat)))

de_gene_keyword_final <- 
  tibble(de_gene_id = c(de_up_gene_keyword_ranked$de_gene_id, 
                        de_down_gene_keyword_ranked$de_gene_id),
         Pathway = c(de_up_gene_keyword_ranked$top_keyword, 
                     de_down_gene_keyword_ranked$second_keyword)) %>% 
  arrange(match(de_gene_id, rownames(normalized_lcpm_fil_DE_mat)))

# plot heatmap, use complex heatmap

## all DE genes heatmap
# palette_value <- circlize::colorRamp2(
#   breaks = c(-10, 0, 10, 30),  # Adjust range to match your data
#   colors = c("blue", "white", "brown", "red")  # Gradient: purple -> white -> yellow
# )
# palette_value <-  circlize::colorRamp2(
#   seq(-5, 5, length.out = 11), 
#   rev(RColorBrewer::brewer.pal(11, "RdBu")))

# Heatmap(normalized_lcpm_fil_DE_mat)


# key genes heatmap
# only add L1PA2 number
# keep <- rownames(normalized_lcpm_fil_DE_mat) %in% all_key_pathway_unique$human_ensembl_gene
# mat_key <- normalized_lcpm_fil_DE_mat[keep,]
mat_key <- normalized_lcpm_fil_DE_mat
# Scale by row (z-score scaling)
mat_key_scaled <- t(apply(mat_key, 1, scale))  # Apply scaling (z-score) across rows
colnames(mat_key_scaled) <- colnames(mat_key)


col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "brown"))
column_ha <- HeatmapAnnotation(
  L1PA2 = anno_barplot(MM_metadata_update$L1PA2,
                       gp = gpar(fill = "black")),
  height = unit(1.5, "cm"),
  na_col = "white"
)


row_ha = rowAnnotation(Pathway = de_gene_keyword_final$Pathway,
                       col = list(Pathway = c("E2F" = "pink", 
                                              "INTERFERON" = "lightblue", 
                                              "P53" = "darkgreen", 
                                              # "NFKB" = "orange", 
                                              "Unclassified"="white")), 
                       annotation_name_side = "top",          
                       annotation_name_gp = gpar(col = NA),
                       width = unit(0.5, "cm"))


# get gene symbol
# gene_id_symbol <- normalized_lcpm_fil_DEtb %>% 
#   dplyr::select(gene_id, gene_id_version, gene_symbol) %>% unique() %>% 
#   dplyr::filter(gene_id %in% rownames(mat_key_scaled)) %>% 
#   arrange(match(gene_id, rownames(mat_key_scaled)))
# rownames(mat_key_scaled) <- gene_id_symbol$gene_symbol


de_up_gene_ids <- DE_up_genes$gene_id %>% unique()
de_down_gene_ids <- DE_down_genes$gene_id %>% unique()


gene_direction <- ifelse(rownames(mat_key_scaled)%in% de_up_gene_ids, 
                         "Upregulated", "Downregulated")
names(gene_direction) <- rownames(mat_key_scaled)

gene_direction <- gene_direction[rownames(mat_key_scaled)]
row_order <- order(
  factor(gene_direction, levels = c("Upregulated", "Downregulated"))
)

pdf(file = paste0(dir, "/DE_heatmap_CH_changeCol_purityAdj_v4_r2.pdf"),
    width = 20/2.54, height = 10/2.54)
Heatmap(
  mat_key_scaled,
  # mat_key,
  name = "Z-Score", 
        cluster_rows = TRUE, 
  row_split = factor(
    gene_direction,
    levels = c("Upregulated", "Downregulated")
  ),column_km = 2,
        cluster_columns = TRUE, 
        # clustering_distance_columns = "manhattan",  #canberra manhattan
        # clustering_method_columns = "complete", #ward.D2
        # clustering_method_rows = "complete",
        clustering_method_columns = "ward.D2",
        top_annotation = column_ha, 
        right_annotation = row_ha, 
        # show_row_names = TRUE, 
        show_row_names = FALSE,
        show_column_names = FALSE, 
        # column_split = 6,
        # column_km = 6,
        show_row_dend = TRUE, 
        row_dend_width = unit(10, "mm"),
        show_column_dend = FALSE, 
        # heatmap_legend_param = list(
        #   legend_direction = "horizontal", 
        #   legend_width = unit(6, "cm")), 
        col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
        row_title = "Genes", column_title = "Samples")

dev.off()

# draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right")

















# =========================
# Inputs / sanity checks
# =========================
stopifnot(is.matrix(mat_key_scaled))
stopifnot(exists("MM_metadata_update"))
stopifnot("L1PA2" %in% colnames(MM_metadata_update))
stopifnot(exists("de_gene_keyword_final"))
stopifnot(all(c("de_gene_id", "Pathway") %in% colnames(de_gene_keyword_final)))

samples <- colnames(mat_key_scaled)
genes   <- rownames(mat_key_scaled)

# =========================
# Helper 1: find sample-id column in metadata automatically
# =========================
find_sample_id_col <- function(meta, samples) {
  cand_cols <- names(meta)[sapply(meta, function(z) is.character(z) || is.factor(z))]
  scores <- list()
  for (cn in cand_cols) {
    vals <- as.character(meta[[cn]])
    scores[[cn]] <- sum(samples %in% vals)
  }
  rn <- rownames(meta)
  if (!is.null(rn)) scores[[".rownames"]] <- sum(samples %in% rn)
  
  sc <- unlist(scores)
  best <- names(sc)[which.max(sc)]
  best_n <- max(sc)
  list(best = best, best_n = best_n, scores = sort(sc, decreasing = TRUE))
}

# =========================
# Helper 2: coerce numeric robustly (handles list-column)
# =========================
as_numeric_vec <- function(x) {
  x <- unlist(x, use.names = FALSE)
  suppressWarnings(as.numeric(as.character(x)))
}

# =========================
# Align metadata -> heatmap columns
# =========================
res_id <- find_sample_id_col(MM_metadata_update, samples)

message("==== DEBUG: sample-id matching ====")
message("Best key: ", res_id$best, "  (matched ", res_id$best_n, " / ", length(samples), ")")
print(head(res_id$scores, 10))

if (res_id$best_n == 0) {
  stop("No sample identifier in MM_metadata_update matches colnames(mat_key_scaled). ",
       "Please check heatmap column names vs metadata sample IDs.")
}

meta2 <- MM_metadata_update
if (res_id$best == ".rownames") {
  meta2$.sample_id <- rownames(meta2)
} else {
  meta2$.sample_id <- as.character(meta2[[res_id$best]])
}

meta2 <- meta2[!is.na(meta2$.sample_id) & meta2$.sample_id != "", , drop = FALSE]
meta2 <- meta2[!duplicated(meta2$.sample_id), , drop = FALSE]
rownames(meta2) <- meta2$.sample_id

# Extract L1PA2 aligned to heatmap samples
l1_vec <- as_numeric_vec(meta2[samples, "L1PA2", drop = TRUE])
names(l1_vec) <- samples

message("==== DEBUG: L1PA2 alignment ====")
message("Non-NA L1PA2: ", sum(!is.na(l1_vec)), " / ", length(l1_vec))
if (anyNA(l1_vec)) {
  bad <- names(l1_vec)[is.na(l1_vec)]
  message("Example NA samples (up to 20):")
  print(head(bad, 20))
  stop("L1PA2 has NA after alignment. Fix sample IDs so metadata matches heatmap columns.")
}

# =========================
# Row annotation: map Pathway to heatmap rows
# =========================
path_map <- setNames(as.character(de_gene_keyword_final$Pathway),
                     as.character(de_gene_keyword_final$de_gene_id))

pathway_vec <- path_map[genes]
matched_n <- sum(!is.na(pathway_vec))
message("==== DEBUG: Pathway mapping ====")
message("Matched rows: ", matched_n, " / ", length(genes),
        " (", round(matched_n/length(genes)*100, 1), "%)")

if (matched_n / length(genes) < 0.2) {
  message("WARNING: Pathway matching rate is low (<20%). ",
          "rownames(mat_key_scaled) likely do NOT match de_gene_keyword_final$de_gene_id.",
          "\nFix: make them the same ID type (both gene symbol or both Ensembl).")
}

pathway_vec[is.na(pathway_vec)] <- "Unclassified"

row_ha <- rowAnnotation(
  Pathway = pathway_vec,
  col = list(Pathway = c(
    "E2F" = "pink",
    "INTERFERON" = "lightblue",
    "P53" = "darkgreen",
    "Unclassified" = "white"
  )),
  annotation_name_side = "top",
  annotation_name_gp = gpar(col = NA),
  width = unit(0.5, "cm")
)

# =========================
# Column annotation: L1PA2 gradient bars (blue -> white -> brown)
# =========================
q <- quantile(l1_vec, probs = c(0.05, 0.50, 0.95), na.rm = TRUE)
col_fun_L1 <- colorRamp2(c(q[1], q[2], q[3]), c("blue", "white", "brown"))
bar_cols <- col_fun_L1(l1_vec)

column_ha <- HeatmapAnnotation(
  L1PA2 = anno_barplot(
    l1_vec,
    gp = gpar(fill = bar_cols, col = NA),
    border = FALSE
  ),
  height = unit(1.5, "cm"),
  na_col = "white"
)

# =========================
# Column clustering + GUARANTEED flip (high L1PA2 on RIGHT)
# =========================
COL_DIST_METHOD   <- "manhattan"
COL_HCLUST_METHOD <- "complete"

col_dist   <- dist(t(mat_key_scaled), method = COL_DIST_METHOD)
col_hclust <- hclust(col_dist, method = COL_HCLUST_METHOD)
col_dend   <- as.dendrogram(col_hclust)

ord <- order.dendrogram(col_dend)
ord_samples <- samples[ord]

k <- max(5, floor(length(ord_samples) * 0.25))
left_set  <- ord_samples[1:k]
right_set <- ord_samples[(length(ord_samples)-k+1):length(ord_samples)]

left_mean  <- mean(l1_vec[left_set],  na.rm = TRUE)
right_mean <- mean(l1_vec[right_set], na.rm = TRUE)

flipped <- FALSE
if (right_mean < left_mean) {
  col_dend <- rev(col_dend)
  flipped <- TRUE
}

message("==== DEBUG: flip check by ORDER ====")
message("Left 25% mean L1PA2  = ", signif(left_mean, 4))
message("Right 25% mean L1PA2 = ", signif(right_mean, 4))
message("Flipped dendrogram? ", flipped)

# =========================
# Heatmap
# =========================

ht <- Heatmap(
  mat_key_scaled, name = "Z-Score",
  cluster_rows = TRUE,
  cluster_columns = col_dend,
  top_annotation = column_ha,
  right_annotation = row_ha,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = TRUE,
  row_dend_width = unit(10, "mm"),
  show_column_dend = FALSE,
  col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
  row_title = "Genes", column_title = "Samples"
)

pdf(file = paste0(dir, "/DE_heatmap_CH_changeCol_2.pdf"),
    width = 16/2.54, height = 16/2.54)
draw(ht)
dev.off()

# Optional save
# pdf(file.path(output_dir, "heatmap_L1PA2_high_on_right.pdf"), width = 12, height = 10)
# draw(ht)
# dev.off()










