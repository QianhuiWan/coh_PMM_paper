# use count matrix from feature counts
# get gene counts stats, MMRF
# get DE topTable
pathway_genes_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2"
# dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR1"
dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_redo_res"
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
  read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_redo_res/all_activeL1PA2_patientsIn4Groups_764.tsv")

MM_metadata_update <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/clinical_data_BMprimary_addID_addSubtypes.tsv")

# patientsIn4Groups <- patientsIn4Groups %>% 
#   dplyr::filter(sampleID %in% MM_metadata$sampleID) %>% 
#   arrange(match(sampleID, MM_metadata$sampleID))

# load topTable output
topTable_dat <- read_tsv(paste0(dir, "/lmVoom_topTable_Q4Q1_dat.tsv")) %>% 
  dplyr::select(gene_id, logFC, FDR)

lcpm <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s1_checkStrandness_R4/lcpm_TMM_764_STARcounts_cpmFilterApplied.tsv")


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

pdf(file = paste0(dir, "/DE_heatmap_CH_changeCol.pdf"),
    width = 16/2.54, height = 16/2.54)
Heatmap(mat_key_scaled, name = "Z-Score", 
        cluster_rows = TRUE, 
        cluster_columns = TRUE, 
        clustering_distance_columns = "manhattan",  #canberra
        # clustering_method_columns = "complete",
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

     




