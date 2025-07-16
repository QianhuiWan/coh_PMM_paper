
# heatmap
# use s6_activationScore_R3 output
# only plot L1PA2

input_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/activationScore_R2/filtered_L1"
output_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/activationScore_R2/filtered_L1"

# .libPaths("/home/qwan/miniconda3/envs/coh/lib/R/library")
library(tidyverse)
library(data.table)
library(dplyr)
library(readr)
library(rlang)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(patchwork)

options(scipen=999)

# read in tsv files
## get path of tsv files
activeL1_path <- list.files(input_dir, recursive = TRUE, 
                            pattern = paste0(".*\\.tsv$"), 
                            full.names = TRUE)

activeL1_list <- list()
for (i in 1:length(activeL1_path)) {
  sample <-  str_remove_all(basename(activeL1_path[i]), "_.+")
  # tsv_name_i <- str_replace_all(basename(activeL1_path[i]), "^[^_]+_", "")
  tsv_name_i <- str_remove_all(basename(activeL1_path[i]), ".tsv")
  tsv_file <- read_tsv(file = activeL1_path[i])
  tsv_file$sampleID <- rep(sample, nrow(tsv_file))
  activeL1_list[[tsv_name_i]] <- tsv_file
}

# Function to match classes to the first data frame
match_classes <- function(df, ref_df) {
  for (col in names(ref_df)) {
    if (col %in% names(df)) {
      df[[col]] <- as(df[[col]], class(ref_df[[col]]))
    }
  }
  return(df)
}

# Apply the function to the list
ref_df <- activeL1_list[[1]] # Reference data frame (the first in the list)
activeL1_list <- lapply(activeL1_list, match_classes, ref_df = ref_df)

activeL1_combined_df <- activeL1_list %>% 
  dplyr::bind_rows(.id = "group") 

# get default situation (named sense) L1 activation
activeL1_combined_sense <- activeL1_combined_df %>% 
  dplyr::filter(str_detect(group, "_sense_L1score"))

write_tsv(activeL1_combined_sense, 
          file = paste0(output_dir, "/sense_L1s.tsv"))

# get antisense L1 activation
activeL1_combined_antisense <- activeL1_combined_df %>% 
  dplyr::filter(str_detect(group, "_antisense_L1score"))

write_tsv(activeL1_combined_antisense, 
          file = paste0(output_dir, "/antisense_L1s.tsv"))



activeL1_combined_sense <- read_tsv(paste0(output_dir, "/sense_L1s.tsv"))
activeL1_combined_antisense <- read_tsv(paste0(output_dir, "/antisense_L1s.tsv"))

# TE regions
## sense (based on stranded library)
te_sense_L1s_mat <- activeL1_combined_sense %>% 
  # dplyr::select(Geneid, sampleID, Ratio_Sense)
  pivot_wider(id_cols = Geneid, names_from = sampleID, 
              values_from = Ratio_Sense) %>% 
  tibble::column_to_rownames(var="Geneid") %>% 
  as.matrix() 
te_sense_L1s_mat[is.na(te_sense_L1s_mat)] <- 0


breaks <-c(seq(min(te_sense_L1s_mat), 1, length.out =10), seq(1,5, length.out =10), 
           seq(5, max(te_sense_L1s_mat), length.out =10)) %>% unique()

custom_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(30)

row_dist <- dist(te_sense_L1s_mat)  # Calculate row distances
row_clust <- hclust(row_dist, method = "complete") 
# Flip the row order (reverse the order of the dendrogram)
row_clust$order <- rev(row_clust$order)

p1 <- pheatmap::pheatmap(te_sense_L1s_mat, #scale = "row",
                        # clustering_method = "ward.D2",
                        cluster_rows  = row_clust,
                        clustering_distance_cols = "euclidean",
                        color = custom_colors, breaks = breaks,
                        treeheight_row = 0, treeheight_col = 15,
                        border_color = NA, fontsize = 8,
                        # fontfamily = "Arial",
                        show_rownames = FALSE)

ggsave(paste0(output_dir, "/sense_L1s_heatmap.pdf"), 
       plot = p1$gtable, width = 5, height = 16/5.29, units = "cm", limitsize = FALSE)


## antisense (based on stranded library)

te_antisense_L1s_mat <- activeL1_combined_antisense %>% 
  pivot_wider(id_cols = Geneid, names_from = sampleID, 
              values_from = Ratio_Antisense) %>% 
  tibble::column_to_rownames(var="Geneid") %>% 
  as.matrix() 
te_antisense_L1s_mat[is.na(te_antisense_L1s_mat)] <- 0

breaks <-c(seq(min(te_antisense_L1s_mat), 1, length.out =10), seq(1,3, length.out =10), seq(2, max(te_antisense_L1s_mat), length.out =10)) %>% unique()

custom_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(30)

row_dist <- dist(te_antisense_L1s_mat)  # Calculate row distances
row_clust <- hclust(row_dist, method = "ward.D2") 
# Flip the row order (reverse the order of the dendrogram)
row_clust$order <- rev(row_clust$order)

col_dist <- dist(t(te_antisense_L1s_mat), method = "canberra")  # Calculate row distances
col_clust <- hclust(col_dist, method = "ward.D2") 
# Flip the row order (reverse the order of the dendrogram)
col_clust$order <- rev(col_clust$order)

p2 <- pheatmap::pheatmap(te_antisense_L1s_mat, #scale = "row", 
                        clustering_method = "ward.D2",
                        # cluster_rows  = row_clust, 
                        cluster_cols = col_clust,
                        clustering_distance_cols = "canberra",
                        treeheight_row = 0, treeheight_col = 15,
                        border_color = NA, fontsize = 9.6,
                        show_rownames = FALSE, color = custom_colors, 
                        breaks = breaks)

ggsave(paste0(output_dir, "/antisense_L1s_heatmap.pdf"), 
       plot = p2$gtable, width = 5, height = 7.7, units = "cm", limitsize = FALSE)




# plot TE regions that have both sense and antisense activation
te_sense_L1s_mat_fil <- 
  te_sense_L1s_mat[rownames(te_sense_L1s_mat) %in% 
                     rownames(te_antisense_L1s_mat), ]

te_antisense_L1s_mat_fil <- 
  te_antisense_L1s_mat[rownames(te_antisense_L1s_mat) %in% 
                         rownames(te_sense_L1s_mat), ]

breaks <- c(seq(min(te_sense_L1s_mat_fil), 1, length.out =10), seq(1,2, length.out =10), seq(2, max(te_sense_L1s_mat_fil), length.out =10)) %>% unique()

row_dist <- dist(te_sense_L1s_mat_fil)  # Calculate row distances
col_dist <- dist(t(te_sense_L1s_mat_fil), method = "euclidean")  # Calculate row distances
col_clust <- hclust(col_dist, method = "ward.D2") 
# Flip the row order (reverse the order of the dendrogram)
col_clust$order <- c(col_clust$order[2:length(col_clust$order)], col_clust$order[1])

p3 <- pheatmap::pheatmap(te_sense_L1s_mat_fil, #scale = "row", 
                         clustering_method = "ward.D2",
                         cluster_cols = col_clust,
                         # clustering_distance_cols = "manhattan",
                         show_rownames = FALSE, cluster_rows = FALSE,
                         treeheight_col = 20, 
                         color = custom_colors, breaks = breaks)
ggsave(paste0(output_dir, "/bothSenseAntisense_L1s_sense_heatmap.pdf"), 
       plot = p3$gtable, width = 5, height = 7, units = "cm", limitsize = FALSE)



breaks <- c(seq(min(te_antisense_L1s_mat_fil), 1, length.out =10), seq(1,2, length.out =10), 
           seq(2, max(te_antisense_L1s_mat_fil), length.out =10)) %>% unique()

# row_dist <- dist(te_antisense_L1s_mat_fil)  # Calculate row distances
# col_dist <- dist(t(te_antisense_L1s_mat_fil), method = "euclidean")  # Calculate row distances
# col_clust <- hclust(col_dist, method = "ward.D2") 
# # Flip the row order (reverse the order of the dendrogram)
# col_clust$order <- c(col_clust$order[2:length(col_clust$order)], col_clust$order[1])

p4 <- pheatmap::pheatmap(te_antisense_L1s_mat_fil, #scale = "row", 
                         clustering_method = "ward.D2",
                         cluster_cols = col_clust,
                         # clustering_distance_cols = "maximum",
                         # clustering_distance_cols = "manhattan",
                         show_rownames = FALSE, cluster_rows = FALSE,
                         treeheight_col = 20, 
                         color = custom_colors, breaks = breaks)
ggsave(paste0(output_dir, "/bothSenseAntisense_L1s_antisense_heatmap.pdf"), 
       plot = p4$gtable, width = 5, height = 7, units = "cm", limitsize = FALSE)


# save te_sense_L1s_mat
write_tsv(rownames_to_column(as.data.frame(te_sense_L1s_mat), var = "te_id"), 
          file = paste0(output_dir, "/te_sense_L1s_mat_fil.tsv"))

# save te_antisense_L1s_mat
write_tsv(rownames_to_column(as.data.frame(te_antisense_L1s_mat), var = "te_id"), 
          file = paste0(output_dir, "/te_antisense_L1s_mat_fil.tsv"))









