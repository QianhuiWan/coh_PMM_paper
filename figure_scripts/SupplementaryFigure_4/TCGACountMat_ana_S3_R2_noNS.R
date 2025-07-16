
# Step two (S3), plot KZFP fold change for all datasets

.libPaths()
options(scipen=999)

library(tidyverse)
library(TCGAbiolinks)
library(pheatmap)
library(ComplexHeatmap)
library(cluster)
library(seriation)
library(stats)
library(RColorBrewer)
library(ggplot2)
library(ggplotify)

outRes_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_tcga"

# load tsv
topTable_kzfps_ads <- read_tsv(paste0(outRes_dir, "/topTable_kzfps_allDataSets.tsv"))

# plot KZFP logFC in heatmap for all datasets (as columns)

topTable_kzfps_ads_mat <- 
  topTable_kzfps_ads %>% 
  # dplyr::filter(!Sig=="NS") %>%
  # dplyr::filter(abs(logFC)>=1) %>%
  dplyr::filter(abs(logFC)>=1 & FDR < 0.05) %>%
  pivot_wider(id_cols = gene_name, names_from = dataset, 
              values_from = logFC, values_fill = 0) %>% 
  tibble::column_to_rownames(var="gene_name") %>% 
  as.matrix()


## annotations for columns and plot
### Get TCGA project metadata
projects <- TCGAbiolinks::getGDCprojects()
### Filter for TCGA projects and select relevant columns, get cancer type names
tcga_metadata <- projects[grepl("TCGA", projects$project_id), 
                          c("project_id", "name")] %>% 
  # Select relevant columns
  dplyr::select(Project = project_id, Cancer_type = name) %>%
  # Filter for projects in tcga_projects
  dplyr::filter(Project %in% colnames(topTable_kzfps_ads_mat)) %>%
  # Reorder to match tcga_projects
  arrange(match(colnames(topTable_kzfps_ads_mat), Project)) 

### add original tisue type to cancer type
tcga_metadata$TissueType <- 
  # c("Endocrine", "Digestive", "Digestive", "Soft Tissue", 
  #   "Reproductive", "Digestive", "Urinary", "Head and Neck", "Respiratory", 
  #   "Respiratory", "Urinary", "Urinary", "Urinary", "Digestive", 
  #   "Reproductive", "Digestive", "Reproductive", "Digestive", "Brain", 
  #   "Endocrine", "Skin", "Reproductive")
  # c("Endocrine", "Lymphatic", "Digestive", "Digestive", "Soft Tissue", 
  #   "Reproductive", "Digestive", "Urinary", "Head and Neck", "Respiratory", 
  #   "Liver", "Respiratory", "Urinary", "Urinary", "Urinary", "Digestive", 
  #   "Reproductive", "Digestive", "Reproductive", "Digestive", "Brain", 
  #   "Endocrine", "Skin", "Reproductive")
  c("Endocrine", "Digestive", "Digestive", "Soft Tissue",
    "Reproductive", "Digestive", "Urinary", "Head and Neck", "Respiratory",
    "Liver", "Respiratory", "Urinary", "Urinary", "Urinary", "Digestive",
    "Reproductive", "Digestive", "Reproductive", "Digestive", "Brain",
    "Endocrine", "Skin", "Reproductive")

tcga_metadata$TissueType <- factor(tcga_metadata$TissueType, 
                                   levels = unique(tcga_metadata$TissueType))

# Column annotations (e.g., Sample Type and Condition)
col_annotation <- HeatmapAnnotation(
  TissueType = tcga_metadata$TissueType,
  col = list(
    TissueType = setNames(
      colorRampPalette(brewer.pal(8, "Set3"))(11)[-2], 
      levels(tcga_metadata$TissueType))
  )
)

# Create the heatmap with column annotations
# register_DendSer()
# get_seriation_method("dist", "DendSer")
methods <- list_seriation_methods()
# for (method in methods$dist[c(1,4:8, 11:length(methods$dist))]) {
# for (method in methods$dist[c(11:length(methods$dist))]) {
for (method in methods$dist[c(1)]) {
  o1 = seriate(dist(topTable_kzfps_ads_mat), method = method)
  o2 = seriate(dist(t(topTable_kzfps_ads_mat)), method = method)
  
  p <- Heatmap(
    topTable_kzfps_ads_mat,
    top_annotation = col_annotation,
    name = "LogFC", 
    # clustering_distance_columns = "spearman",
    # clustering_method_columns = "complete",
    # clustering_distance_rows = "maximum",
    # clustering_method_rows= "ward",
    row_order = get_order(o1),
    column_order = get_order(o2),
    # cluster_rows = as.dendrogram(o1[[1]]), 
    # cluster_columns = as.dendrogram(o2[[1]]),
    column_names_centered = TRUE, 
    show_row_names = FALSE
    # cluster_rows = diana(topTable_kzfps_ads_mat),
    # cluster_columns = agnes(t(topTable_kzfps_ads_mat))
  )
  
  # Convert to ggplot object
  ggplot_heatmap <- as.ggplot(p)
  ggsave(paste0(outRes_dir, "/pdf_test/kzfps_FC_TCGAds_", method, ".pdf"), 
         plot = ggplot_heatmap,
         width = 20, height = 15, units = "cm", limitsize = FALSE)
}


library(circlize)
color_scale <- colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
o1 = seriate(dist(topTable_kzfps_ads_mat), method = "OLO_ward")
o2 = seriate(dist(t(topTable_kzfps_ads_mat)), method = "OLO_ward")
# o1 = seriate(dist(topTable_kzfps_ads_mat), method = "R2E")
# o2 = seriate(dist(t(topTable_kzfps_ads_mat)), method = "R2E")
p <- Heatmap(
  topTable_kzfps_ads_mat,
  top_annotation = col_annotation,
  name = "LogFC", 
  # clustering_distance_columns = "spearman",
  # clustering_method_columns = "complete",
  # clustering_distance_rows = "maximum",
  # clustering_method_rows= "complete",
  # row_order = get_order(o1),
  # column_order = get_order(o2),
  cluster_rows = as.dendrogram(o1[[1]]),
  cluster_columns = rev(as.dendrogram(o2[[1]])),
  column_names_centered = TRUE, 
  show_row_names = FALSE,
  # col = color_scale
  # cluster_rows = diana(topTable_kzfps_ads_mat),
  # cluster_columns = agnes(t(topTable_kzfps_ads_mat))
)

# Convert to ggplot object
ggplot_heatmap <- as.ggplot(p)
ggsave(paste0(outRes_dir, "/pdf_test/kzfps_FCabove1_FDR05_TCGAds_OLO_ward.pdf"), 
       plot = ggplot_heatmap,
       width = 16, height = 10, units = "cm", limitsize = FALSE)


# get row clusters
# Convert dendrogram to hclust
hc <- as.hclust(as.dendrogram(o1[[1]]))

# Cut into k clusters (e.g., 2)
clusters <- cutree(hc, k = 2)


topKZFPinCommon <- o1[[1]]$labels[o1[[1]]$order][1:24]


p2 <- Heatmap(
  topTable_kzfps_ads_mat[rownames(topTable_kzfps_ads_mat)%in%topKZFPinCommon,],
  top_annotation = col_annotation,
  name = "LogFC", 
  # clustering_distance_columns = "spearman",
  # clustering_method_columns = "complete",
  # clustering_distance_rows = "maximum",
  # clustering_method_rows= "complete",
  # row_order = get_order(o1),
  # column_order = get_order(o2),
  # cluster_rows = as.dendrogram(o1[[1]]),
  cluster_columns = rev(as.dendrogram(o2[[1]])),
  column_names_centered = TRUE, 
  show_row_names = TRUE,
  col = color_scale
  # cluster_rows = diana(topTable_kzfps_ads_mat),
  # cluster_columns = agnes(t(topTable_kzfps_ads_mat))
)

# Convert to ggplot object
ggplot_heatmap_p2 <- as.ggplot(p2)
ggsave(paste0(outRes_dir, "/pdf_test/kzfps_FCabove1_FDR05_TCGAds_OLO_ward_topCluster.pdf"), 
       plot = ggplot_heatmap_p2,
       width = 22, height = 15, units = "cm", limitsize = FALSE)




