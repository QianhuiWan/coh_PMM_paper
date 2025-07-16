# heatmap of all L1 binding KZFPs
# retrieve all KZFPs from biomart in this round4
# USE TCGA count matrix

# load packages
library(tidyverse)
library(magrittr)
library(stringr)
library(tidyr)
# library(biomaRt)
library(RColorBrewer)
# laod packages
library(rtracklayer)
library(ggplot2)
library(pheatmap)
library(data.table)
library(edgeR)
library(plyranges)
options(scipen=999)


# cor plot between ZNF and TEtx ################################################
# load files
data <- import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf")
data <- as.data.frame(data, stringsAsFactors = FALSE)
data <- data[data$type == "transcript", ]

# L1binding_KZFPs <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/all_overlaps_L1_binding_KZFPs_PMM.tsv")
L1binding_KZFPs <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/all_overlaps_L1_binding_KZFPs_MM_R4.tsv")
L1binding_KZFPs <- L1binding_KZFPs %>%
  group_by(KZFP_name) %>%       # Group by KZFP_name
  filter(n() > 10) %>%
  ungroup()
L1PA2binding_KZFPs <- L1binding_KZFPs %>% dplyr::filter(repName == "L1PA2")


# KZFP age
KZFP_age <- read_tsv(paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/KZFP_evolutional_age.tsv"))
# KZFP cluster
KZFP_cluster_chr19 <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/KZFP_clusters_coord_hg38_chr19_noStrand.tsv")
KZFP_cluster_chr19_GR <- KZFP_cluster_chr19 %>% as_granges()
kzfp_genes_GR <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/kzfps_hg38_all_geneCoordinates.tsv") %>% 
  as_granges()
seqlevels(kzfp_genes_GR) <- paste0("chr", seqlevels(kzfp_genes_GR))

cluster_info <- kzfp_genes_GR %>% 
  join_overlap_left(KZFP_cluster_chr19_GR) 

# calculate repeat transcript number
activeL1_combined_sense <- 
  read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_R1res/pred_sense_L1s_step2R2_S6R3.tsv")
activeL1_combined_antisense <- 
  read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_R1res/pred_antisense_L1s_step2R2_S6R3.tsv")

MM_tetx <- rbind(activeL1_combined_sense, activeL1_combined_antisense) 

activeL1_all_forSheet <- MM_tetx %>% 
  dplyr::filter(repFamily == "L1") %>% 
  dplyr::count(sampleID) %>% 
  right_join(MM_tetx, by = "sampleID") %>%
  dplyr::mutate(L1 = n) %>% 
  dplyr::group_by(repName) %>% 
  dplyr::count(sampleID, .drop = FALSE)
# summarise(count = n(), .groups = "drop")  # Summarise counts

activeL1PA2_forSheet <- activeL1_all_forSheet[activeL1_all_forSheet$repName=="L1PA2",] %>% 
  dplyr::mutate(L1PA2 = n) %>% 
  as.data.frame() %>% 
  dplyr::select(sampleID, L1PA2)

activeL1PA3_forSheet <- activeL1_all_forSheet[activeL1_all_forSheet$repName=="L1PA3",] %>% 
  dplyr::mutate(L1PA3 = n) %>% 
  as.data.frame() %>% 
  dplyr::select(sampleID, L1PA3)

activeL1_all_forSheet_save <- MM_tetx %>% 
  dplyr::filter(repFamily == "L1") %>% 
  dplyr::count(sampleID) %>% 
  right_join(MM_tetx, by = "sampleID") %>%
  dplyr::mutate(L1 = n) %>% 
  as.data.frame() %>% 
  dplyr::select(-n) %>% 
  left_join(activeL1PA2_forSheet, by = "sampleID") %>% 
  left_join(activeL1PA3_forSheet, by = "sampleID")

write_tsv(activeL1_all_forSheet_save, 
          file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR1/num_active_all_L1s_L1PA2_L1PA3_mm.tsv")

L1_TEtx_num <- activeL1_all_forSheet_save %>% 
  dplyr::select(sample = sampleID, L1, L1PA2, L1PA3) %>% unique()


# get TCGA count mat, match barcode to clinial data
raw_counts <- read_tsv(paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/raw_counts_geneCounts_MMRF_BMprimary.tsv")) %>%
  tibble::column_to_rownames(var = "gene_id")
clinical_data_addIDs <- read_tsv(paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess", "/clinical_data_BMprimary_addID.tsv")) #764

raw_counts_df_BMprimary <- raw_counts
## create DGElist
### get group info.
if (unique(clinical_data_addIDs$barcode == colnames(raw_counts_df_BMprimary))) {
  groups <- clinical_data_addIDs$L1_quantile_group
  groups[is.na(groups)] <- "nogroup"
}

dge_rawCounts <- DGEList(counts=raw_counts_df_BMprimary,
                         # genes=raw_counts_df_BMprimary[,1:3],
                         #remove.zeros = TRUE,
                         group = groups)

## normalize counts with TMM method
dge_rawCounts <- calcNormFactors(dge_rawCounts, method = "TMM")
normalized_lcpm <- cpm(dge_rawCounts, normalized.lib.sizes = TRUE, log = TRUE)
lcpm_df <- normalized_lcpm
colnames(lcpm_df) <- clinical_data_addIDs$sampleID

deTopTable_dat <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR1/lmVoom_topTable_Q4Q1_dat.tsv") %>% 
  dplyr::mutate(gene_name = str_replace_all(gene_name, ".+_", ""))

## annotate ZNF with chr, annotate repeat with repeat age

kzfps <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/kzfps_hg38_allHaveGeneSymbol.tsv")

lcpm_df %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_name") %>% 
  pivot_longer(cols = -c(gene_name), 
               names_to = "sample", values_to = "lcpm") %>% 
  mutate(lcpm = ifelse(lcpm < 0, 0, lcpm)) %>% 
  write_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR1/lcpm_TMM_not_fil.tsv")


ZNFs_lcpm_L1_tetx <- lcpm_df %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_name") %>% 
  # dplyr::mutate(gene_name = str_replace_all(gene_name, ".+_", "")) %>% 
  dplyr::mutate(gene_name = str_replace_all(gene_name, "_.+", "")) %>% 
  dplyr::filter(gene_name %in% kzfps$hgnc_symbol) %>% 
  dplyr::filter(gene_name %in% unique(L1binding_KZFPs$KZFP_name)) %>%
  # dplyr::filter(gene_name %in% unique(L1PA2binding_KZFPs$KZFP_name)) %>%
  pivot_longer(cols = -c(gene_name), 
               names_to = "sample", values_to = "lcpm") %>% 
  mutate(lcpm = ifelse(lcpm < 0, 0, lcpm)) %>% 
  # dplyr::mutate(sample = str_replace_all(sample, "_hg38", "")) %>% 
  pivot_wider(names_from = gene_name, values_from = lcpm,
              id_cols = c(sample)) %>% 
  left_join(L1_TEtx_num[,c("sample", "L1PA2")], by = "sample") %>%  
  # dplyr::filter(str_detect(sampleID, "PMM")) %>% 
  # replace_na(list(TEtranscript_num = 0)) %>% 
  column_to_rownames(var = "sample") %>% 
  dplyr::filter(!is.na(L1PA2))

ZNFs_lcpm_L1_tetx %>%  rownames_to_column(var="sampleName") %>%
write_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR1/KZFPs_lcpm_L1_tetx.tsv")

# Identify columns with zero variance (all zeros or constant values)
zero_variance_cols <- sapply(ZNFs_lcpm_L1_tetx, function(x) var(x, na.rm = TRUE) == 0)

# View columns with zero variance
names(zero_variance_cols[zero_variance_cols])
# Remove columns with zero variance
ZNFs_lcpm_L1_tetx_clean <- ZNFs_lcpm_L1_tetx[, !zero_variance_cols]

# Calculate the correlation matrix
# cor_matrix <- cor(ZNFs_lcpm_L1_tetx_clean)

# Basic correlation plot
# library(ggcorrplot)


## get ZNF rows and ZNF columns from cor matrix
matrix_subset_ZNFs <- ZNFs_lcpm_L1_tetx_clean %>% 
  # dplyr::select(ZFP69B:ZNF75D) %>% 
  dplyr::select(-L1PA2) %>% 
  as.matrix() %>% 
  t()


# function to get ZNF gene annotations
# cor_matrix_subset <- matrix_subset_ZNFs
ZNF_geneAnno_fun <- function(cor_matrix_subset){
  ZNF_geneAnno <- data %>% 
    dplyr::select(seqnames, start, end, strand, gene_name) %>% 
    unique() %>% 
    dplyr::filter(gene_name %in% kzfps$hgnc_symbol) %>% 
    plyranges::as_granges() %>% 
    plyranges::group_by(gene_name) %>% 
    plyranges::reduce_ranges() %>% 
    as.data.frame() %>% 
    dplyr::filter(gene_name %in% rownames(cor_matrix_subset)) %>% 
    dplyr::arrange(match(gene_name, rownames(cor_matrix_subset))) %>% 
    dplyr::left_join(deTopTable_dat[, 1:6], by = "gene_name") %>% 
    na.omit()
  
  cluster_info <- kzfp_genes_GR %>% 
    join_overlap_left(KZFP_cluster_chr19_GR) %>% 
    plyranges::filter(seqnames %in% ZNF_geneAnno$seqnames) %>% 
    plyranges::filter(gene_name %in% ZNF_geneAnno$gene_name) %>% 
    plyranges::arrange(match(gene_name, ZNF_geneAnno$gene_name)) %>% 
    as.data.frame() %>% 
    dplyr::select(gene_name, clusterID)
  
  ZNF_geneAnno <- ZNF_geneAnno %>% 
    dplyr::left_join(cluster_info, by = "gene_name")
  return(ZNF_geneAnno)
}



# cor plot for ZNFs ############################################################
# Calculate the correlation matrix
matrix_subset_ZNFsX <- cor(t(matrix_subset_ZNFs))
ZNF_geneAnno <- ZNF_geneAnno_fun(matrix_subset_ZNFsX)
## get BMPC ZNFs
# ZNF_geneAnno <- ZNF_geneAnno[ZNF_geneAnno$BMPC_mean_lcpm>1 & 
#                                ZNF_geneAnno$logFC<(-0.5) & ZNF_geneAnno$PMM_SD_mean_lcpm>0.1,]
# ZNF_geneAnno <- ZNF_geneAnno[ZNF_geneAnno$logFC<(-0.5),]
# ZNF_geneAnno <- ZNF_geneAnno[ZNF_geneAnno$FDR<0.05,]

keep_ZNF <- rownames(matrix_subset_ZNFsX) %in% ZNF_geneAnno$gene_name
matrix_subset_ZNFsX <- matrix_subset_ZNFsX[keep_ZNF, keep_ZNF]
table(rownames(matrix_subset_ZNFsX) == ZNF_geneAnno$gene_name)
table(colnames(matrix_subset_ZNFsX) == ZNF_geneAnno$gene_name)

#### Row annotation: Chromosome information
annotation_row <- data.frame(
  chromosome = ZNF_geneAnno$seqnames,
  # `logFCinBMPCabove1` = ZNF_geneAnno$logFCinBMPCabove1,
  logFC = ZNF_geneAnno$logFC,
  # `var_above0.2_PMM` = ZNF_geneAnno$VarInPMMabove0.2,
  mean_lcpm_MM = ZNF_geneAnno$AveExpr,
  Chr19_Cluster = ZNF_geneAnno$clusterID,
  row.names = rownames(matrix_subset_ZNFsX)
)

#### color annotation
logFC_colors <-  colorRampPalette(c("darkgreen", "white", "orange2"))(100)
# SD_mean_lcpm_colors <- colorRampPalette(c("white", "blue"))(100)
mean_lcpm_PMM_colors <- colorRampPalette(c("white", "orange", "orange2"))(100)

annotation_colors = list(
  chromosome = setNames(c("#FFD92F", "#E5C494", "#66C2A5","#B3B3B3","#F08F6D", 
                          "#8DA0CB", "#D2A29F",  "#E78AC3", "#A6D854"   ), 
                        unique(annotation_row$chromosome)),
  logFC = logFC_colors,
  # `logFCinBMPCabove1` = c("Yes" = "darkgreen", "No" = "white"),
  # SD_mean_lcpm = SD_mean_lcpm_colors,
  # `var_above0.2_PMM` = c("Yes" = "blue", "No" = "grey"),
  mean_lcpm_MM = mean_lcpm_PMM_colors,
  Chr19_Cluster = setNames(c("#66C2A5", "#FFFFFF","#AB98C8", "#C6B18B", "#E1D83B",
                           "#E9C783", "#B3B3B3"),
  unique(annotation_row$Chr19_Cluster))
)

# Define symmetric breaks to center 0 to white
max_val <- max(abs(matrix_subset_ZNFsX))  # Find the maximum absolute value
breaks <- seq(-max_val, max_val, length.out = 101)  # Ensure 0 is at the center
# Define custom color palette with white at 0
custom_colors <- colorRampPalette(c("blue", "white", "red"))(100)

p1 <- pheatmap::pheatmap(matrix_subset_ZNFsX, 
                         color = custom_colors, 
                         breaks = breaks, cutree_rows = 4,cutree_cols = 4,
                         treeheight_col = 20,
                         # annotation_names_row = FALSE,
                         annotation_row = annotation_row,
                         annotation_colors = annotation_colors, 
                         treeheight_row = 0)

ggsave("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR1/MM_KZFP_Cor_plot_R5.pdf", 
       plot = p1,
       width = 9, height = 7, units = "in", limitsize = FALSE)



# cor plot
## only select Q1 and Q4
clinical_data_addIDs_Q1Q4 <- clinical_data_addIDs %>% 
  dplyr::filter(L1_quantile_group %in% c("Q1", "Q4"))
  

p2 <- ZNFs_lcpm_L1_tetx_clean %>% 
  dplyr::filter(rownames(.) %in% clinical_data_addIDs_Q1Q4$sampleID) %>% 
  dplyr::select(ZNF324, ZNF136, 
                # ZNF93, ZNF680,
                ZNF425, ZNF671, ZNF17, ZNF765,
                L1PA2) %>% 
  tibble::rownames_to_column(var="sample") %>% 
  pivot_longer(
    cols = -c(sample, L1PA2), # Exclude "Sample" and "L1PA2" from being melted
    names_to = "gene",     # New column for variable names
    values_to = "lcpm"       # New column for values
  ) %>% 
  ggplot(aes(x = lcpm, y = L1PA2)) + 
  geom_point() + geom_smooth(method = "lm")+ 
  facet_wrap(~gene, scales = "free_x")+
  # ylim(c(-20,150))+
  ggpubr::stat_cor(label.y=c(135,135), method = "pearson",size=3) +
  # ggrepel::geom_text_repel(aes(label = sample), max.overlaps = 100, 
  # box.padding = 0.5, size = 3)+
  ggpubr::theme_pubr(base_size = 9.6)

ggsave("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/ZNF382_ZNF141_L1PA2_Cor.pdf", 
       plot = p2,
       width = 16, height = 8, units = "cm", limitsize = FALSE)




# heatmap for KZFP expression ##################################################
matrix_subset_ZNFsX <- matrix_subset_ZNFs
ZNF_geneAnno <- ZNF_geneAnno_fun(matrix_subset_ZNFsX)

# ZNF_geneAnno <- ZNF_geneAnno[ZNF_geneAnno$BMPC_mean_lcpm>1 & 
#                                ZNF_geneAnno$logFC<(-0.5) & ZNF_geneAnno$PMM_SD_mean_lcpm>0.1,]
# ZNF_geneAnno <- ZNF_geneAnno[ZNF_geneAnno$BMPC_mean_lcpm>1 & 
#                                ZNF_geneAnno$logFC<(-0.5),]

# ZNF_geneAnno <- ZNF_geneAnno[ZNF_geneAnno$logFC<0,]
# ZNF_geneAnno <- ZNF_geneAnno[ZNF_geneAnno$FDR<0.05,]
keep_ZNF <- rownames(matrix_subset_ZNFsX) %in% ZNF_geneAnno$gene_name
matrix_subset_ZNFsX <- matrix_subset_ZNFsX[keep_ZNF,]

cor_c1and3_KZFPs <- p1$tree_row$labels[p1$tree_row$order][c(1:3, 8:11)]
keep_ZNF <- rownames(matrix_subset_ZNFsX) %in% cor_c1and3_KZFPs
ZNF_geneAnno <- ZNF_geneAnno[keep_ZNF,]
matrix_subset_ZNFsX <- matrix_subset_ZNFsX[keep_ZNF,]
table(rownames(matrix_subset_ZNFsX) == ZNF_geneAnno$gene_name)

subset_KZFP_age <- KZFP_age %>% 
  dplyr::filter(KZFP %in% rownames(matrix_subset_ZNFsX)) %>% 
  arrange(match(KZFP, rownames(matrix_subset_ZNFsX)))

#### Row annotation: Chromosome information
annotation_row <- data.frame(
  Chromosome = ZNF_geneAnno$seqnames,
  KZFP_age = subset_KZFP_age$KZFP_age,
  Chr19_Cluster = ZNF_geneAnno$clusterID,
  row.names = rownames(matrix_subset_ZNFsX)
)

annotation_col <- data.frame(
  L1PA2_num = ZNFs_lcpm_L1_tetx_clean$L1PA2,
  row.names = colnames(matrix_subset_ZNFsX)
)


#### color annotation
annotation_colors = list(
  Chromosome = setNames(c("#FFD92F", "#8DA0CB"), 
                        unique(annotation_row$Chromosome)),
  KZFP_age = colorRampPalette(c("yellow", "brown"))(100),
  L1PA2_num = colorRampPalette(c("grey", "darkgreen"))(100),
  Chr19_Cluster = setNames(c("#66C2A5", "#AB98C8",
                             "#E9C783", "#B3B3B3", "#FFFFFF"),
                           unique(annotation_row$Chr19_Cluster))
)

# Define symmetric breaks to center 0 to white
scaled_lcpm <- t(apply(matrix_subset_ZNFsX, 1, scale)) 
colnames(scaled_lcpm) <- colnames(matrix_subset_ZNFsX)
matrix_subset_ZNFsX <- scaled_lcpm
max_val <- max(abs(matrix_subset_ZNFsX))  # Find the maximum absolute value
breaks <- seq(-max_val, max_val, length.out = 101)  # Ensure 0 is at the center
# Define custom color palette with white at 0
custom_colors <- colorRampPalette(c("blue", "white", "red"))(100)

p3 <- pheatmap::pheatmap(matrix_subset_ZNFsX, 
                         # scale = "row",
                         color = custom_colors, 
                         breaks = breaks, 
                         # clustering_method = "ward.D2",
                         # clustering_distance_cols = "canberra",
                         annotation_row = annotation_row, 
                         show_colnames = FALSE,
                         annotation_col = annotation_col, treeheight_col = 20,
                         annotation_colors = annotation_colors, fontsize = 12,
                         treeheight_row = 20
)


ggsave("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR1/L1_KZFP_lcpm_plot_addL1PA2num_MM.pdf", 
       plot = p3,
       width = 9, height = 6, units = "in", limitsize = FALSE)


# test plot
sumval <- colSums(matrix_subset_ZNFs[rownames(matrix_subset_ZNFs)%in% rownames(matrix_subset_ZNFsX), ])
plot(sumval ~ ZNFs_lcpm_L1_tetx_clean$L1PA2)
# Fit a linear regression model
lm_fit <- lm(sumval ~ ZNFs_lcpm_L1_tetx_clean$L1PA2)
summary(lm_fit)
# Add the regression line
abline(lm_fit, col = "red", lwd = 2)



