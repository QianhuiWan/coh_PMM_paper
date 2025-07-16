# previous round did not have all KZFPs
# retrieve all KZFPs from biomart in this round4
# Note: output from same pipeline as MMRF analysis

# load packages
library(tidyverse)
library(magrittr)
library(stringr)
library(tidyr)
library(biomaRt)
library(RColorBrewer)
# laod packages
library(rtracklayer)
library(ggplot2)
library(pheatmap)
library(data.table)
library(plyranges)
options(scipen=999)


# cor plot between ZNF and TEtx ################################################
# load files
# tetx_df <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_star_allo_round4/counts_tx_round3/tetx.tsv")

data <- import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf")
data <- as.data.frame(data, stringsAsFactors = FALSE)
data <- data[data$type == "transcript", ]

# L1binding_KZFPs <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/all_overlaps_L1_binding_KZFPs_PMM.tsv")
# L1binding_KZFPs <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/all_overlaps_L1_binding_KZFPs_PMM_R3.tsv")
L1binding_KZFPs <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/all_overlaps_L1_binding_KZFPs_PMM_R4.tsv")

L1binding_KZFPs <- L1binding_KZFPs %>%
  group_by(KZFP_name) %>%       # Group by KZFP_name
  filter(n() >= 10) %>%
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
L1_tetx_num <- read_tsv(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/activationScore_R2/filtered_L1/num_active_all_L1s_L1PA2_L1PA3_pmm.tsv")

L1_TEtx_num <- L1_tetx_num %>% 
  dplyr::select(sample = sampleID, L1, L1PA2, L1PA3) %>% unique()


# get annotation file: repName and repFam match
# tetx_teAnno <- is only expression repName and repeatFam


# KZFPs
# lcpm_df <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_star_allo_round4/counts_tx_round4/lcpm_TMM.tsv")
lcpm_PMM <- readr::read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/lcpm_TMM_PMM_STARcounts.tsv")
deTopTable_dat <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/topTable_res.tsv")

## annotate ZNF with chr, annotate repeat with repeat age
### retrive KZFP
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") #i.e. Human genes (GRCh38.p14)
### get C2H2 ZFP and KZFP genes anno
### Retrieve genes with the C2H2-type zinc finger domain (InterPro ID: IPR007087)
c2h2_zfps <- getBM(
  attributes = c("hgnc_symbol", 
                 "ensembl_gene_id", "ensembl_gene_id_version",
                 "gene_biotype", 
                 "interpro", 
                 "description",
                 "external_gene_name"
  ),
  filters = "interpro",
  values = "IPR013087", # this is the C2H2 domain
  mart = ensembl
)

### Retrieve all genes with KRAB domain, then filter for KZFP
krab_genes <- getBM(
  attributes = c("hgnc_symbol", 
                 "ensembl_gene_id", "ensembl_gene_id_version",
                 "gene_biotype", 
                 "interpro", "description", "external_gene_name"),
  filters = "interpro",
  values = "IPR001909", # this is the KRAB domain
  mart = ensembl
)
### filter KZFPs from c2h2_zfps
kzfps <- c2h2_zfps[c2h2_zfps$ensembl_gene_id%in%krab_genes$ensembl_gene_id,] #410
kzfps <- kzfps[!(kzfps$hgnc_symbol==""),]


ZNFs_lcpm_L1_tetx <- lcpm_PMM %>% 
  dplyr::mutate(gene_name_original = gene_name) %>% 
  dplyr::mutate(gene_name = str_replace_all(gene_name_original, "_.+", "")) %>% 
  dplyr::filter(gene_name %in% kzfps$hgnc_symbol) %>% 
  dplyr::filter(gene_name %in% unique(L1binding_KZFPs$KZFP_name)) %>%
  # dplyr::filter(gene_name %in% unique(L1PA2binding_KZFPs$KZFP_name)) %>% 
  mutate(lcpm = ifelse(lcpm < 0, 0, lcpm)) %>% 
  pivot_wider(names_from = gene_name, values_from = lcpm,
              id_cols = c(sample)) %>% 
  left_join(L1_TEtx_num[,c("sample", "L1PA2")], "sample") %>%  
  # dplyr::filter(str_detect(sampleID, "PMM")) %>% 
  # replace_na(list(TEtranscript_num = 0)) %>% 
  column_to_rownames(var = "sample")


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
  dplyr::select(ZNF141:ZNF182) %>% 
  as.matrix() %>% 
  t()


# function to get ZNF gene annotations
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
    dplyr::left_join(deTopTable_dat, by = c("gene_name"="gene_symbol")) %>% 
    dplyr::mutate(moreInBMPC = ifelse(logFC < 0, "Yes", "No")) %>% 
    dplyr::mutate(`logFCinBMPCabove1` = ifelse(logFC < -1, "Yes", "No"))
  
  # add PMM lcpm SD and SD/mean 
  CV_cutoff <- lcpm_PMM %>% 
    dplyr::mutate(gene_name_original = gene_name) %>% 
    dplyr::mutate(gene_name = str_replace_all(gene_name_original, "_.+", "")) %>% 
    mutate(lcpm = ifelse(lcpm < 0, 0, lcpm)) %>% 
    dplyr::filter(str_detect(sample, "PMM")) %>% 
    dplyr::group_by(gene_name) %>% 
    dplyr::summarise(SD_lcpm = sd(lcpm, na.rm =TRUE),
                     mean_lcpm = mean(lcpm, na.rm =TRUE),
                     SD_mean_lcpm = SD_lcpm/mean_lcpm)
  CV_cutoff <- median(CV_cutoff$SD_mean_lcpm, na.rm = TRUE)
  
  pmm_lcpm_SD_mean <- lcpm_PMM %>% 
    dplyr::mutate(gene_name_original = gene_name) %>% 
    dplyr::mutate(gene_name = str_replace_all(gene_name_original, "_.+", "")) %>%
    dplyr::filter(gene_name %in% ZNF_geneAnno$gene_name) %>% 
    dplyr::arrange(match(gene_name, ZNF_geneAnno$gene_name)) %>% 
    mutate(PMM_lcpm = ifelse(lcpm < 0, 0, lcpm)) %>% 
    dplyr::filter(str_detect(sample, "PMM")) %>% 
    dplyr::group_by(gene_name) %>% 
    dplyr::summarise(PMM_SD_lcpm = sd(PMM_lcpm, na.rm =TRUE),
                     PMM_mean_lcpm = mean(PMM_lcpm, na.rm =TRUE),
                     PMM_SD_mean_lcpm = PMM_SD_lcpm/PMM_mean_lcpm)
  
  BMPC_lcpm_SD_mean <- lcpm_PMM %>% 
    dplyr::mutate(gene_name_original = gene_name) %>% 
    dplyr::mutate(gene_name = str_replace_all(gene_name_original, "_.+", "")) %>%
    dplyr::filter(gene_name %in% ZNF_geneAnno$gene_name) %>% 
    dplyr::arrange(match(gene_name, ZNF_geneAnno$gene_name)) %>% 
    mutate(BMPC_lcpm = ifelse(lcpm < 0, 0, lcpm)) %>% 
    dplyr::filter(str_detect(sample, "BMPC")) %>% 
    dplyr::group_by(gene_name) %>% 
    dplyr::summarise(BMPC_SD_lcpm = sd(BMPC_lcpm, na.rm =TRUE),
                     BMPC_mean_lcpm = mean(BMPC_lcpm, na.rm =TRUE),
                     BMPC_SD_mean_lcpm = BMPC_SD_lcpm/BMPC_mean_lcpm)
  
  cluster_info <- kzfp_genes_GR %>% 
    join_overlap_left(KZFP_cluster_chr19_GR) %>% 
    plyranges::filter(seqnames %in% ZNF_geneAnno$seqnames) %>% 
    plyranges::filter(gene_name %in% ZNF_geneAnno$gene_name) %>% 
    plyranges::arrange(match(gene_name, ZNF_geneAnno$gene_name)) %>% 
    as.data.frame() %>% 
    dplyr::select(gene_name, clusterID)
  
  ZNF_geneAnno <- ZNF_geneAnno %>% 
    dplyr::left_join(pmm_lcpm_SD_mean, by = "gene_name") %>% 
    dplyr::left_join(BMPC_lcpm_SD_mean, by = "gene_name") %>% 
    dplyr::left_join(cluster_info, by = "gene_name") %>% 
    dplyr::mutate(`VarInPMMabove0.2` = ifelse(PMM_SD_mean_lcpm > 0.2, "Yes", "No"))
    # dplyr::mutate(highVarInPMM = ifelse(SD_mean_lcpm > CV_cutoff, "Yes", "No"))
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

# keep_ZNF <- rownames(matrix_subset_ZNFsX) %in% ZNF_geneAnno$gene_name
# matrix_subset_ZNFsX <- matrix_subset_ZNFsX[keep_ZNF, keep_ZNF]
# table(rownames(matrix_subset_ZNFsX) == ZNF_geneAnno$gene_name)
table(colnames(matrix_subset_ZNFsX) == ZNF_geneAnno$gene_name)

#### Row annotation: Chromosome information
annotation_row <- data.frame(
  chromosome = ZNF_geneAnno$seqnames,
  # `logFCinBMPCabove1` = ZNF_geneAnno$logFCinBMPCabove1,
  logFC = ZNF_geneAnno$logFC,
  `var_above0.2_PMM` = ZNF_geneAnno$VarInPMMabove0.2,
  mean_lcpm_PMM = ZNF_geneAnno$PMM_mean_lcpm,
  # Chr19_Cluster = ZNF_geneAnno$clusterID,
  row.names = rownames(matrix_subset_ZNFsX)
)

#### color annotation
logFC_colors <-  colorRampPalette(c("darkgreen", "white", "orange2"))(100)
# SD_mean_lcpm_colors <- colorRampPalette(c("white", "blue"))(100)
mean_lcpm_PMM_colors <- colorRampPalette(c("white", "orange", "orange2"))(100)

annotation_colors = list(
  # chromosome = setNames(colorRampPalette(brewer.pal(8, "Set3"))(25), 
  #                       levels(annotation_row$Chromosome)),
  chromosome = setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(unique(annotation_row$chromosome))), 
                        unique(annotation_row$chromosome)),
  logFC = logFC_colors,
  # `logFCinBMPCabove1` = c("Yes" = "darkgreen", "No" = "white"),
  # SD_mean_lcpm = SD_mean_lcpm_colors,
  `var_above0.2_PMM` = c("Yes" = "blue", "No" = "grey"),
  mean_lcpm_PMM = mean_lcpm_PMM_colors
  # Chr19_Cluster = setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(unique(annotation_row$Chr19_Cluster))), 
                           # unique(annotation_row$Chr19_Cluster))
)

# Define symmetric breaks to center 0 to white
max_val <- max(abs(matrix_subset_ZNFsX))  # Find the maximum absolute value
breaks <- seq(-max_val, max_val, length.out = 101)  # Ensure 0 is at the center
# Define custom color palette with white at 0
custom_colors <- colorRampPalette(c("blue", "white", "red"))(100)

p1 <- pheatmap::pheatmap(matrix_subset_ZNFsX, 
                         color = custom_colors, 
                         breaks = breaks, cutree_rows = 3,cutree_cols = 3,
                         # annotation_names_row = FALSE,
                         annotation_row = annotation_row, 
                         annotation_colors = annotation_colors, 
                         treeheight_row = 0)

ggsave("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/L1PA2binding_KZFP_Cor_plot.pdf", 
       plot = p1,
       width = 8, height = 6.5, units = "in", limitsize = FALSE)


p2 <- ZNFs_lcpm_L1_tetx_clean %>% 
  dplyr::select(ZNF141, ZNF382, L1PA2) %>% 
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
  ggrepel::geom_text_repel(aes(label = sample), max.overlaps = 100, 
                           box.padding = 0.5, size = 3)+
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

# ZNF_geneAnno <- ZNF_geneAnno[ZNF_geneAnno$logFC<(-0.5),]
# ZNF_geneAnno <- ZNF_geneAnno[ZNF_geneAnno$PMM_SD_mean_lcpm>0.1,]

# keep_ZNF <- rownames(matrix_subset_ZNFsX) %in% ZNF_geneAnno$gene_name
# matrix_subset_ZNFsX <- matrix_subset_ZNFsX[keep_ZNF,]
# table(rownames(matrix_subset_ZNFsX) == ZNF_geneAnno$gene_name)


# cluster3 KZFPs
cor_c3_KZFPs <- p1$tree_row$labels[p1$tree_row$order][c(15, 19:22)]
keep_ZNF <- rownames(matrix_subset_ZNFsX) %in% cor_c3_KZFPs
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
  row.names = rownames(matrix_subset_ZNFsX)
)

annotation_col <- data.frame(
  L1PA2_num = ZNFs_lcpm_L1_tetx_clean$L1PA2,
  row.names = colnames(matrix_subset_ZNFsX)
)


#### color annotation
annotation_colors = list(
  Chromosome = setNames(c("#66C2A5", "#8DA0CB", "#FFD92F"), 
                        unique(annotation_row$Chromosome)),
  KZFP_age = colorRampPalette(c("yellow", "brown"))(100),
  L1PA2_num = colorRampPalette(c("grey", "darkgreen"))(100)
)

# Define the breaks, ensuring 0 is exactly in the middle
breaks <- c(seq(-3, 0, length.out = 11), seq(0.3, 3, length.out = 10))  # 0 included only once

# Create a custom palette that ensures white is in the center
custom_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(length(breaks) - 1)

row_dist <- dist(matrix_subset_ZNFsX)  # Calculate row distances
row_clust <- hclust(row_dist, method = "ward.D2") 

col_dist <- dist(t(matrix_subset_ZNFsX), method = "canberra")  # Calculate row distances
col_clust <- hclust(col_dist, method = "complete") 
# Flip the row order (reverse the order of the dendrogram)
col_clust$order <- rev(col_clust$order)


p3 <- pheatmap::pheatmap(matrix_subset_ZNFsX, scale = "row", 
                        # clustering_method = "ward.D2",
                        # clustering_distance_cols = "canberra",
                        # clustering_distance_cols = "manhattan",
                        annotation_row = annotation_row, 
                        annotation_col = annotation_col,
                        annotation_colors = annotation_colors, fontsize = 12,
                        breaks = breaks, color = custom_colors,
                        # cluster_cols = col_clust,
                        treeheight_row = 0)

ggsave("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/L1_KZFP_lcpm_plot_addL1PA2num.pdf", 
       plot = p3,
       width = 6.5, height = 5, units = "in", limitsize = FALSE)




# check pathway genes
pathwayGenes_lcpm_L1_tetx <- lcpm_PMM %>% 
  dplyr::mutate(gene_name_original = gene_name) %>% 
  dplyr::mutate(gene_name = str_replace_all(gene_name_original, "_.+", "")) %>% 
  dplyr::filter(gene_name %in% c("NFKB1", "CDKN1A")) %>% 
  mutate(lcpm = ifelse(lcpm < 0, 0, lcpm)) %>% 
  pivot_wider(names_from = gene_name, values_from = lcpm,
              id_cols = c(sample)) %>% 
  left_join(L1_TEtx_num, "sample") %>%  
  # dplyr::filter(str_detect(sampleID, "PMM")) %>% 
  # replace_na(list(TEtranscript_num = 0)) %>% 
  column_to_rownames(var = "sample")

pathwayGeneCor <- pathwayGenes_lcpm_L1_tetx %>% 
  dplyr::select(NFKB1, CDKN1A, L1PA2) %>% 
  tibble::rownames_to_column(var="sample") %>% 
  pivot_longer(
    cols = -c(sample, L1PA2), # Exclude "Sample" and "L1PA2" from being melted
    names_to = "gene",     # New column for variable names
    values_to = "lcpm"       # New column for values
  ) %>% 
  ggplot(aes(x = lcpm, y = L1PA2)) + 
  geom_point() + geom_smooth(method = "lm")+ 
  facet_wrap(~gene, scales = "free_x")+
  ggpubr::stat_cor(label.y=c(135,135), method = "pearson",size=3) +
  ggrepel::geom_text_repel(aes(label = sample), max.overlaps = 10, box.padding = 0.5, size = 3)+
  theme_set(ggpubr::theme_pubr(base_size = 12))

ggsave("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/NFKB1_CDKN1A_L1PA2_Cor.pdf", 
       plot = pathwayGeneCor,
       width = 15, height = 8, units = "cm", limitsize = FALSE)





