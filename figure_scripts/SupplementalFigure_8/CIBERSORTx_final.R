
# devtools::install_github('andygxzeng/BoneMarrowMap')

library(tidyverse)
library(stringr)
library(Seurat)
library(Matrix)
library(purrr)
library(tibble)
# library(BoneMarrowMap)

# prepare mean counts per BM cell type (single cell)
bm_sc_ref <- readRDS("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/cibersort/references/BoneMarrowMap_Annotated_Dataset_expandedFeatures.rds")

## get mean counts for 24 CellType_Broad types
# Extract raw counts
counts <- GetAssayData(
  bm_sc_ref,
  assay = "RNA",
  layer = "counts")

# Optional: remove genes with zero expression across all cell types
counts <- counts[
  Matrix::rowSums(counts) > 0,]

# Save for CIBERSORTx
cell_type <- bm_sc_ref@meta.data[colnames(counts), "CellType_Broad"]
cell_type <- gsub(" ", "_", cell_type)
cell_type <- gsub("\\+", "pos", cell_type)
cell_type <- gsub("-", "_", cell_type)

names(cell_type) <- colnames(counts)

set.seed(123)
max_cells <- 100
# max_cells <- 50

selected_cells <- unname(unlist(tapply(names(cell_type), cell_type, function(x) {
  sample(x, min(length(x), max_cells))
})))

counts_small <- counts[, selected_cells, drop = FALSE]

# remove genes not expressed in selected cells
# keep_genes <- Matrix::rowSums(counts_small > 0) >= 24 # 0.01 * 2400

lib_size <- Matrix::colSums(counts_small)
counts_cpm <- t(t(counts_small) / lib_size * 1e6)
log_cpm <- log1p(counts_cpm)

gene_var <- Matrix::rowMeans(log_cpm^2) - Matrix::rowMeans(log_cpm)^2
keep_genes <- gene_var >= quantile(gene_var, 0.75)

counts_small <- counts_small[keep_genes, , drop = FALSE]

ref_counts <- data.frame(
  GeneSymbol = rownames(counts_small),
  as.data.frame(as.matrix(counts_small), check.names = FALSE),
  check.names = FALSE
)

out_file <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/cibersort/references/BM_reference_downsampled_100cells_per_typeBroad_raw_counts.txt"
# out_file <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/cibersort/references/BM_reference_downsampled_50cells_per_typeBroad_raw_counts.txt"

header_labels <- c("GeneSymbol", cell_type[selected_cells])

writeLines(
  paste(header_labels, collapse = "\t"),
  con = out_file
)

write.table(
  ref_counts,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  append = TRUE)

# There is a 50mb file size limit for online cibersortx sw

# prepare counts per MM and PMM bulk sample
mm_counts <- read.table("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s1_checkStrandness_R4_gemini_uniqueMap/gene_count_matrix_STARcounts_MM_764.tsv", 
                        header=T,row.names=1,sep="\t", check.names = FALSE)

pmm_counts <- read.table("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/gene_counts/gene_count_matrix_STARcounts_PMM_R3_addBMPC.tsv", 
                         header=T,row.names=1,sep="\t", check.names = FALSE)
pmm_counts <- pmm_counts[, !str_detect(colnames(pmm_counts), "^CD19\\+")]
map <- c(
  "CD19-GSM2171212" = "BMPC4",
  "CD19-GSM2171213" = "BMPC5",
  "CD19-GSM2171214" = "BMPC6",
  "CD19-GSM2171215" = "BMPC7")

old_names <- colnames(pmm_counts)

new_names <- old_names

for (old_id in names(map)) {
  new_names <- str_replace(
    new_names,
    paste0("^", old_id),
    map[[old_id]]
  )}

colnames(pmm_counts) <- new_names


# get mm+PMM count matrix and remove redundant geneSymbols
bulk_counts <- 
  mm_counts %>% 
  tibble::rownames_to_column(var = "GeneSymbol") %>% 
  left_join(tibble::rownames_to_column(pmm_counts, var = "GeneSymbol"), 
            by = "GeneSymbol") %>% 
  dplyr::mutate(GeneSymbol = str_replace_all(GeneSymbol, "ENSG..+_", ""))
# sum duplicated genes
bulk_counts_sum <- rowsum(
  as.matrix(bulk_counts %>% dplyr::select(-GeneSymbol)),
  group = bulk_counts$GeneSymbol,
  reorder = FALSE)

bulk_counts_sum <- as.data.frame(bulk_counts_sum) %>%
  tibble::rownames_to_column("GeneSymbol") %>% 
  dplyr::filter(GeneSymbol %in% ref_counts$GeneSymbol)

write.table(bulk_counts_sum,
  file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/cibersort/mm_pmm_counts/mm_pmm_raw_counts.txt",
  sep = "\t",
  quote = FALSE, row.names = FALSE,
  col.names = TRUE)

# 800 sample too much, split into 4 batches ####################################
batch_size <- 200

sample_batches <- split(
  colnames(bulk_counts_sum)[-1],
  ceiling(seq_along(colnames(bulk_counts_sum)[-1]) / batch_size)
)

purrr::iwalk(sample_batches, function(samples, i) {
  bulk_counts_sum[, c("GeneSymbol", samples), drop = FALSE] %>%
    write.table(
      file = file.path("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/cibersort/mm_pmm_counts", paste0("PMM_MM_bulk_counts_batch_", i, ".txt")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
})


# CIBERSORTx cell proportion prediction #######################
# read in the CIBERSORTx result 
mm_batch1_prop <- read_csv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/cibersort/proportion/CIBERSORTx_Job19_Results_batch1_200.csv")
mm_batch2_prop <- read_csv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/cibersort/proportion/CIBERSORTx_Job20_Results_batch2_200.csv")
mm_batch3_prop <- read_csv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/cibersort/proportion/CIBERSORTx_Job21_Results_batch3_200.csv")
mm_batch4_prop <- read_csv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/cibersort/proportion/CIBERSORTx_Job22_Results_batch4_200.csv")

# manipulate, plot and save the final cell type proportion result
mm_pmm_prop <- rbind(mm_batch1_prop, mm_batch2_prop, 
                     mm_batch3_prop, mm_batch4_prop) %>% 
  dplyr::mutate(group = ifelse(str_detect(Mixture, "SRR"), "MM", 
                               ifelse(str_detect(Mixture, "PMM"), "PMM", "BMPC")))
write_tsv(mm_pmm_prop, file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/cibersort/proportion/estimated_mm_pmm_prop.tsv")

cell_prop <- mm_pmm_prop %>% 
  dplyr::select(-`P-value`, -Correlation, -RMSE) %>% 
  pivot_longer(
    cols = -c(Mixture, group),
    names_to = "cell_type",
    values_to = "Proportion") %>% 
  ggplot(aes(x=cell_type, y=Proportion, fill = group))+
  geom_boxplot(outlier.shape = NA)+
  ggpubr::theme_pubr(base_size = 9.6)+
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  labs(
    x = "Cell type",
    y = "Estimated proportion"
  )

ggsave(paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/cibersort/proportion", "/estimated_cell_proportion.pdf"), 
       plot = cell_prop, width = 15, height = 10, units = "cm", 
       limitsize = FALSE)

# perform analysis for cell type proportion and update DE for MM
## L1 number and plasma cell proportion correlation
mm_pmm_l1num <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/cibersort/L1_num_cor/mm_pmm_L1_num.tsv")

plot_df <- mm_pmm_prop %>% 
  dplyr::left_join(mm_pmm_l1num, by = c("Mixture" = "sampleID")) %>% 
  dplyr::mutate(group = group.y) %>% 
  dplyr::mutate(
    group = factor(
      group,
      levels = c("BMPC", "MM (Q1)", "MM (Q2)", "MM (Q3)", "MM (Q4)", "PMM")
    )
  )

my_comparisons <- combn(
    levels(plot_df$group),
    2,
    simplify = FALSE
  )

plasma_MMPMM_plot <- ggplot(plot_df, aes(x = group, y = Plasma_Cell, fill = group)) +
  geom_boxplot(
    outlier.size = 0.4
  ) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     method.args = list(exact = FALSE),
                     # method = "t.test", 
                     p.adjust.method = "BH",
                     size = 2.5,
                     label = "p.signif")+
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.35))
  ) +
  guides(fill = "none") +
  theme_pubr(base_size = 9.6, legend = c(0.8, 0.3)) +  # Set base font size
  theme(
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.key.size = unit(0.4, "lines"),     # make the color square smaller
    legend.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9, face = "bold")  # or just size = 9
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Estimated plasma cell proportion"
  )

ggsave(paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/cibersort/L1_num_cor", "/estimatedPlasmaProportion_MMPMM_plot_forSUPP.pdf"), 
       plot = , width = 10, height = 9, units = "cm", 
       limitsize = FALSE)

## DE after adjusting for plasma cell proportion



## plasma cell proportion is significantly enriched in Q1-Q4?
### see /home/qwan/githubRepo/coh_PMM/RNAseq_scripts/MMRF_analysis/mmrf_subtypes/mm_meta_analysis.R







