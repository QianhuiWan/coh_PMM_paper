
# load pkg and data needed
library(tidyverse)
library(magrittr)
library(bsseq)
library(methrix)
library(plyranges)
library(BiocParallel)
library(mixOmics)
library(gridExtra)
library(biomaRt)

input_dir <- '/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/methylKit_output'
outDir <- '/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/methylKit_output/DMR_res'
KZFP_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs"
bsseq_obj <- readRDS("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/methylKit_output/bismark_methylKit_bsseq_obj.rds")


# 1. get DMR DNAm mean for each sample
# step1: load DMR and CpGislands
dmrs_raw <- readRDS(paste0(outDir, "/bismark_methylKit_bsseq_DMRseq_dmrs.rds"))
dmrs <- dmrs_raw %>% 
  plyranges::filter(qval<0.05) %>% 
  plyranges::mutate(hyper_hypo = ifelse(stat >0, "hyper", "hypo"))

# also want to have all DNAm GR
meth_mat <- getMeth(bsseq_obj, type = "raw")  # type = "raw" gives proportions
meth_mat[is.nan(meth_mat)] <- NA
dnam_GR <- granges(bsseq_obj)  # returns a GRanges object

# attach methylation levels as metadata columns
for (i in seq_len(ncol(meth_mat))) {
  mcols(dnam_GR)[[colnames(meth_mat)[i]]] <- meth_mat[, i]
}


# get KZFP promoter overlapped CpG islands
## load KZFP coordinates  #378
kzfp_ranges <- read_tsv(file = paste0(KZFP_dir, "/kzfps_hg38_all_geneCoordinates.tsv")) %>% 
  as_granges() %>% 
  plyranges::filter(seqnames %in%  c(paste0(c(1:22, "X", "Y")))) %>% 
  as.data.frame() %>% 
  plyranges::mutate(seqnames = paste0("chr", seqnames)) %>% 
  as_granges()

## load KZFP promoters #
promoter_ranges <- read_tsv(file = paste0(KZFP_dir, "/kzfps_hg38_all_TSS2000.tsv")) %>% 
  as_granges() %>% 
  plyranges::filter(seqnames %in%  c(paste0(c(1:22, "X", "Y")))) %>% 
  as.data.frame() %>% 
  plyranges::mutate(seqnames = paste0("chr", seqnames)) %>% 
  as_granges() 

## load CpG islands (unmasked)
# cpgIsland_GR <- rtracklayer::import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/cpgIsland/UCSC_CpGislands_unmask.bed") %>%
#   plyranges::filter(seqnames %in%  c(paste0("chr", c(1:22, "X", "Y")))) %>% 
#   plyranges::mutate(cpgIsland_id = paste("cpg", 1:length(.), sep = "_"))

## load CpG islands (masked)
cpgIsland_GR <- rtracklayer::import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/cpgIsland/UCSC_CpGislands.bed") %>%
  plyranges::filter(seqnames %in%  c(paste0("chr", c(1:22, "X", "Y")))) %>% 
  plyranges::mutate(cpgIsland_id = paste("cpg", 1:length(.), sep = "_"))

# get KZFP promoter overlapped CpG islands
kzfp_promoter_CpGisland <- find_overlaps(cpgIsland_GR, promoter_ranges) #325


# 2. get gene expression for each sample

## PMM gene expression
### load sample sheet with L1 subfam numbers
PMM_ss <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/PMM_sampleSheet_R2.tsv")
PMM_ss$group <- str_replace_all(PMM_ss$group, "BMPC", "control")
PMM_ss_L1PAnum <- PMM_ss %>% 
  dplyr::select(sampleName, L1PA2, L1PA3) %>% 
  dplyr::mutate(sampleName = str_replace_all(sampleName, "BMPC", "C")) %>% 
  column_to_rownames(var = "sampleName")

### load lcpm PMM
lcpm_pmm <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/lcpm_TMM_PMM_STARcounts_R2.tsv")

### down-regulated KZFPs
topTable_dat_PMM <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/topTable_res_R2.tsv")

normalized_lcpm_fil_pmmDownKZFPtb <- lcpm_pmm %>% 
  dplyr::mutate(gene_id_version = str_replace_all(gene_name, ".+_", "")) %>% 
  dplyr::mutate(gene_id = str_remove(gene_id_version, "\\..*")) %>% 
  dplyr::mutate(gene_symbol = str_replace_all(gene_name, "_.+", "")) %>% 
  dplyr::select(-gene_id_version, -gene_symbol, -gene_name) %>% 
  dplyr::left_join(topTable_dat_PMM, by = "gene_id") %>% 
  # dplyr::filter(FDR<0.05 & logFC< -0.5) %>%
  dplyr::filter(logFC< -0.5) %>%
  dplyr::filter(gene_id_version %in% kzfp_ranges$gene_id) %>% 
  dplyr::mutate(lcpm = ifelse(lcpm < 0, 0, lcpm))

pmm_down_KZFP_ids <- normalized_lcpm_fil_pmmDownKZFPtb$gene_id_version %>%
  unique() #24. 96

pmm_lcpm_SD_mean <- lcpm_pmm %>%
    dplyr::mutate(gene_id_version = str_replace_all(gene_name, ".+_", "")) %>%
  dplyr::mutate(gene_id = str_remove(gene_id_version, "\\..*")) %>%
  dplyr::mutate(gene_symbol = str_replace_all(gene_name, "_.+", "")) %>%
  dplyr::filter(gene_id_version %in% pmm_down_KZFP_ids) %>%
  # dplyr::arrange(match(gene_name, ZNF_geneAnno$gene_name)) %>%
  mutate(PMM_lcpm = ifelse(lcpm < 0, 0, lcpm)) %>%
  dplyr::filter(str_detect(sample, "PMM")) %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarise(PMM_SD_lcpm = sd(PMM_lcpm, na.rm =TRUE),
                   PMM_mean_lcpm = mean(PMM_lcpm, na.rm =TRUE),
                   PMM_SD_mean_lcpm = PMM_SD_lcpm/PMM_mean_lcpm)

normalized_lcpm_fil_pmmDownKZFPtb_fil <- normalized_lcpm_fil_pmmDownKZFPtb %>%
  dplyr::left_join(pmm_lcpm_SD_mean, by = "gene_symbol") %>%
  dplyr::filter(PMM_SD_mean_lcpm > 0.1)

pmm_down_KZFP_ids <- normalized_lcpm_fil_pmmDownKZFPtb_fil$gene_id_version %>%
  unique() #72

## MM gene expression
# load topTable output
topTable_dat <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_redo_res/lmVoom_topTable_Q4Q1_dat.tsv") %>% 
  dplyr::select(gene_id, logFC, FDR)
lcpm_mm <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s1_checkStrandness_R4/lcpm_TMM_764_STARcounts_cpmFilterApplied.tsv")

# load TMM lcpm matrix and get DE genes lcpm
normalized_lcpm_fil_mmDEtb <- lcpm_mm %>% 
  dplyr::mutate(gene_id_version = str_replace_all(gene_name, "_.+", "")) %>% 
  dplyr::mutate(gene_id = str_remove(gene_id_version, "\\..*")) %>% 
  dplyr::mutate(gene_symbol = str_replace_all(gene_name, ".+_", "")) %>% 
  dplyr::left_join(topTable_dat, by = "gene_id") %>% 
  dplyr::filter(FDR<0.05 & abs(logFC)>1) %>%
  dplyr::mutate(lcpm = ifelse(lcpm < 0, 0, lcpm))

mm_DE_gene_ids <- normalized_lcpm_fil_mmDEtb$gene_id_version %>% unique() #720


# get MM DE gene overlapped CpG islands
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                   # use when https://www.ensembl.org is not working
                   host = "https://nov2020.archive.ensembl.org") 
# get all KZFP GRanges
mm_DE_genes <- getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version", 
                                   "external_gene_name", "description",
                                   "chromosome_name",
                                   "start_position", "end_position", "strand"),
                    filters = "ensembl_gene_id_version",
                    values = mm_DE_gene_ids,  # Example: genes with "ZNF" or "KZFP" in their names
                    mart = ensembl)

# Convert to GRanges
mm_DE_gene_ranges <- GRanges(
  seqnames = mm_DE_genes$chromosome_name,
  ranges = IRanges(start = mm_DE_genes$start_position, end = mm_DE_genes$end_position),
  strand = ifelse(mm_DE_genes$strand == 1, "+", "-"),
  gene_name = mm_DE_genes$external_gene_name,
  gene_id = mm_DE_genes$ensembl_gene_id_version,
  description = mm_DE_genes$description
)

# Define promoter regions around the TSS
mm_DE_gene_promoter_ranges <- 
  promoters(mm_DE_gene_ranges, upstream = 2000, downstream = 2000) %>% 
  plyranges::filter(seqnames %in%  c(paste0(c(1:22, "X", "Y")))) %>% 
  as.data.frame() %>% 
  plyranges::mutate(seqnames = paste0("chr", seqnames)) %>% 
  as_granges() 

mm_DE_gene_promoter_CpGisland <- find_overlaps(cpgIsland_GR, mm_DE_gene_promoter_ranges) #260



# 3. plot cor
# cor between DMR overlapped KZFP promoter overlapped CpG islands and KZFP gene expression
## need matrix list: lcpm, matched DNAm and L1PA2, L1PA3 numbers

KZFP_lcpm <- lcpm_pmm %>% 
  dplyr::mutate(gene_name_original = gene_name) %>% 
  dplyr::mutate(gene_name = str_replace_all(gene_name_original, "_.+", "")) %>% 
  dplyr::mutate(gene_id_version = str_replace_all(gene_name_original, ".+_", "")) %>% 
  dplyr::filter(gene_id_version %in% pmm_down_KZFP_ids) %>%
  # dplyr::filter(gene_name %in% c("ZNF425", "ZNF141", "ZNF382","ZNF680", "ZFP28")) %>% 
  # mutate(lcpm = ifelse(lcpm < 0, 0, lcpm)) %>%
  pivot_wider(names_from = gene_name, values_from = lcpm,
              id_cols = c(sample)) %>% 
  dplyr::mutate(sample = str_replace_all(sample, "BMPC", "C")) %>% 
  arrange(sample, rownames(PMM_ss_L1PAnum)) %>% 
  column_to_rownames(var = "sample")

zero_variance_cols <- sapply(KZFP_lcpm, function(x) var(x, na.rm = TRUE) == 0)
# View columns with zero variance
names(zero_variance_cols[zero_variance_cols])
# Remove columns with zero variance
KZFP_lcpm <- KZFP_lcpm[, !zero_variance_cols]

### load CpG meth and get mean DNAm for each DMR at KZFP promoter regions
nms <- colnames(values(dnam_GR))
cols <- lapply(nms, function(.) { col <- ensym(.); quo(mean(!!col, na.rm=TRUE)) })

cols <- lapply(nms, function(col_name) {
  col_sym <- sym(col_name)
  quo(mean(!!col_sym, na.rm = TRUE))
})
# the column names for the summarize op
names(cols) <- paste0("mean_", nms)

dmr_KZFPpromoterCGI_df <- kzfp_promoter_CpGisland %>%  
  dplyr::filter(gene_id %in% pmm_down_KZFP_ids) %>% 
  dplyr::filter(gene_name %in% colnames(KZFP_lcpm)) %>% 
  find_overlaps(dnam_GR) %>% 
  plyranges::group_by(gene_name) %>% 
  plyranges::summarise(!!! cols) %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "gene_name") %>% 
  t() %>% 
  `rownames<-`(str_replace_all(rownames(.), "mean_", "")) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  dplyr::mutate(sample = str_replace_all(sample, "_rest", "")) %>% 
  dplyr::mutate(sample = str_replace_all(sample, "B", "C")) %>% 
  dplyr::filter(sample %in%  rownames(PMM_ss_L1PAnum)) %>% 
  arrange(sample, rownames(PMM_ss_L1PAnum)) %>% 
  column_to_rownames(var = "sample")

# just plot DNAm against gene expression for down-regulated KZFPs
# KZFP_lcpm_forplot <- KZFP_lcpm[, colnames(KZFP_lcpm) %in% colnames(dmr_KZFPpromoterCGI_df)]
KZFP_lcpm_forplot <- 
  KZFP_lcpm[, match(colnames(dmr_KZFPpromoterCGI_df), colnames(KZFP_lcpm))] %>% 
  dplyr::mutate(expression = rep("lcpm", nrow(.))) %>% 
  rownames_to_column("sample")

KZFPpromoterCGI_df_forplot <- dmr_KZFPpromoterCGI_df %>% 
  # dplyr::select(ZNF141,ZNF382) %>% 
  dplyr::mutate(expression = rep("DNAm", nrow(.))) %>% 
  rownames_to_column("sample") %>%
  rbind(KZFP_lcpm_forplot) %>% 
  pivot_longer(
    cols = -c(sample, expression),
    names_to = "KZFP",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = expression,
    values_from = value
  ) %>% 
  dplyr::mutate(DNAm = ifelse(is.nan(DNAm), NA, DNAm))

KZFPpromoterCGI_DNAm_lcpm_plot <- 
  KZFPpromoterCGI_df_forplot %>% 
  dplyr::filter(!str_detect(sample, "C")) %>%
  ggplot(aes(x = lcpm, y = DNAm)) + 
  geom_point() + geom_smooth(method = "lm")+ 
  facet_wrap(~sample, scales = "free_x")+
  # ylim(c(-20,150))+
  ggpubr::stat_cor(label.y=c(0.8,0.8), method = "pearson",size=3) +
  # ggrepel::geom_text_repel(aes(label = sample), max.overlaps = 100,
  #                          box.padding = 0.5, size = 3)+
  ggpubr::theme_pubr(base_size = 9.6)


# ggsave(paste0(outDir, "/sigDown_KZFPpromoterCGI_DNAm_lcpm_Cor_facetSample.pdf"), 
ggsave(paste0(outDir, "/KZFPpromoterCGI_DNAm_lcpm_Cor_facetSample.pdf"), 
       plot = KZFPpromoterCGI_DNAm_lcpm_plot,
       width = 16, height = 12, units = "cm", limitsize = FALSE)


KZFPpromoterCGI_DNAm_lcpm_plot <- 
  KZFPpromoterCGI_df_forplot %>% 
  dplyr::filter(!str_detect(sample, "C")) %>%
  ggplot(aes(x = lcpm, y = DNAm)) + 
  geom_point() + geom_smooth(method = "lm")+ 
  # facet_wrap(~sample, scales = "free_x")+
  # ylim(c(-20,150))+
  ggpubr::stat_cor(label.y=c(0.8,0.8), method = "pearson",size=3) +
  # ggrepel::geom_text_repel(aes(label = sample), max.overlaps = 100,
  #                          box.padding = 0.5, size = 3)+
  ggpubr::theme_pubr(base_size = 9.6)


# ggsave(paste0(outDir, "/sigDown_KZFPpromoterCGI_DNAm_lcpm_Cor_allSample.pdf"), 
ggsave(paste0(outDir, "/KZFPpromoterCGI_DNAm_lcpm_Cor_allSample.pdf"), 
       plot = KZFPpromoterCGI_DNAm_lcpm_plot,
       width = 16, height = 8, units = "cm", limitsize = FALSE)


# plot heatmap for MM DE genes: DE gene promoter DMR DNAm
## generate matrix: rows as DE genes with promoter CGI, columns are MM samples
nms <- colnames(values(dnam_GR))
cols <- lapply(nms, function(.) { col <- ensym(.); quo(mean(!!col, na.rm=TRUE)) })

cols <- lapply(nms, function(col_name) {
  col_sym <- sym(col_name)
  quo(mean(!!col_sym, na.rm = TRUE))
})
# the column names for the summarize op
names(cols) <- paste0("mean_", nms)

dmr_mm_DE_promoterCGI_mat <- mm_DE_gene_promoter_CpGisland %>% 
  find_overlaps(dnam_GR) %>% 
  plyranges::group_by(gene_name) %>% 
  plyranges::summarise(!!! cols) %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "gene_name") %>% 
  `colnames<-`(str_replace_all(colnames(.), "mean_", "")) %>% 
  as.matrix()

dmr_mm_DE_promoterCGI_mat[is.nan(dmr_mm_DE_promoterCGI_mat)] <- 0


library(ComplexHeatmap)
library(circlize)

# Set color scale from 0 (white) to 1 (blue)
col_fun <- colorRamp2(c(0, 1), c("white", "blue"))

pdf(file = paste0(outDir, "/mm_DE_DNAm_heatmap_CH.pdf"),
    width = 16/2.54, height = 10/2.54)
# Plot heatmap
Heatmap(
  dmr_mm_DE_promoterCGI_mat,
  show_row_names = FALSE,
  clustering_distance_columns = "manhattan", 
  clustering_method_columns = "complete",
  col = col_fun,
  heatmap_legend_param = list(
    title = "DNAm",
    title_gp = gpar(fontsize = 9.6, fontface = "bold"),  # optional styling
    labels_gp = gpar(fontsize = 10)                     # optional styling
  )
)
dev.off()


## complexheatmap









# Use DIABLO, i.e. plsda
# KZFP_lcpm <- KZFP_lcpm[colnames(KZFP_lcpm) %in% 
#                          c("ZNF425", "ZNF141", "ZNF382","ZNF680", "ZFP28")]

X <- list(mRNA = as.matrix(KZFP_lcpm), 
          DNAm = as.matrix(dmr_KZFPpromoterCGI_df),
          TE = as.matrix(PMM_ss_L1PAnum))

MyResult.pca <- spca(X$mRNA, ncomp = 2, scale = TRUE)
plotIndiv(MyResult.pca, group = PMM_ss$group, ind.names = FALSE,
          legend = TRUE, 
          title = 'KZFPs, PCA comp 1 - 2')

Y <- factor(PMM_ss$group, levels = c("control", "PMM"))
summary(Y)

design <- matrix(0.1, ncol = length(X), nrow = length(X), 
                 dimnames = list(names(X), names(X)))
diag(design) <- 0
design <- cbind(design, Y = c(1, 1, 1))
design <- rbind(design, Y = c(1, 1, 1, 0))
diag(design) <- 0

# check cor. between comp1 in RNA and DNAm matrix
pls.res <- pls(X$mRNA, X$DNAm, ncomp = 1)
cor(pls.res$variates$X, pls.res$variates$Y)


# optimise the number of PC to use
diablo.res <- block.splsda(X, Y, ncomp = 2, design = design)

selectVar(diablo.res, block = 'mRNA', comp = 1)

p2 <- plotDiablo(diablo.res, ncomp = 1)

pdf(paste0(outDir, "/kzfp_dataIntegration_cor.pdf"))
plotDiablo(diablo.res, ncomp = 1)
dev.off()


p3 <- circosPlot(diablo.res, cutoff = 0.7, line = TRUE, 
                 color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
                 color.cor = c("chocolate3","grey20"), size.labels = 1.5)


pdf(paste0(outDir, "/kzfp_dataIntegration_circos.pdf"))
circosPlot(diablo.res, cutoff = 0.7, line = TRUE, 
           color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), size.labels = 2, 
           size.variables = 0.6)
dev.off()

network(diablo.res, blocks = c(1,2,3), 
        cutoff = 0.4,
        color.node = c('darkorchid', 'brown1', 'lightgreen'),
        # To save the plot, uncomment below line
        #save = 'png', name.save = 'diablo-network'
)







