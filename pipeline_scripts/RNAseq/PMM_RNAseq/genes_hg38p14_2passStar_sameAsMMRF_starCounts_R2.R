
# R2: featureCount L1 number used in this round

# .libPaths()
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# use count matrix from feature counts
# get gene counts stats, PMM
dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe"
output_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts"

# laod packages
library(tidyverse)
library(magrittr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(rtracklayer)
library(ggpubr)
library(edgeR)
library(limma)
library(msigdbr)

# getNum function
getNum <- function(path, pattern){
  # get path and sample names
  counts_path <- list.files(path, recursive = TRUE, pattern = pattern, 
                            full.names = TRUE)
  sample_names <- str_remove_all(basename(counts_path), pattern = "_..+")
  
  # get ref. genome version
  refGenome = "hg38"
  
  # load data
  df_list <- list()
  for (i in 1:length(sample_names)) {
    if (str_detect(counts_path[i], "TEs.bed")) {
      # load TEs.bed including FPKM and TPM
      data <- read.table(counts_path[i], 
                         header = FALSE, sep = "\t", 
                         stringsAsFactors = FALSE)
      colnames(data) <- c("seqnames", "start", "end", 
                          "gene_id", "transcript_id", "strand", "repeat", 
                          "cov", "FPKM", "TPM")
      data$cov = as.numeric(data$cov)
      data$FPKM = as.numeric(data$FPKM)
      data$TPM = as.numeric(data$TPM)
      
    }else if (str_detect(counts_path[i], "ballgown.gtf")) {
      # load tx.gtf including FPKM and TPM
      data <- import(counts_path[i])
      data <- as.data.frame(data, stringsAsFactors = FALSE)
      data <- data[data$type == "transcript",]
      data$cov = as.numeric(data$cov)
      data$FPKM = as.numeric(data$FPKM)
      data$TPM = as.numeric(data$TPM)
      
    }else if (str_detect(counts_path[i], "alignedReads.txt")) {
      # load aligned reads number
      data <- read.table(file = counts_path[i], stringsAsFactors = FALSE)
      colnames(data) <- "num_aligned_reads"
    }else if (str_detect(counts_path[i], "gene_counts.txt")) {
      # load aligned reads number
      data <- read.table(file = counts_path[i], stringsAsFactors = FALSE, 
                         header = TRUE)
      # colnames(data) <- c("Geneid", "Chr","Start", "End", "Strand",
      #                     "Length", "crick", "watson")
      # data <- data %>% dplyr::mutate(totalCount = crick + watson)
      colnames(data) <- c("Geneid", "Chr","Start", "End", "Strand",
                          "Length", "totalCount")
    }else if (str_detect(counts_path[i], "_ReadsPerGene.out.tab")){
      # load aligned reads number
      data <- read.table(file = counts_path[i], stringsAsFactors = FALSE, 
                         header = FALSE, skip = 4)
      colnames(data) <- c("Geneid", "totalCount", "watson", "crick")
    }else if (str_detect(counts_path[i], "transcript_count_matrix.csv")) {
      # load csv files
      data <- read_csv(counts_path[i])
      
    }else if (str_detect(counts_path[i], "gene_count_matrix.csv")) {
      data <- read_csv(counts_path[i])
      
    }else if (str_detect(counts_path[i], "TEtxs_count_matrix.csv")) {
      data <- read.table(file = counts_path[i], stringsAsFactors = FALSE)
      colnames(data) <- c("transcript_id", "counts", "repeat")
    }
    
    # add data df to df_list
    df_list[[i]] <- as.data.frame(data)
  }
  
  # name df_list
  names(df_list) <- paste0(sample_names, "_", refGenome)
  
  # retun the df_list
  return(df_list)
}


## use getNum FUN to get all transcripts #######################################

df_list <- getNum(path = paste0(dir),
                  pattern = "_ReadsPerGene.out.tab")
# saveRDS(df_list, file =  paste0(output_dir, "/allSample_genes_list_STARcounts.rds"))
# df_list <- readRDS(file =  paste0(output_dir, "/allSample_genes_list_STARcounts.rds"))

# bind list for gene numbers
## only keep total count > 0 and plot the library size for rnaseq samples
g_count_list <- lapply(df_list, function(x) {nrow(x[x[,2]>0,])})
num_totalG <- do.call("rbind", g_count_list)
num_totalG <- as.data.frame(num_totalG, stringsAsFactors = FALSE)
colnames(num_totalG) <- "num_G"

num_totalG %>% rownames_to_column(var = "samples") %>% 
  ggplot(aes(x=samples, y=num_G))+
  geom_bar(stat = 'identity')+
  labs(title = "detected genes", x="", y = "number of total genes")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# bind list for gene matrix
## gene list to df and add id column
gene_df <- df_list %>% 
  bind_rows(.id = "sampleID")
## save all gene counts to tsv
readr::write_tsv(gene_df,
                 file = paste0(output_dir, "/genes_starCounts_PMM_R2.tsv"))

# load gtf annotation files for human genes
data <- rtracklayer::import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/teAnno_round3/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf")
data <- as.data.frame(data, stringsAsFactors = FALSE)
data <- data[data$type == "transcript",]

# get unique gene_id and gene_name
genes_anno <- data %>% 
  dplyr::filter(gene_id %in% unique(gene_df$Geneid)) %>% 
  dplyr::select(gene_id, gene_name) %>% 
  unique()


# gene counts analysis
count_matrix <- gene_df %>%
  dplyr::left_join(genes_anno, by = c("Geneid"="gene_id")) %>% 
  # dplyr::mutate(Geneid = paste0(gene_name, "_", Geneid)) %>% 
  dplyr::mutate(Geneid = paste0(Geneid, "_", gene_name)) %>% 
  # dplyr::select(-watson, -crick) %>% #unique() %>% 
  pivot_wider(names_from = sampleID, values_from = totalCount,
              id_cols = Geneid, values_fill = 0) %>% 
  column_to_rownames(var="Geneid") %>% 
  as.matrix()
# count_matrix[is.na(count_matrix)] <- 0
colnames(count_matrix) <- str_replace_all(colnames(count_matrix), "_..+", "")

# Save as a .tsv file
write.table(count_matrix, 
            file =  paste0(output_dir, "/gene_count_matrix_STARcounts_PMM_R2.tsv"), 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


# get sample meta data
## load L1 number info
L1_tetx_num <- read_tsv(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/activationScore_R2/filtered_L1/num_active_all_L1s_L1PA2_L1PA3_pmm.tsv")
L1_tetx_num <- L1_tetx_num %>% 
  dplyr::select(sampleID, L1, L1PA2, L1PA3) %>% unique()

PMM_sampleSheet <- 
  colnames(count_matrix) %>% 
  as.data.frame() %>%
  `colnames<-`(c("sampleName")) %>% 
  dplyr::mutate(group=str_replace_all(sampleName, "[:digit:]", "")) %>% 
  left_join(L1_tetx_num, by = c("sampleName"="sampleID"))

PMM_sampleSheet %>% write_tsv(paste0(output_dir, "/PMM_sampleSheet_R2.tsv"))



# normalize counts with TMM method
count_matrix <- read.table(paste0(output_dir, "/gene_count_matrix_STARcounts_PMM.tsv"))
PMM_sampleSheet <- read_tsv(paste0(output_dir, "/PMM_sampleSheet.tsv"))
table(colnames(count_matrix)==PMM_sampleSheet$sampleName)

dge <- DGEList(counts = count_matrix, group = PMM_sampleSheet$group)

# Define a threshold for filtering
keep <- rowSums(cpm(dge, normalized.lib.sizes = TRUE) >=1 ) >= 3
# keep <- filterByExpr(dge, group=dge$samples$group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
normalized_lcpm <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)

lcpm_df <- normalized_lcpm %>% as.data.frame() %>% 
  rownames_to_column(var = "gene_name") %>% 
  pivot_longer(cols = c(starts_with("BMPC"), starts_with("PMM")),
               names_to = "sample", 
               values_to = "lcpm") %>% 
  mutate(lcpm = ifelse(lcpm < 0, 0, lcpm))

readr::write_tsv(lcpm_df, file = paste0(output_dir,"/lcpm_TMM_PMM_STARcounts_R2.tsv"))


# DE 
## design matrix
design <- model.matrix(~0 + group, data=PMM_sampleSheet)
## make contrast
con <- makeContrasts(group = groupPMM-groupBMPC,
                     levels=design) 
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

## get topTable output
dat <- topTable(v, coef = "group", number=Inf) %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  dplyr::mutate(gene_id_version = str_replace_all(gene_name, ".+_", "")) %>% 
  dplyr::mutate(gene_id = str_remove(gene_id_version, "\\..*")) %>% 
  dplyr::mutate(gene_symbol = str_replace_all(gene_name, "_.+", "")) %>% 
  # distinct(gene_name, .keep_all = TRUE) %>%
  dplyr::select(-t, -B, FDR = `adj.P.Val`) %>%
  mutate(Sig = ifelse(FDR < 0.05 & logFC > 1, "Up", 
                      ifelse(FDR < 0.05 & logFC < -1, "Down", "NS"))) 

write_tsv(dat, file = paste0(output_dir,"/topTable_res_R2.tsv"))

## GSEA
library(msigdbr)
as.data.frame(msigdbr_collections())
m_df = msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- m_df[,c("gs_name", "ensembl_gene")] # "human_gene_symbol",

### Create a ranked list
ranked_list <- dat %>%
  dplyr::select(gene_id, logFC) %>% 
  arrange(desc(logFC)) %>%
  deframe()  # Convert to named vector (gene names as names, logFC as values)

## rank result by p values
dum <- dat$FDR
## Avoid Inf when taking log of 0 
# dum[dum == 0] <- 1e-300
score = -log10(dum)
#It is important to differentiate up and down regulation
score <- ifelse(dat$logFC >= 0, score, -score)
names(score) <- dat$gene_id
score <- sort(score, decreasing = TRUE)
## Add jitter to avoid ties
set.seed(52)
score <- score + rnorm(length(score), mean = 0, sd = 1e-6)
score <- sort(score, decreasing = TRUE)

## perform GSEA
library(clusterProfiler)
em_1 <- GSEA(ranked_list, TERM2GENE=hallmark, eps = 0, 
           pvalueCutoff=0.05, pAdjustMethod="BH")

em_2 <- GSEA(score, TERM2GENE=hallmark, eps = 0, 
           pvalueCutoff=0.05, pAdjustMethod="BH")

# Create a dotplot of top 30 pathways
dotplot(em_1, orderBy="NES", x="NES", showCategory = 20)
dotplot(em_2, orderBy="NES", x="NES", showCategory = 30)

dotplot(em_2, orderBy = "NES", x = "NES", showCategory = 25) +
  # coord_flip() +
  scale_color_gradient(low = "red", high = "blue") +
  scale_size(range = c(2, 10)) +
  labs(
    title = "Top Enriched Pathways",
    x = "Normalized Enrichment Score (NES)",
    y = "",
    color = "Adjusted p-value",
    size = "Gene Count"
  ) +
  theme_pubr(base_size = 12, legend = "right")+
  theme(axis.text.y = element_text(size = 9.6))


ggsave(paste0(output_dir, "/lmVoomLfc_gsea_NES_top20.pdf"), width = 8, height = 12)




# plot lcpm of some genes:
p1 <- 
  lcpm_df %>% 
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>% 
  dplyr::filter(gene_name %in% c("ZNF141", "MPHOSPH8", 
                                 "ATF7IP", "SIRT1",
                                 "TASOR", "MORC2", "PPHLN1", "SETDB1", 
                                 "TASOR2", "IRF2",
                                 "DNMT1", "DNMT3A")) %>% 
  # distinct() %>%
  dplyr::left_join(PMM_sampleSheet, by = c("sample"="sampleName")) %>% 
  # ggplot(aes(x=sample_index, y=lcpm, fill=gene_name))+
  # ggplot(aes(x=reorder(sample, L1PA2), y=lcpm))+
  ggplot(aes(x=L1PA2, y=lcpm))+
  geom_point()+
  # geom_bar(stat = "identity")+
  geom_smooth(method = "lm")+
  ggpubr::stat_cor(method = "pearson")+
  # geom_bar(stat = "identity", position = "dodge")+
  labs(title = "HUSH complex genes", x="")+
  facet_wrap(~gene_name)+
  # facet_wrap(~gene_name+transcript_id, scales = "free_y")+
  ggpubr::theme_pubr()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(paste0(output_dir, "/lcpm_HUSHgenes_R2_STARcounts.pdf"), plot = p1,
       width = 25, height = 20, units = "cm", limitsize = FALSE)


# Create a data frame with genes and their annotations
key_gene_df <- data.frame(
  Gene = c("ZNF141", "ZNF382", "ZNF429", "ZNF483", "ZNF860", 
           "NFKB1", "JAK1", "IFIH1", "IRF1",
           "SNAI2", "TWIST1",
           "E2F1", "E2F2", "E2F5", "E2F7", "E2F8", "RB1", "CDK1",
           "SUZ12", "EZH2", "PRC1", "EED",
           "TP53", "MDM2"),
  Category = c("Young L1 ZNF", "Young L1 ZNF", "Young L1 ZNF", 
               "Stem Cell ZNF", "Sig.Down ZNF in PMM", 
               "NFKB", "NFKB", "IFN", "IFN",
               "EMT", "EMT",
               "E2F", "E2F", "E2F", "E2F", "E2F", "E2F", "E2F",
               "PRC2", "PRC2", "PRC2", "PRC2",
               "TP53 Pathway", "TP53 Pathway")
)

#NFKB(NFKB -> proinflammatory cytokines), IFN(IRF3 -> type I IFNs)
# PRC2, EZH2 bind to young L1s
p2 <- 
  lcpm_df %>% 
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>% 
  dplyr::filter(gene_name %in% key_gene_df$Gene) %>% 
  dplyr::left_join(key_gene_df, by = c("gene_name"="Gene")) %>% 
  dplyr::left_join(PMM_sampleSheet, by = c("sample"="sampleName")) %>% 
  # ggplot(aes(x=sample_index, y=lcpm, fill=gene_name))+
  ggplot(aes(x=L1PA2, y=lcpm))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson")+
  labs(title = "key genes (ZNF, NFKB, EMT, E2F, PRC and TP53)", x="")+
  facet_wrap(~Category+gene_name, scales = "free_y", ncol = 5)+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(paste0(output_dir, "/lcpm_keygenes_R2_STARcounts.pdf"), plot = p2,
       width = 25, height = 35, units = "cm", limitsize = FALSE)


