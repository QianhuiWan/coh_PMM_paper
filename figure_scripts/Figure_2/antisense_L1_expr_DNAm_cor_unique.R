library(tidyverse)
library(edgeR)
library(GenomicRanges)
library(plyranges)

# get antisense L1 up 10kb DNAm
mean_DNAm_files <- list.files("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/methylKit_output/region_mean", 
                              pattern = "DNAm.tsv", full.names = TRUE)

mean_DNAm_sample_ls <- list()
for (i in 1:length(mean_DNAm_files)) {
  sample_name <- str_remove(basename(mean_DNAm_files[i]), "_10kb_avg_DNAm..+")
  s <- read_tsv(mean_DNAm_files[i], col_names = FALSE) %>% select(X1, X6)
  colnames(s) <- c("te_id", "mean_DNAm")
  mean_DNAm_sample_ls[[sample_name]] <- s
}

mean_DNAm_sample <- mean_DNAm_sample_ls %>% 
  bind_rows(.id = "sample")


# get antisense L1 up 10kb expression
mean_expr_files <- list.files("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/antisenseTE_expr/up10kb_antisense_expr", 
                              pattern = "up10kb_counts.txt$", full.names = TRUE)

mean_expr_sample_ls <- list()
for (j in 1:length(mean_expr_files)) {
  sample_name <- str_remove(basename(mean_expr_files[j]), "_antisenseTE_up..+")
  s_g <- read.table(mean_expr_files[j], stringsAsFactors = FALSE, 
                  header = TRUE)
  colnames(s_g) <- c("te_id", "Chr","Start", "End", "Strand",
                      "Length", "crick", "watson")
  s_g <- s_g %>% dplyr::mutate(totalCount = crick + watson) %>% 
    dplyr::select(te_id:Length, totalCount)
  mean_expr_sample_ls[[sample_name]] <- s_g
}

mean_expr_sample <- mean_expr_sample_ls %>% 
  bind_rows(.id = "sample")

## counts to lcpm
mean_expr_matrix <- mean_expr_sample %>% 
  pivot_wider(id_cols = te_id, names_from = sample, 
              values_from = totalCount) %>% 
  column_to_rownames(var = "te_id") %>% 
  as.matrix()

cpm_mean_expr_matrix <- cpm(mean_expr_matrix)

lcpm_mean_expr_matrix <- log2(cpm_mean_expr_matrix + 1)

lcpm_mean_expr_long <- lcpm_mean_expr_matrix %>% 
  as.data.frame() %>% rownames_to_column(var = "te_id") %>% 
  pivot_longer(cols = -c(te_id), 
               names_to = "sample", values_to = "lcpm") 

mean_expr_sample <- mean_expr_sample %>% 
  left_join(lcpm_mean_expr_long, by = c("sample", "te_id"))

## correlation
sample_mean <- mean_expr_sample %>% 
  dplyr::mutate(sample = str_replace_all(sample, 
                                         c("BMPC1" = "B1_rest", 
                                           "BMPC2" = "B2_rest", 
                                           "BMPC3" = "B3_rest"))) %>% 
  left_join(mean_DNAm_sample, c("sample", "te_id"))

sample_mean_cor <- 
  sample_mean %>% 
  ggplot2::ggplot(aes(x = lcpm, y = mean_DNAm))+
  geom_point()+
  geom_smooth(method = "lm")+ 
  ggpubr::stat_cor(method = "pearson", size=4, label.y = 1.1) +
  # ggrepel::geom_text_repel(aes(label = sample), size = 3)+
  labs(x= "Mean expression (lcpm)", y= "Mean DNA methylation")+
  theme_set(ggpubr::theme_pubr(base_size = 12))


# ggsave("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/antisenseTE_expr/sample_DNAm_expr_mean_Cor.pdf",
#        plot = sample_mean_cor,
#        device = "pdf", width = 12, height = 10, units = "cm")



facet_sample_mean_cor <- 
  sample_mean %>% 
  dplyr::filter(str_detect(sample, "PMM")) %>%
  dplyr::filter(sample != "PMM15") %>%
  ggplot2::ggplot(aes(x = lcpm, y = mean_DNAm))+
  geom_point()+
  geom_smooth(method = "lm")+ 
  ggpubr::stat_cor(method = "pearson", size=4, label.y = 1.1) +
  # ggrepel::geom_text_repel(aes(label = sample), size = 3)+
  facet_wrap(~sample,  ncol = 4)+ #scales = "free_x",
  labs(x= "Mean expression (lcpm)", y= "Mean DNA methylation")+
  theme_set(ggpubr::theme_pubr(base_size = 12))

ggsave("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/antisenseTE_expr/facet_PMMsample_DNAm_expr_mean_Cor.pdf",
       plot = facet_sample_mean_cor,
       device = "pdf", width = 22, height = 15, units = "cm")


PMM15_sample_mean_cor <- 
  sample_mean %>% 
  dplyr::filter(sample == "PMM15") %>% 
  ggplot2::ggplot(aes(x = lcpm, y = mean_DNAm))+
  geom_point()+
  geom_smooth(method = "lm")+ 
  ggpubr::stat_cor(method = "pearson", size=4, label.y = 1.1) +
  # ggrepel::geom_text_repel(aes(label = sample), size = 3)+
  facet_wrap(~sample, scales = "free_x")+
  labs(x= "Mean expression (lcpm)", y= "Mean DNA methylation")+
  theme_set(ggpubr::theme_pubr(base_size = 12))


ggsave("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/antisenseTE_expr/facet_samplePMM15_DNAm_expr_mean_Cor.pdf",
       plot = PMM15_sample_mean_cor,
       device = "pdf", width = 10, height = 9, units = "cm")


# B rest DNAm and BMPC expression

facet_sample_mean_cor <- 
  sample_mean %>% 
  dplyr::filter(str_detect(sample, "PMM")) %>%
  dplyr::filter(sample != "PMM15") %>%
  ggplot2::ggplot(aes(x = lcpm, y = mean_DNAm))+
  geom_point()+
  geom_smooth(method = "lm")+ 
  ggpubr::stat_cor(method = "pearson", size=4, label.y = 1.1) +
  # ggrepel::geom_text_repel(aes(label = sample), size = 3)+
  facet_wrap(~sample,  ncol = 4)+ #scales = "free_x",
  labs(x= "Mean expression (lcpm)", y= "Mean DNA methylation")+
  theme_set(ggpubr::theme_pubr(base_size = 12))

ggsave("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/antisenseTE_expr/facet_PMMsample_DNAm_expr_mean_Cor.pdf",
       plot = facet_sample_mean_cor,
       device = "pdf", width = 22, height = 15, units = "cm")


# inactive L1 plots ############################################################

# get antisense L1 up 10kb DNAm
mean_DNAm_files <- list.files("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_bowtie1_filter/methylKit_output/region_mean", 
                              pattern = "DNAm_inactiveL1.tsv", 
                              full.names = TRUE)

mean_DNAm_sample_ls <- list()
for (i in 1:length(mean_DNAm_files)) {
  sample_name <- str_remove(basename(mean_DNAm_files[i]), "_10kb_avg_DNAm..+")
  s <- read_tsv(mean_DNAm_files[i], col_names = FALSE) %>% select(X1, X6)
  colnames(s) <- c("te_id", "mean_DNAm")
  mean_DNAm_sample_ls[[sample_name]] <- s
}

mean_DNAm_sample <- mean_DNAm_sample_ls %>% 
  bind_rows(.id = "sample")


# get antisense L1 up 10kb expression
mean_expr_files <- list.files("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/antisenseTE_expr/up10kb_inactiveL1_expr", 
                              pattern = "up10kb_counts.txt$", full.names = TRUE)

mean_expr_sample_ls <- list()
for (j in 1:length(mean_expr_files)) {
  sample_name <- str_remove(basename(mean_expr_files[j]), "_antisenseTE_up..+")
  s_g <- read.table(mean_expr_files[j], stringsAsFactors = FALSE, 
                    header = TRUE)
  colnames(s_g) <- c("te_id", "Chr","Start", "End", "Strand",
                     "Length", "crick", "watson")
  s_g <- s_g %>% dplyr::mutate(totalCount = crick + watson) %>% 
    dplyr::select(te_id:Length, totalCount)
  mean_expr_sample_ls[[sample_name]] <- s_g
}

mean_expr_sample <- mean_expr_sample_ls %>% 
  bind_rows(.id = "sample")

## counts to lcpm
mean_expr_matrix <- mean_expr_sample %>% 
  pivot_wider(id_cols = te_id, names_from = sample, 
              values_from = totalCount) %>% 
  column_to_rownames(var = "te_id") %>% 
  as.matrix()

cpm_mean_expr_matrix <- cpm(mean_expr_matrix)

lcpm_mean_expr_matrix <- log2(cpm_mean_expr_matrix + 1)

lcpm_mean_expr_long <- lcpm_mean_expr_matrix %>% 
  as.data.frame() %>% rownames_to_column(var = "te_id") %>% 
  pivot_longer(cols = -c(te_id), 
               names_to = "sample", values_to = "lcpm") 

mean_expr_sample <- mean_expr_sample %>% 
  left_join(lcpm_mean_expr_long, by = c("sample", "te_id"))

## correlation
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene_gr <- genes(txdb, single.strand.genes.only=FALSE) 
gene_gr <- reduce(gene_gr, ignore.strand = TRUE)

inactiveL1 <- read_bed("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/activationScore_R2/filtered_L1/inactive_TE_up10kb.bed")

inactiveL1_fil <- inactiveL1[!overlapsAny(inactiveL1, gene_gr, ignore.strand = TRUE)]


sample_mean <- mean_expr_sample %>% 
  left_join(mean_DNAm_sample, c("sample", "te_id")) %>% 
  na.omit() %>%
  dplyr::filter(te_id %in% inactiveL1_fil$name) %>% 
  dplyr::filter(lcpm <1)


sample_mean_cor <- 
  sample_mean %>% 
  ggplot2::ggplot(aes(x = lcpm, y = mean_DNAm))+
  geom_point()+
  geom_smooth(method = "lm")+ 
  ggpubr::stat_cor(method = "pearson", size=4, label.y = 1.1) +
  # ggrepel::geom_text_repel(aes(label = sample), size = 3)+
  labs(x= "Mean expression (lcpm)", y= "Mean DNA methylation")+
  theme_set(ggpubr::theme_pubr(base_size = 12))


# ggsave("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/antisenseTE_expr/sample_DNAm_expr_mean_Cor.pdf",
#        plot = sample_mean_cor,
#        device = "pdf", width = 12, height = 10, units = "cm")



facet_sample_mean_cor <- 
  sample_mean %>% 
  dplyr::filter(str_detect(sample, "PMM")) %>% 
  dplyr::filter(sample != "PMM15") %>% 
  ggplot2::ggplot(aes(x = lcpm, y = mean_DNAm))+
  geom_point()+
  geom_smooth(method = "lm")+ 
  ggpubr::stat_cor(method = "pearson", size=4, label.y = 1.1) +
  # ggrepel::geom_text_repel(aes(label = sample), size = 3)+
  facet_wrap(~sample,  ncol = 4)+ #scales = "free_x",
  labs(x= "Mean expression (lcpm)", y= "Mean DNA methylation")+
  theme_set(ggpubr::theme_pubr(base_size = 12))

ggsave("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/antisenseTE_expr/facet_PMMsample_DNAm_expr_mean_Cor.pdf",
       plot = facet_sample_mean_cor,
       device = "pdf", width = 22, height = 15, units = "cm")


PMM15_sample_mean_cor <- 
  sample_mean %>% 
  dplyr::filter(sample == "PMM15") %>% 
  ggplot2::ggplot(aes(x = lcpm, y = mean_DNAm))+
  geom_point()+
  geom_smooth(method = "lm")+ 
  ggpubr::stat_cor(method = "pearson", size=4, label.y = 1.1) +
  # ggrepel::geom_text_repel(aes(label = sample), size = 3)+
  facet_wrap(~sample, scales = "free_x")+
  labs(x= "Mean expression (lcpm)", y= "Mean DNA methylation")+
  theme_set(ggpubr::theme_pubr(base_size = 12))


ggsave("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/antisenseTE_expr/facet_samplePMM15_DNAm_expr_mean_Cor.pdf",
       plot = PMM15_sample_mean_cor,
       device = "pdf", width = 10, height = 9, units = "cm")


