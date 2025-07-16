
outDir <- '/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1'

library(tidyverse)
library(data.table)
library(plyranges)
library(rtracklayer)
PMM_metadata <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/PMM_sampleSheet.tsv")
# PMM_metadata$group <- str_replace_all(PMM_metadata$group, "BMPC", "control")
PMM_ss_L1PAnum <- PMM_metadata %>% 
  dplyr::select(sampleName, L1PA2, L1PA3) %>% 
  dplyr::filter(str_detect(sampleName, "PMM"))

# PMD info
pmd_PMMs <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/pmd_allPMM_samples/allPMM_samples_pmd.bed")
pmd_PMMs_GR <- import(pmd_PMMs) %>% 
  plyranges::filter(seqnames %in%  c(paste0("chr", c(1:22, "X", "Y")))) %>% 
  plyranges::mutate(pmd_id = paste0("pmd_", 1:length(.)))

# S1: get enhancer for H1 and H9
# fantom5_enh <- fread("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/fantom5_enhancer/F5.hg38.enhancers.expression.usage.matrix.gz", 
#                      sep = "\t", header = TRUE, nThread = 1)
fantom5_enh <- read.table("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/fantom5_enhancer/F5.hg38.enhancers.expression.usage.matrix.gz",
                          sep = "\t", header = TRUE, row.names = 1)

fantom5_enh_ESC <- 
  fantom5_enh %>% 
  dplyr::select(CNhs14067, CNhs14068, CNhs13964, 
                CNhs11917, CNhs12824, CNhs12837) %>%  # H1 and H9
  rownames_to_column(var = "enhancer_coord")

## filter enhancer in use (enhancer usage ==1) ==3
fantom5_enh_ESC_H1 <- fantom5_enh_ESC %>% 
  dplyr::select(enhancer_coord, CNhs14067, CNhs14068, CNhs13964) %>% 
  rowwise() %>%
  filter(sum(c_across(-enhancer_coord) == 1) == 3)%>% 
  as.data.frame()

write_tsv(fantom5_enh_ESC_H1, 
          file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/fantom5_enhancer/fantom5_enh_inUse_ESC_H1.tsv")

fantom5_enh_ESC_H9 <- fantom5_enh_ESC %>% 
  dplyr::select(enhancer_coord, CNhs11917, CNhs12824, CNhs12837) %>% 
  rowwise() %>%
  filter(sum(c_across(-enhancer_coord) == 1) == 3) %>% 
  as.data.frame()

write_tsv(fantom5_enh_ESC_H9, 
          file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/fantom5_enhancer/fantom5_enh_inUse_ESC_H9.tsv")


fantom5_enh_ESC_H1 <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/fantom5_enhancer/fantom5_enh_inUse_ESC_H1.tsv")
fantom5_enh_ESC_H9 <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/fantom5_enhancer/fantom5_enh_inUse_ESC_H9.tsv")

## get ehn GRs
### H1
fantom5_enh_ESC_H1_GR <- fantom5_enh_ESC_H1 %>% 
  dplyr::mutate(seqnames = str_extract_all(enhancer_coord, "^[^:]+"),
                start = as.integer(str_extract_all(enhancer_coord, "(?<=:)[0-9]+")),
                end =  as.integer(str_extract_all(enhancer_coord, "(?<=-)[0-9]+"))
                ) %>% 
  as_granges()

### H9
fantom5_enh_ESC_H9_GR <- fantom5_enh_ESC_H9 %>% 
  dplyr::mutate(seqnames = str_extract_all(enhancer_coord, "^[^:]+"),
                start = as.integer(str_extract_all(enhancer_coord, "(?<=:)[0-9]+")),
                end =  as.integer(str_extract_all(enhancer_coord, "(?<=-)[0-9]+"))
  ) %>% 
  as_granges()




# S2: overlap with TE and L1 regions respectively
## Load RepeatMasker BED file
repeat_bed_path <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_repeatMasker_hg38_addTEid.tsv"
repeat_bed_df <- read_tsv(repeat_bed_path, 
                          col_types = cols(
                            genoName = col_character(),
                            genoStart = col_integer(),
                            genoEnd = col_integer(),
                            strand = col_character(),
                            te_id = col_character(),
                            repName = col_character(),
                            repClass = col_character(),
                            repFamily = col_character(),
                            age = col_double()
                          ))

## Convert to GRanges
repeat_ranges <- 
  repeat_bed_df %>%
  dplyr::mutate(seqnames = genoName, 
                start = genoStart + 1, end = genoEnd) %>%
  dplyr::select(-c(genoName:genoEnd)) %>% 
  as_granges()

## get L1s
repeat_ranges_L1 <- 
  repeat_ranges %>% 
  dplyr::filter(repFamily == "L1")
## get L1HS
repeat_ranges_L1HS <- 
  repeat_ranges %>% 
  dplyr::filter(repName == "L1HS")
## get L1PA2s
repeat_ranges_L1PA2 <- 
  repeat_ranges %>% 
  dplyr::filter(repName == "L1PA2")
## get L1PA3s
repeat_ranges_L1PA3 <- 
  repeat_ranges %>% 
  dplyr::filter(repName == "L1PA3")

## find overlaps
fantom5_enh_ESC_H1_TE_GR <- find_overlaps(fantom5_enh_ESC_H1_GR, repeat_ranges) #276
fantom5_enh_ESC_H1_L1_GR <- find_overlaps(fantom5_enh_ESC_H1_GR, repeat_ranges_L1) #12
find_overlaps(fantom5_enh_ESC_H1_TE_GR, pmd_PMMs_GR) #26
find_overlaps(fantom5_enh_ESC_H1_L1_GR, pmd_PMMs_GR) #0

fantom5_enh_ESC_H9_TE_GR <- find_overlaps(fantom5_enh_ESC_H9_GR, repeat_ranges) #475
fantom5_enh_ESC_H9_L1_GR <- find_overlaps(fantom5_enh_ESC_H9_GR, repeat_ranges_L1) #25
find_overlaps(fantom5_enh_ESC_H9_TE_GR, pmd_PMMs_GR) #76
find_overlaps(fantom5_enh_ESC_H9_L1_GR, pmd_PMMs_GR) #4

# S4: get DNA methylation mean at overalpped enhancers

meth_snps_filGR <- readRDS("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1/meth_snps_cov_fil_matGR_withXY.rds")
nms <- colnames(values(meth_snps_filGR))
cols <- lapply(nms, function(.) { col <- ensym(.); quo(mean(!!col)) })
# the column names for the summarize op
names(cols) <- paste0("mean_", nms)


fantom5_TEenh_ESC_H9_DNAm_mat <- 
  fantom5_enh_ESC_H9_TE_GR %>% 
  find_overlaps(meth_snps_filGR) %>% 
  plyranges::group_by(enhancer_coord) %>% 
  plyranges::summarise(!!! cols) %>% 
  as.data.frame() %>% 
  `colnames<-`(str_replace_all(colnames(.), "mean_", "")) %>% 
  as.data.frame() %>% 
  dplyr::select(-starts_with("B")) %>%
  # pivot_longer(cols = c(starts_with("PMM")),
  #              names_to = "sample",
  #              values_to = "enhancer_DNAm") %>%
  column_to_rownames(var = "enhancer_coord") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "sampleName") %>% 
  dplyr::filter(sampleName %in%  PMM_ss_L1PAnum$sampleName) %>% 
  left_join(PMM_ss_L1PAnum, by = "sampleName") %>% 
  column_to_rownames(var = "sampleName") %>% 
  as.matrix()


# Compute correlation & p-values
cor_results <- rcorr(fantom5_enh_ESC_H9_DNAm_mat)

# Extract correlation matrix and p-values
cor_results$r[, "L1PA2"]  # Correlation coefficients

keep <- cor_results$P[, "L1PA2"] < 0.05

enhancer_coord_sigCorL1PA2 <- cor_results$P[, "L1PA2"][keep]
enhancer_coord_sigCorL1PA2_r <- cor_results$r[, "L1PA2"][keep]

enhancer_coord_sigCorL1PA2_df <- data.frame(
  cor_p = enhancer_coord_sigCorL1PA2[1:(length(enhancer_coord_sigCorL1PA2)-2)],
  cor = enhancer_coord_sigCorL1PA2_r[1:(length(enhancer_coord_sigCorL1PA2_r)-2)]
) %>% 
  rownames_to_column(var = "enhancer_coord") %>% 
  dplyr::mutate(group = ifelse(cor >0.7, "pos", "neg"))



# cor between these enhancer DNAm and L1PA2 numbers
# sig: APP_CGI39, SECTM1_CGI33
p1 <- 
  fantom5_TEenh_ESC_H9_DNAm_mat %>% as.data.frame() %>% 
  rownames_to_column(var = "sampleName") %>% 
  pivot_longer(cols = c(starts_with("chr")),
               names_to = "enhancer_coord",
               values_to = "enhancer_DNAm") %>%
  dplyr::filter(enhancer_coord %in% 
                  enhancer_coord_sigCorL1PA2_df$enhancer_coord) %>% 
  left_join(enhancer_coord_sigCorL1PA2_df, by = "enhancer_coord") %>% 
  ggplot(aes(x = L1PA2, y = enhancer_DNAm, color = group)) + 
  geom_point() + geom_smooth(method = "lm")+ 
  # facet_wrap(~enhancer_coord, scales = "free_x")+
  ggpubr::stat_cor(method = "pearson", size=3, label.y = c(1,1.2)) +
  # ggrepel::geom_text_repel(aes(label = sample), size = 3)+
  labs(x = "number of L1PA2", y="DNA methylation")+
  theme_set(ggpubr::theme_pubr(base_size = 12))


ggsave(paste0(outDir, "/fantom5_enh_ESC_H9_DNAm_L1PA2_cor_test.pdf"),
       plot = p1, 
       device = "pdf", 
       width = 18, height = 10,
       # width = 22, height = 35, 
       units = "cm")  


fantom5_TEenh_ESC_H1_DNAm_mat <- 
  fantom5_enh_ESC_H1_TE_GR %>% 
  find_overlaps(meth_snps_filGR) %>% 
  plyranges::group_by(enhancer_coord) %>% 
  plyranges::summarise(!!! cols) %>% 
  as.data.frame() %>% 
  `colnames<-`(str_replace_all(colnames(.), "mean_", "")) %>% 
  as.data.frame() %>% 
  dplyr::select(-starts_with("B")) %>%
  # pivot_longer(cols = c(starts_with("PMM")),
  #              names_to = "sample",
  #              values_to = "enhancer_DNAm") %>%
  column_to_rownames(var = "enhancer_coord") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "sampleName") %>% 
  dplyr::filter(sampleName %in%  PMM_ss_L1PAnum$sampleName) %>% 
  left_join(PMM_ss_L1PAnum, by = "sampleName") %>% 
  column_to_rownames(var = "sampleName") %>% 
  as.matrix()


# Compute correlation & p-values
cor_results <- rcorr(fantom5_TEenh_ESC_H1_DNAm_mat)

# Extract correlation matrix and p-values
cor_results$r[, "L1PA2"]  # Correlation coefficients

keep <- cor_results$P[, "L1PA2"] < 0.05

enhancer_coord_sigCorL1PA2 <- cor_results$P[, "L1PA2"][keep]
enhancer_coord_sigCorL1PA2_r <- cor_results$r[, "L1PA2"][keep]

enhancer_coord_sigCorL1PA2_df <- data.frame(
  cor_p = enhancer_coord_sigCorL1PA2[1:(length(enhancer_coord_sigCorL1PA2)-2)],
  cor = enhancer_coord_sigCorL1PA2_r[1:(length(enhancer_coord_sigCorL1PA2_r)-2)]
) %>% 
  rownames_to_column(var = "enhancer_coord") %>% 
  dplyr::mutate(group = ifelse(cor >0.7, "pos", "neg"))


# cor between these enhancer DNAm and L1PA2 numbers
# sig: APP_CGI39, SECTM1_CGI33
p2 <- 
  fantom5_TEenh_ESC_H1_DNAm_mat %>% as.data.frame() %>% 
  rownames_to_column(var = "sampleName") %>% 
  pivot_longer(cols = c(starts_with("chr")),
               names_to = "enhancer_coord",
               values_to = "enhancer_DNAm") %>%
  dplyr::filter(enhancer_coord %in% 
                  enhancer_coord_sigCorL1PA2_df$enhancer_coord) %>% 
  left_join(enhancer_coord_sigCorL1PA2_df, by = "enhancer_coord") %>% 
  ggplot(aes(x = L1PA2, y = enhancer_DNAm, color = group)) + 
  geom_point() + geom_smooth(method = "lm")+ 
  # facet_wrap(~enhancer_coord, scales = "free_x")+
  ggpubr::stat_cor(method = "pearson", size=3, label.y = c(1.22,1.5)) +
  # ggrepel::geom_text_repel(aes(label = sample), size = 3)+
  labs(x = "number of L1PA2", y="DNA methylation")+
  theme_set(ggpubr::theme_pubr(base_size = 12))


ggsave(paste0(outDir, "/fantom5_enh_ESC_H1_DNAm_L1PA2_cor_test.pdf"),
       plot = p2, 
       device = "pdf", 
       width = 18, height = 10,
       # width = 22, height = 35, 
       units = "cm")  


# any related genes??
# H9: ZNF287
# H1: te_2386985, both enhancer and KLF6 binding site on chr7

# plot te_2386985
p3 <- 
  fantom5_TEenh_ESC_H1_DNAm_mat %>% as.data.frame() %>% 
  rownames_to_column(var = "sampleName") %>% 
  pivot_longer(cols = c(starts_with("chr")),
               names_to = "enhancer_coord",
               values_to = "enhancer_DNAm") %>%
  dplyr::filter(enhancer_coord %in% "chr7:138243527-138244000") %>% 
  left_join(enhancer_coord_sigCorL1PA2_df, by = "enhancer_coord") %>% 
  ggplot(aes(x = L1PA2, y = enhancer_DNAm, color = group)) + 
  geom_point() + geom_smooth(method = "lm")+ 
  facet_wrap(~enhancer_coord, scales = "free_x")+
  ggpubr::stat_cor(method = "pearson", size=3, label.y = c(1,1.2)) +
  ggrepel::geom_text_repel(aes(label = sampleName), size = 3)+
  labs(x = "number of L1PA2", y="DNA methylation")+
  theme_set(ggpubr::theme_pubr(base_size = 12))


ggsave(paste0(outDir, "/fantom5_enh_ESC_H1_DNAm_L1PA2_cor_L1ME5.pdf"),
       plot = p3, 
       device = "pdf", 
       width = 18, height = 10,
       # width = 22, height = 35, 
       units = "cm")  



