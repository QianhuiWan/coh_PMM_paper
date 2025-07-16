
library(tidyverse)
library(data.table)
library(magrittr)
library(plyranges)
library(ggplot2)
library(ggpubr)

input_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_methylGrapher_gemini/ins_supporting_reads"

################################################################################
# step 1: load INS supporting reads for PMMs and get CPM 
################################################################################

# INS support reads for all samples:
ins_files <- list.files(input_dir, pattern = "_ins_reads_nodes.txt$", recursive = TRUE, full.names = TRUE)

# get sample name
samples <- sub("_ins_reads_nodes.txt", "", basename(ins_files))

# calculate INS supporting reads number

ins_read_counts <- rbindlist(lapply(seq_along(ins_files), function(i) {
  f <- ins_files[i]
  s <- samples[i]
  df <- fread(f)
  
  # Count supporting reads per INS_id
  ins_support <- df[, .(ins_reads_num = uniqueN(read_id)), by = ins_id]
  
  # Add sample name
  ins_support[, sample := s]
  
  # Reorder columns for clarity
  ins_support[, .(sample, ins_id, ins_reads_num)]
}))


# load total reads
total_reads <- rbindlist(lapply(samples, function(s) {
  f <- file.path(input_dir, s, paste0(s, "_total_reads.txt"))
  if (file.exists(f)) {
    n <- as.integer(readLines(f, warn = FALSE))
    data.table(sample = s, total_reads = n)
  } else {
    data.table(sample = s, total_reads = NA_integer_)
  }
}))

# merge
merged <- merge(ins_read_counts, total_reads, by = "sample")

# calculate CPM
merged[, INS_reads_CPM := (ins_reads_num / total_reads) * 1e6]

# order by CPM
merged <- merged[order(-INS_reads_CPM)]

# save output
fwrite(merged, paste0(input_dir, "/INS_supporting_reads_CPM.tsv"), sep = "\t")

# print first 6 output
print(head(merged))

################################################################################
# step 2: filter INS and get PMM INS
################################################################################

merged <- fread("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_methylGrapher_gemini/ins_supporting_reads/INS_supporting_reads_CPM.tsv")

control_ins_ids <- unique(merged[grepl("rest", sample), ins_id]) #309911

# ins annotation #310716
ins_anno <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/methylGrapher_ref/pangenomeInsertion/pan_ins_CHM13graph_linearCoords_REF_ALT_Length50_nodes_ids.tsv", col_names = FALSE) 
ins_anno <- janitor::row_to_names(ins_anno, row_number = 1)
# ins_anno_6000 <- ins_anno[ins_anno$ins_len <7000 & ins_anno$ins_len > 5900,]
ins_anno_6000 <- ins_anno[ins_anno$ins_len <7000 & ins_anno$ins_len > 1,] #299968
# ins_anno_6000 <- ins_anno[ins_anno$ins_len <1000,]

aligned_ins <- read.table("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/methylGrapher_ref/pangenomeInsertion/CHM13_INS_plus_flank1000bp_to_graph_gaf_uniqINSids.txt")
colnames(aligned_ins) <- "ins_id" #308641

PMM_ins <- merged %>% 
  # dplyr::filter(INS_reads_CPM > quantile(INS_reads_CPM, probs = 0.25)) %>%
  dplyr::filter(!ins_id %in% control_ins_ids) %>% 
  dplyr::filter(ins_id %in% aligned_ins$ins_id) %>%
  dplyr::filter(INS_reads_CPM > quantile(INS_reads_CPM, probs = 0.1)) %>%
  dplyr::filter(ins_id %in% ins_anno_6000$ins_id)

fwrite(PMM_ins, paste0(input_dir, "/INS_supporting_reads_inPMM.tsv"), sep = "\t")

ins_counts <- PMM_ins %>% 
  group_by(sample) %>%
  summarise(active_ins = n_distinct(ins_id))

arrange(ins_counts, -active_ins)
fwrite(arrange(ins_counts, -active_ins), 
       paste0(input_dir, "/INS_num_inPMM.tsv"), sep = "\t")

################################################################################
# step 3: plot PMM INS number against PMM L1 number
################################################################################

pmm_L1num <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/activationScore_R2/filtered_L1/num_active_all_L1s_L1PA2_L1PA3_pmm.tsv")

pmm_L1num <- pmm_L1num %>% 
  dplyr::select(sample = sampleID, L1) %>% 
  unique() %>% 
  left_join(ins_counts, by = c("sample")) %>% 
  na.omit()
  # mutate(active_ins = replace_na(active_ins, 0))

# fit model
fit <- lm(active_ins ~ L1, data = pmm_L1num)
summary_fit <- summary(fit)

# get R² and p value
r2 <- summary_fit$r.squared
pval <- summary_fit$coefficients[2,4]

# plot
p <- ggplot(pmm_L1num, aes(x = L1, y = active_ins)) +
  geom_point(size = 3, color = "#0072B2") +  
  geom_smooth(method = "lm", se = FALSE, color = "#D55E00", linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = sample), size = 3)+
  annotate("text",
           x = 100, y = 620,
           label = paste0("R² = ", signif(r2, 3), "\np = ", signif(pval, 3)),
           hjust = 1, vjust = 1, size = 4.5, color = "black") +
  labs(
    # title = "Active L1 vs. Read-supported Insertions",
    x = "Number of active L1s",
    y = "Number of read-supported insertions") +
  theme_pubr(base_size = 9.6)

ggsave(paste0(input_dir,"/L1_ins_cor_rmControls.pdf"),
       plot = p,
       device = "pdf", width = 12, height = 10, units = "cm")


################################################################################
# step 4: get MEI and overlap with PMM INS
################################################################################

## get PMM_ins_GR
all_PMM_restB_ins <- merged 
all_PMM_restB_ins_GR <- ins_anno %>% 
  dplyr::filter(ins_id %in% all_PMM_restB_ins$ins_id) %>%
  dplyr::mutate(start = as.integer(start), end = start) %>% 
  dplyr::select(seqnames, start, end, ins_id, ins_len) %>% 
  as_granges()
  
all_MEI <- read_csv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/t2t_CHM13/MEI_HGSVC3/MEI_Callset_T2T-CHM13.ALL.20241211.csv.gz")

all_MEI_GR <- all_MEI %>% 
  dplyr::mutate(seqnames = CHROM, start = as.integer(POS), end = start) %>% 
  dplyr::select(seqnames, start, end, ID, TE_Designation) %>% 
  as_granges()

all_PMM_restB_ins_MEI_GR <- find_overlaps(all_PMM_restB_ins_GR, all_MEI_GR)
all_PMM_restB_ins_MEI_df <- all_PMM_restB_ins_MEI_GR %>% as.data.frame() %>% 
  dplyr::filter(ID != "chr9-141236608-INS-1578")

all_PMM_restB_ins_addMEIinfo <- all_PMM_restB_ins %>% 
  dplyr::filter(ins_id %in% all_PMM_restB_ins_MEI_df$ins_id) %>% 
  dplyr::left_join(all_PMM_restB_ins_MEI_df, by = "ins_id")


################################################################################
# step 5: plot PMM MEI number against PMM L1 number
################################################################################


ins_counts <- all_PMM_restB_ins_addMEIinfo %>% 
  mutate(group = str_replace_all(sample, "[:digit:]", "")) %>%
  dplyr::filter(group != "B_EBV") %>% 
  dplyr::filter(INS_reads_CPM > quantile(INS_reads_CPM, prob = 0.5)) %>%
  # dplyr::filter(INS_reads_CPM > 0.7) %>% 
  dplyr::filter(ins_reads_num > 1) %>%
  group_by(sample, TE_Designation) %>%
  summarise(active_ins = n_distinct(ins_id))

arrange(ins_counts, -active_ins)

ins_counts %>%
  mutate(group = str_replace_all(sample, "[:digit:]", "")) %>%
  dplyr::filter(group != "B_EBV") %>% 
  ggplot(aes(x = sample, y = active_ins)) +
  geom_point(aes(fill = group), alpha = 0.3) +
  facet_wrap(~TE_Designation)+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ins_counts %>%
  mutate(group = str_replace_all(sample, "[:digit:]", "")) %>%
  dplyr::filter(group != "B_EBV") %>% 
  ggplot(aes(x = group, y = active_ins)) +
  geom_point(aes(fill = group), alpha = 0.3) +
  facet_wrap(~TE_Designation)+
  stat_compare_means(method = "wilcox.test", label = "p.signif")+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


pmm_L1num <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/activationScore_R2/filtered_L1/num_active_all_L1s_L1PA2_L1PA3_pmm.tsv")

pmm_L1num <- pmm_L1num %>% 
  dplyr::select(sample = sampleID, L1) %>% 
  unique() %>% 
  left_join(ins_counts, by = c("sample")) %>% 
  na.omit()
# mutate(active_ins = replace_na(active_ins, 0))

# fit model
fit <- lm(active_ins ~ L1, data = pmm_L1num[pmm_L1num$TE_Designation=="LINE/L1",])
fit <- lm(active_ins ~ L1, data = pmm_L1num[pmm_L1num$TE_Designation=="Retroposon/SVA",])
summary_fit <- summary(fit)

# get R² and p value
r2 <- summary_fit$r.squared
pval <- summary_fit$coefficients[2,4]

# plot

p_facet <- ggplot(pmm_L1num, aes(x = L1, y = active_ins)) +
  geom_point(size = 3, color = "#0072B2") +  
  geom_smooth(method = "lm", se = FALSE, color = "#D55E00", linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = sample), size = 3) +
  # scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  stat_cor(method = "pearson", label.y = 1500
           # label.x.npc = "left",
           # label.y.npc = 0.99   # numeric, between 0 and 1 (1 = very top)
           )+  
  facet_wrap(~TE_Designation)+
  labs(
    x = "Number of active L1s",
    y = "Number of insertions") +
  theme_pubr(base_size = 9.6)

ggsave(paste0(input_dir,"/L1_allMEIins_cor_facet.pdf"),
       plot = p_facet,
       device = "pdf", width = 22, height = 12, units = "cm")

p_sva <- pmm_L1num %>% 
  dplyr::filter(TE_Designation == "Retroposon/SVA") %>% 
  ggplot(aes(x = L1, y = active_ins)) +
  geom_point(size = 3, color = "#0072B2") +  
  geom_smooth(method = "lm", se = FALSE, color = "#D55E00", linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = sample), size = 3)+
  # annotate("text",
  #          x = 100, y = 260,
  #          label = paste0("R² = ", signif(r2, 3), "\np = ", signif(pval, 2)), 
  #          hjust = 1, vjust = 1, size = 4.5, color = "black") +
  stat_cor(method = "pearson", label.y = 260)+  
  labs(
    # title = "Active L1 vs. Read-supported Insertions",
    x = "Number of active L1s",
    y = "Number of SVA insertions") +
  theme_pubr(base_size = 9.6)

ggsave(paste0(input_dir,"/L1_svaMEIins_cor.pdf"),
       plot = p_sva,
       device = "pdf", width = 12, height = 10, units = "cm")



################################################################################
# step 6: other plot exploration
################################################################################

all_PMM_restB_ins_addMEIinfo %>%
  mutate(group = str_replace_all(sample, "[:digit:]", "")) %>%
  dplyr::filter(group != "B_EBV") %>% 
  dplyr::filter(INS_reads_CPM > 1) %>% 
  dplyr::filter(ins_reads_num > 10) %>% 
  ggplot(aes(x = group, y = INS_reads_CPM)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_point(position = position_jitter(width = 0.15), alpha = 0.3, size = 1) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "red", fatten = 2) +
  facet_wrap(~TE_Designation) +
  stat_compare_means(method = "wilcox.test", label = "p.signif")

all_PMM_restB_ins_addMEIinfo %>%
  mutate(group = str_replace_all(sample, "[:digit:]", "")) %>%
  dplyr::filter(group != "B_EBV") %>% 
  dplyr::filter(INS_reads_CPM > 0.5) %>% 
  dplyr::filter(ins_reads_num > 5) %>% 
  ggplot(aes(x = sample, y = INS_reads_CPM)) +
  geom_violin(aes(fill = group), alpha = 0.3, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1)+
  facet_wrap(~TE_Designation)+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))





