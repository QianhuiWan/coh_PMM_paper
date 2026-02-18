# sample hyper methyl regions in PMDs
# + check CGI and CGI +1000bp CGI-1000bp regions (overlap with hyper methyl regions)
# L1 not filtered by 5900bp

library(tidyverse)
library(data.table)
library(magrittr)
library(rtracklayer)
library(plyranges)
library(ggplot2)
library(ggpubr)

# load PMDs identified using all PMM sample ####################################
pmd_PMMs <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/pmd_allPMM_samples/allPMM_samples_pmd.bed")
pmd_PMMs_GR <- import(pmd_PMMs) %>%
  plyranges::filter(seqnames %in% c(paste0("chr", c(1:22, "X", "Y")))) %>%
  plyranges::mutate(pmd_id = paste0("pmd_", 1:length(.)))

# load DNAm mat for all samples (PMM + B cells) ################################
meth_mat_Dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1"
meth_snps_cov_fil_mat <- readRDS(paste0(meth_mat_Dir, "/meth_snps_cov_fil_matGR_withXY.rds"))

meth_snps_fil_mat_addBMPC <- readRDS("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1/meth_snps_fil_matGR_withXY_addBMPC.rds")
meth_snps_cov_fil_mat <- meth_snps_cov_fil_mat %>%
  plyranges::join_overlap_left(meth_snps_fil_mat_addBMPC)

# load hyper methyl regions, for each samples in for loop
## first get sample id and group ###############################################
# file.list <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/*/05.HyperMR/*_CpG.hypermr")
## remove EBV samples
# file.list <- file.list[!str_detect(file.list, "EBV")]
outdir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/pmd_res"

## repeat and CpG annotations ##################################################
rmsk_hg38 <- fread("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/repeatMasker_hg38")
rmsk_hg38_te_ids <- fread("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/hg38_rmsk_teCoord_ids.txt",
  col.names = c("genoName", "genoStart", "genoEnd", "strand", "te_id")
)
rmsk_hg38_te_ids <- rmsk_hg38_te_ids %>% mutate(genoEnd = genoEnd - 1)
### get repeat GRange ##########################################################
repeat_GR <- rmsk_hg38 %>%
  # dplyr::mutate(div = milliDiv/1000) %>% ## Convert percent divergence of RepeatMasker output to divergence
  # dplyr::mutate(div_jc = -(0.75*log(1 - (4*div)/3))) %>% ## Apply Juke's Cantor correction to divergence distance
  # dplyr::mutate(age=(div_jc*1000) / 2.2) %>%
  # dplyr::filter(age<=40) %>%
  dplyr::left_join(rmsk_hg38_te_ids,
    by = c(
      "genoName" = "genoName", "genoStart" = "genoStart",
      "genoEnd" = "genoEnd", "strand" = "strand"
    )
  ) %>%
  dplyr::select(
    seqnames = genoName, start = genoStart, end = genoEnd,
    strand, te_id, repName, repClass, repFamily
  ) %>%
  dplyr::mutate(repStrand = strand) %>%
  dplyr::mutate(start = start + 1) %>%
  plyranges::filter(seqnames %in% c(paste0("chr", c(1:22, "X", "Y")))) %>%
  as_granges()

### get CpGisland GRange #######################################################
cpgIsland_GR <- rtracklayer::import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/cpgIsland/UCSC_CpGislands_unmask.bed") %>%
  plyranges::filter(seqnames %in% c(paste0("chr", c(1:22, "X", "Y")))) %>%
  plyranges::mutate(cpgIsland_id = paste("cpg", 1:length(.), sep = "_")) %>%
  mutate(CGI_start = start, CGI_end = end)

PMM_sampleSheet <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/gene_counts/PMM_sampleSheet_R3_addBMPC.tsv")

# for loop to get hyper regions in PMD #########################################
## filter hypermr for PMD and L1 repeat regions for each sample one by one
# PMD_hypermr_dfl <- list()
# L1_5P_hypermr_dfl <- list()
# L1_CpGisland_info_dfl <- list()
# for (i in 1:ncol(mcols(meth_snps_cov_fil_mat))) {
#   sampleID <- colnames(mcols(meth_snps_cov_fil_mat))[i]
#
#   # 3. get L1 CpGislands, to consider L1 CpG islands
#   ## all L1 CpG islands
#   repeat_L1_GR <- repeat_GR[repeat_GR$repFamily=="L1",]
#   # repeat_L1_GR_full <- repeat_L1_GR %>% plyranges::filter(width(.) >= 5900)
#   cpgIsland_L1_GR <- find_overlaps(cpgIsland_GR, repeat_L1_GR)
#   ## L1 CpG islands in PMM PMDs
#   cpgIsland_L1_PMD_GR <- find_overlaps(cpgIsland_L1_GR, pmd_PMMs_GR)
#
#   # 5. calculate mean DNAm for L1 CpG islands
#   cpgIsland_L1_DNAm <-
#     find_overlaps(cpgIsland_L1_GR, meth_snps_cov_fil_mat[, sampleID]) %>%
#     plyranges::group_by(cpgIsland_id, te_id) %>%
#     plyranges::summarise(cpgIsland_DNAm = mean(!!sym(sampleID))) %>%
#     as.data.frame() %>%
#     left_join(as.data.frame(mcols(repeat_L1_GR)), by = "te_id")
#
#   # 7. join mean DNAm of L1 CpG island with hypermr info.
#   ## Extract metadata as data frames
#   cpgIsland_L1_PMD_GR_metadata <-
#     as.data.frame(mcols(cpgIsland_L1_PMD_GR)) %>%
#     dplyr::select(cpgIsland_id, pmd_id) %>%
#     unique()
#   ## Perform a left_join on data
#   cpgIsland_L1_DNAm_addmeta <- cpgIsland_L1_DNAm %>%  # 2282 unique L1 CpG island
#     dplyr::left_join(cpgIsland_L1_PMD_GR_metadata,
#                      by = c("cpgIsland_id")) %>%
#     dplyr::mutate(overlap_PMD = ifelse(is.na(pmd_id), "No", "Yes"))
#
#   # 8. write each sample res. to the list
#   L1_CpGisland_info_dfl[[sampleID]] <- cpgIsland_L1_DNAm_addmeta
# }

# save the for loop outputs ####################################################
# saveRDS(PMD_hypermr_dfl, file = paste0(outdir, "/PMD_hypermr_dfl_addBMPC.rds"))
# saveRDS(L1_5P_hypermr_dfl, file = paste0(outdir, "/L1_5P_hypermr_dfl_addBMPC.rds"))
# saveRDS(L1_CpGisland_info_dfl, file = paste0(outdir, "/L1_CpGisland_info_dfl_addBMPC.rds"))

# PMD_hypermr_dfl <- readRDS(paste0(outdir, "/PMD_hypermr_dfl_addBMPC.rds"))
# L1_5P_hypermr_dfl <- readRDS(paste0(outdir, "/L1_5P_hypermr_dfl_addBMPC.rds"))
L1_CpGisland_info_dfl <- readRDS(paste0(outdir, "/L1_CpGisland_info_dfl_addBMPC.rds"))


# plot for all L1 subfamily and the subfamilies with DNAm difference
# plot L1CpGisland DNAm against hypermr number
cpgIsland_L1_meanDNAm_df <- L1_CpGisland_info_dfl %>%
  lapply(function(x) {
    x <- dplyr::select(
      x, cpgIsland_id, cpgIsland_DNAm, te_id, repName,
      repFamily, repClass, repStrand, pmd_id, overlap_PMD
    )
    # x = dplyr::select(x, cpgIsland_id, cpgIsland_DNAm, cpgIsland_3p1000_DNAm,
    #                   te_id, repName, repFamily, repClass,repStrand,
    #                   hypermr_num, hypermr_meanDNAm, pmd_id, overlap_PMD)
    x <- unique(x)
  }) %>%
  dplyr::bind_rows(.id = "sample") %>%
  dplyr::mutate(group = str_replace_all(sample, "[:digit:]", ""))


## plot all samples

# filter for L1 CpG island have DNAm low in PMM
cpgIsland_L1_meanDNAm_delta <-
  cpgIsland_L1_meanDNAm_df %>%
  dplyr::group_by(cpgIsland_id, group) %>%
  dplyr::summarise(
    cpgIsland_DNAm_group = mean(cpgIsland_DNAm),
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols = cpgIsland_id,
    names_from = group,
    values_from = cpgIsland_DNAm_group
  ) %>%
  dplyr::mutate(delta_DNAm = PMM - B_rest)

## filter 5900bp

activeL1_sense <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/activationScore/filtered_L1/sense_L1s.tsv")
activeL1_antisense <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/activationScore/filtered_L1/antisense_L1s.tsv")


p3 <- cpgIsland_L1_meanDNAm_df %>%
  # dplyr::filter(te_id %in% activeL1_antisense$Geneid) %>%
  ggplot2::ggplot(aes(
    x = reorder(sample, cpgIsland_DNAm,
      FUN = median,
      decreasing = FALSE
    ),
    y = cpgIsland_DNAm, fill = group
  )) +
  # geom_bar(stat = "identity")+
  # geom_boxplot()+
  # geom_smooth(method = "lm")+
  facet_wrap(~overlap_PMD,
    labeller = as_labeller(c(
      "Yes" = "In PMDs",
      "No" = "Not in PMDs"
    ))
  ) +
  # stat_boxplot(geom = "errorbar",
  #              width = 0.15, linetype = 2, # Line type
  #              lwd = 0.5)+
  geom_boxplot(outlier.shape = NA) +
  theme_pubr(base_size = 9) + # Set base font size
  guides(colour = "none") +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.key.size = unit(0.4, "lines"), # make the color square smaller
    legend.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9, face = "bold", colour = "black"), # or just size = 9
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "", y = "DNAm of L1 CGIs", fill = "group")

ggsave(paste0(outdir, "/L1CGI_boxplot_R4_similarToFig1PlotBC_addBMPC.pdf"),
  plot = p3, width = 12, height = 5, units = "cm",
  limitsize = FALSE
)

p4 <- cpgIsland_L1_meanDNAm_df %>%
  dplyr::filter(te_id %in%
    c(activeL1_antisense$Geneid, activeL1_sense$Geneid)) %>%
  ggplot2::ggplot(aes(
    x = reorder(sample, cpgIsland_DNAm),
    y = cpgIsland_DNAm, colour = group
  )) +
  facet_wrap(~overlap_PMD,
    labeller = as_labeller(c(
      "Yes" = "In PMDs",
      "No" = "Not in PMDs"
    ))
  ) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.15, linetype = 2, # Line type
    lwd = 0.5
  ) +
  geom_boxplot(outlier.shape = NA) +
  theme_pubr(base_size = 9) + # Set base font size
  guides(colour = "none") +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.key.size = unit(0.4, "lines"), # make the color square smaller
    legend.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9, face = "bold", colour = "black"), # or just size = 9
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "", y = "DNAm of L1 CGIs", fill = "group")

ggsave(paste0(outdir, "/activeL1CGI_boxplot_R4_similarToFig1PlotBC_addBMPC.pdf"),
  plot = p4, width = 12, height = 5, units = "cm",
  limitsize = FALSE
)


# plot together to add p-values
## L1 CGI
p3_addP <- cpgIsland_L1_meanDNAm_df %>%
  # dplyr::filter(te_id %in% activeL1_antisense$Geneid) %>%
  # ggplot2::ggplot(aes(x = reorder(sample, cpgIsland_DNAm),
  #                     y=cpgIsland_DNAm, fill=overlap_PMD))+
  ggplot2::ggplot(aes(
    x = reorder(sample, cpgIsland_DNAm,
      FUN = function(x) median(x, na.rm = FALSE),
      decreasing = FALSE
    ),
    y = cpgIsland_DNAm, fill = overlap_PMD
  )) +
  # facet_wrap(~overlap_PMD,
  #            labeller = as_labeller(c(
  #              "Yes" = "In PMDs",
  #              "No" = "Not in PMDs"
  #            )))+
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  theme_pubr(base_size = 7) + # Set base font size
  theme(
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.4, "lines"), # make the color square smaller
    legend.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 7, face = "bold", colour = "black"), # or just size = 9
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "", y = "DNAm of L1 CGIs", fill = "In PMDs:")

ggsave(paste0(outdir, "/L1CGI_boxplot_R4_similarToFig1PlotBC_addP_addBMPC.pdf"),
  plot = p3_addP, width = 8, height = 7, units = "cm",
  limitsize = FALSE
)

all_activeL1 <- rbind(activeL1_sense, activeL1_antisense)
all_activeL1$active <- rep("Yes", nrow(all_activeL1))
all_activeL1$activationScore <- ifelse(str_detect(all_activeL1$group, "sense"),
  all_activeL1$Ratio_Sense,
  all_activeL1$Ratio_Antisense
)

p4_addP <- cpgIsland_L1_meanDNAm_df %>%
  dplyr::left_join(all_activeL1[, c("sampleID", "Geneid", "active")],
    by = c(
      "sample" = "sampleID",
      "te_id" = "Geneid"
    )
  ) %>%
  # dplyr::filter(te_id %in%
  #                 c(activeL1_antisense$Geneid, activeL1_sense$Geneid)) %>%
  dplyr::mutate(active = ifelse(is.na(active), "No", active)) %>%
  ggplot2::ggplot(aes(
    x = reorder(sample, cpgIsland_DNAm),
    y = cpgIsland_DNAm, fill = active
  )) +
  # ggplot2::ggplot(aes(x = reorder(group, -cpgIsland_DNAm),
  #                     y=cpgIsland_DNAm, fill=active))+
  # facet_wrap(~overlap_PMD,
  #            labeller = as_labeller(c(
  #              "Yes" = "In PMDs",
  #              "No" = "Not in PMDs"
  #            )))+
  geom_boxplot(outlier.shape = NA) +
  theme_pubr(base_size = 9) + # Set base font size
  # guides(colour = "In PMDs:") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  theme(
    # legend.position = "none",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.key.size = unit(0.4, "lines"), # make the color square smaller
    # legend.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9, face = "bold", colour = "black"), # or just size = 9
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "", y = "DNAm of L1 CGIs")

ggsave(paste0(outdir, "/activeL1CGI_boxplot_R4_similarToFig1PlotBC_addP_addBMPC.pdf"),
  plot = p4_addP, width = 8, height = 7, units = "cm",
  limitsize = FALSE
)


dnAm_col <- "cpgIsland_DNAm" # 或 "cpgIsland_3p1000_DNAm"
anno2 <- cpgIsland_L1_meanDNAm_df %>%
  dplyr::left_join(all_activeL1[, c("sampleID", "Geneid", "active")],
    by = c(
      "sample" = "sampleID",
      "te_id" = "Geneid"
    )
  ) %>%
  dplyr::mutate(active = ifelse(is.na(active), "No", active)) %>%
  dplyr::filter(group == "PMM", overlap_PMD %in% c("No", "Yes")) %>%
  dplyr::group_by(overlap_PMD, active) %>%
  dplyr::summarise(m = mean(.data[[dnAm_col]], na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = active, values_from = m) %>%
  dplyr::mutate(
    label = sprintf("%.3f", `No` - `Yes`),
    x = "PMM",
    y = 1.1
  )


p_add_addP <- cpgIsland_L1_meanDNAm_df %>%
  dplyr::left_join(all_activeL1[, c("sampleID", "Geneid", "active")],
    by = c(
      "sample" = "sampleID",
      "te_id" = "Geneid"
    )
  ) %>%
  # dplyr::filter(te_id %in%
  #                 c(activeL1_antisense$Geneid, activeL1_sense$Geneid)) %>%
  dplyr::mutate(active = ifelse(is.na(active), "No", active)) %>%
  # dplyr::filter(active == "Yes") %>%
  ggplot2::ggplot(aes(
    x = reorder(group, cpgIsland_DNAm),
    y = cpgIsland_DNAm, fill = active
  )) +
  # ggplot2::ggplot(aes(x = reorder(group, -cpgIsland_DNAm),
  #                     y=cpgIsland_DNAm, fill=active))+
  facet_wrap(~overlap_PMD,
    labeller = as_labeller(c(
      "Yes" = "In PMDs",
      "No" = "Not in PMDs"
    ))
  ) +
  geom_boxplot(outlier.shape = NA) +
  theme_pubr(base_size = 9) + # Set base font size
  # guides(colour = "In PMDs:") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  geom_text(
    data = anno2, aes(x = x, y = y, label = label),
    inherit.aes = FALSE, size = 3
  ) +
  theme(
    # legend.position = "none",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.key.size = unit(0.4, "lines"), # make the color square smaller
    # legend.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9, face = "bold", colour = "black"), # or just size = 9
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "", y = "DNAm of L1 CGIs") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.09)))


ggsave(paste0(outdir, "/activeL1CGI_boxplot_R4_similarToFig1PlotBC_addP_addDelta_addBMPC.pdf"),
  plot = p_add_addP, width = 8, height = 7, units = "cm",
  limitsize = FALSE
)

# L1 CGI DNAm and L1 number cor: ###############################################
p4_corP <- cpgIsland_L1_meanDNAm_df %>%
  dplyr::left_join(
    all_activeL1[, c(
      "sampleID", "Geneid",
      "active", "activationScore"
    )],
    by = c(
      "sample" = "sampleID",
      "te_id" = "Geneid"
    )
  ) %>%
  # dplyr::filter(te_id %in%
  #                 c(activeL1_antisense$Geneid, activeL1_sense$Geneid)) %>%
  dplyr::mutate(active = ifelse(is.na(active), "No", active)) %>%
  dplyr::filter(active == "Yes") %>%
  dplyr::left_join(subset(PMM_sampleSheet, select = -group),
    by = c("sample" = "sampleName")
  ) %>%
  dplyr::filter(sample %in% PMM_sampleSheet$sampleName) %>%
  group_by(sample) %>%
  summarise(
    cpgIsland_DNAm = mean(cpgIsland_DNAm),
    L1 = mean(L1),
    activationScore = mean(activationScore)
  ) %>%
  ggplot2::ggplot(aes(
    x = L1,
    # x = activationScore,
    y = cpgIsland_DNAm
  )) +
  # ggplot2::ggplot(aes(x = reorder(group, -cpgIsland_DNAm),
  #                     y=cpgIsland_DNAm, fill=active))+
  # facet_wrap(~overlap_PMD,
  #            labeller = as_labeller(c(
  #              "Yes" = "In PMDs",
  #              "No" = "Not in PMDs"
  #            )))+
  # geom_boxplot(outlier.shape = NA)+
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y = c(1, 1), size = 3) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3) +
  theme_pubr(base_size = 9.6) + # Set base font size
  # guides(colour = "In PMDs:") +
  # stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE)+
  theme(
    # legend.position = "none",
    legend.text = element_text(size = 9.6),
    legend.title = element_text(size = 9.6),
    legend.key.size = unit(0.4, "lines"), # make the color square smaller
    # legend.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9.6, face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "Number of active L1s", y = "Mean DNAm of L1 CGIs")

ggsave(paste0(outdir, "/activeL1CGI_L1num_cor_addP_addBMPC.pdf"),
  plot = p4_corP, width = 8, height = 7, units = "cm",
  limitsize = FALSE
)
