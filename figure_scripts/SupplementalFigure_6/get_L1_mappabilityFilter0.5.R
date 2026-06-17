library(tidyverse)
library(plyranges)
library(rtracklayer)
library(data.table)
library(ggplot2)
library(ggpubr)

# 1. read in activated and inactivated L1s ========
activeL1_sense <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/activationScore/filtered_L1/sense_L1s.tsv")
activeL1_antisense <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/activationScore/filtered_L1/antisense_L1s.tsv")
all_activeL1 <- rbind(activeL1_sense, activeL1_antisense) %>% 
  dplyr::mutate(direction = str_remove_all(group, "_L1score")) %>% 
  dplyr::mutate(direction = str_remove_all(direction, "..+_"))
all_activeL1$active <- rep("Yes", nrow(all_activeL1))

map <- c(
  "CD19-GSM2171212" = "BMPC4",
  "CD19-GSM2171213" = "BMPC5",
  "CD19-GSM2171214" = "BMPC6",
  "CD19-GSM2171215" = "BMPC7")
all_activeL1$sampleID <- str_replace_all(all_activeL1$sampleID, map)

# filter for activated L1 in high mappability regions
mappability_FL_L1 <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/hg38_mappability/fulllength_L1s_mean_mappability.tsv")

mappability_FL_L1_fil <- mappability_FL_L1 %>% 
  dplyr::filter(str_detect(track, "umap")) %>% 
  dplyr::filter(mean_mappability>0.5)

# Reorder sampleID by total count before plotting
activeL1_all_forPlot <- all_activeL1 %>%
  dplyr::filter(Geneid %in% mappability_FL_L1_fil$te_id) %>% 
  count(sampleID) %>%
  right_join(all_activeL1, by = "sampleID") %>%
  mutate(sampleID = fct_reorder(sampleID, n))  %>% # or reorder(sampleID, n)
  dplyr::mutate(
    repName_less = ifelse(
      repName %in% c("L1HS", "L1PA2", "L1PA3", "L1PA4", 
                     "L1PA5", "L1PA6", "L1PA7", "L1PA8"), 
      repName, "others"))

fig1_plot5 <- activeL1_all_forPlot %>%
  ggplot(aes(x = sampleID, fill = repName_less)) +  
  geom_bar() +
  # Make facet labels UPPERCASE
  # facet_wrap(~direction, labeller = labeller(direction = tools::toTitleCase)) +  
  coord_flip() +
  scale_fill_viridis_d(option = "turbo") +
  guides(fill = guide_legend(ncol = 2)) +   #2-column legend
  theme_pubr(base_size = 9, legend = c(0.8, 0.3)) +  # Set base font size
  theme(
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.key.size = unit(0.4, "lines"),     # make the color square smaller
    legend.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9, face = "bold")  # or just size = 9
  ) +
  labs(x = "", y = "Num. Active", fill = "Repeat names")

ggsave(paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/hg38_mappability", "/num_active_all_L1s_pmm_barplot_fig1E_mappFil.pdf"), 
       plot = fig1_plot5, width = 10, height = 8.5, units = "cm", 
       limitsize = FALSE)

fig2_plotB <- activeL1_all_forPlot %>%
  ggplot(aes(x = sampleID, fill = repName_less)) +  
  geom_bar() +
  # Make facet labels UPPERCASE
  facet_wrap(~direction, labeller = labeller(direction = tools::toTitleCase)) +
  coord_flip() +
  scale_fill_viridis_d(option = "turbo") +
  guides(fill = guide_legend(ncol = 2)) +   #2-column legend
  theme_pubr(base_size = 9, legend = c(0.8, 0.3)) +  # Set base font size
  theme(
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.key.size = unit(0.4, "lines"),     # make the color square smaller
    legend.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9, face = "bold")  # or just size = 9
  ) +
  labs(x = "", y = "Num. Active", fill = "Repeat names")

ggsave(paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/hg38_mappability", "/num_active_all_L1s_pmm_barplot_fig2B_mappaFil.pdf"), 
       plot = fig2_plotB, width = 10, height = 8.5, units = "cm", 
       limitsize = FALSE)

# supFig2_A
S2_plotA <- activeL1_all_forPlot %>%
  dplyr::mutate(group = ifelse(str_detect(sampleID, "BMPC"), "BMPC", "PMM")) %>% 
  dplyr::mutate(activationScore = ifelse(direction=="sense", Ratio_Sense, Ratio_Antisense)) %>% 
  # ggplot(aes(x = repName_less , y = activationScore, fill = group)) +  
  ggplot(aes(x = group , y = activationScore, fill = group)) +  
  # geom_boxplot(outlier.shape = NA)+
  geom_boxplot()+
  stat_compare_means(method = "wilcox.test", label = "p.signif")+
  # Make facet labels UPPERCASE
  guides(fill = "none") +   #2-column legend
  theme_pubr(base_size = 9, legend = c(0.8, 0.3)) +  # Set base font size
  theme(
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.key.size = unit(0.4, "lines"),     # make the color square smaller
    legend.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9, face = "bold")  # or just size = 9
  ) +
  labs(x = "", y = "Activation score")

ggsave(paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/hg38_mappability", "/num_active_all_L1s_pmm_boxplot_figS2A_mappaFil.pdf"), 
       plot = S2_plotA, width = 10, height = 8.5, units = "cm", 
       limitsize = FALSE)




# also filter for >0.5 mappability active L1s
all_activeL1_PMM <- all_activeL1 %>%
  dplyr::filter(Geneid %in% mappability_FL_L1_fil$te_id) %>% 
  dplyr::filter(str_detect(sampleID, "PMM")) %>%
  dplyr::select(
    te_id = Geneid, repName, repClass, repFamily,
    seqnames = genoName, start = genoStart, end = genoEnd,
    strand) %>%
  unique() %>% 
  as_granges() %>% 
  plyranges::filter(repName %in% c("L1PA2"))
  # plyranges::filter(repName %in% c("L1HS", "L1PA2", "L1PA3")) 


# all L1 repeat:
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
allTE <- 
  repeat_bed_df %>%
  dplyr::mutate(seqnames = genoName,
                start = genoStart + 1, end = genoEnd) %>%
  dplyr::select(-c(genoName:genoEnd)) %>%
  unique() %>%
  as_granges() 

repeat_GR_inact_L1 <-
  repeat_bed_df %>%
  dplyr::mutate(seqnames = genoName,
    start = genoStart + 1, end = genoEnd) %>%
  dplyr::select(-c(genoName:genoEnd)) %>%
  unique() %>%
  as_granges() %>%
  plyranges::filter(repFamily == "L1") %>%
  plyranges::filter(repName %in% c("L1PA2")) %>%
  # plyranges::filter(repName %in% c("L1HS", "L1PA2", "L1PA3")) %>% 
  plyranges::filter(width(.) > 5900) %>%
  plyranges::filter(!(te_id %in% all_activeL1_PMM$te_id)) %>% 
  dplyr::filter(te_id %in% mappability_FL_L1_fil$te_id) 

# get the L1 age cut off
repeat_GR_inact_L1_PMMage <-
  repeat_bed_df %>%
  dplyr::mutate(
    seqnames = genoName,
    start = genoStart + 1, end = genoEnd
  ) %>%
  dplyr::select(-c(genoName:genoEnd)) %>%
  unique() %>%
  as_granges() %>%
  plyranges::filter(repName %in% c("L1PA2")) %>%
  # plyranges::filter(repName %in% c("L1HS", "L1PA2", "L1PA3")) %>%
  plyranges::filter(te_id %in% all_activeL1_PMM$te_id)

repeat_GR_inact_L1 <- repeat_GR_inact_L1[repeat_GR_inact_L1$age <= max(repeat_GR_inact_L1_PMMage$age),]

# 2. open mappability bw files as indexed handles (no full-genome load) =============
umap_bw_path   <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_mappability/hg38_umap_k150_out/hg38_umap_k150_multi_read.bw"
bismap_bw_path <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_mappability/hg38_bismap_k150_out/hg38_bismap_k150_multi_read.bw"

umap_bw   <- rtracklayer::BigWigFile(umap_bw_path)
bismap_bw <- rtracklayer::BigWigFile(bismap_bw_path)

# 3. compute mean mappability per full-length L1 =============
#    import.bw(which=gr) loads only the queried L1 intervals (fast, low-memory)
#    then overlap-aggregate scores back to each L1 element
calc_mean_mappability <- function(bw_file, gr) {
  bw_data <- rtracklayer::import.bw(bw_file, which = gr)          # only reads relevant regions
  hits    <- findOverlaps(gr, bw_data)
  agg     <- tapply(bw_data$score[subjectHits(hits)],
                    queryHits(hits),
                    mean, na.rm = TRUE)
  result  <- rep(NA_real_, length(gr))
  result[as.integer(names(agg))] <- as.numeric(agg)              # preserve original L1 order
  result
}

message("Computing mean umap mappability for active L1s ...")
umap_active <- calc_mean_mappability(umap_bw, all_activeL1_PMM)

message("Computing mean umap mappability for inactive full-length L1s ...")
umap_inact  <- calc_mean_mappability(umap_bw, repeat_GR_inact_L1)
# umap_inact  <- calc_mean_mappability(umap_bw, repeat_GR_inact_L1_2000)

message("Computing mean bismap mappability for active L1s ...")
bismap_active <- calc_mean_mappability(bismap_bw, all_activeL1_PMM)

message("Computing mean bismap mappability for inactive full-length L1s ...")
bismap_inact  <- calc_mean_mappability(bismap_bw, repeat_GR_inact_L1)
# bismap_inact  <- calc_mean_mappability(bismap_bw, repeat_GR_inact_L1_2000)

# combine into tidy long table
df_active <- data.frame(
  te_id       = all_activeL1_PMM$te_id,
  active      = "Yes",
  umap_mean   = umap_active,
  bismap_mean = bismap_active
)

df_inact <- data.frame(
  te_id       = repeat_GR_inact_L1$te_id,
  active      = "No",
  umap_mean   = umap_inact,
  bismap_mean = bismap_inact
)

df_mappability <- rbind(df_active, df_inact)

# pivot to long for plotting (one row per L1 × track)
df_long <- df_mappability %>%
  tidyr::pivot_longer(
    cols      = c(umap_mean, bismap_mean),
    names_to  = "track",
    values_to = "mean_mappability"
  ) %>%
  dplyr::mutate(
    track = dplyr::recode(track,
      umap_mean   = "umap (RNAseq, k150)",
      bismap_mean = "bismap (WGBS, k150)"
    )
  )

message("Active L1s   : ", nrow(df_active),   " elements")
message("Inactive L1s : ", nrow(df_inact),    " elements")

# =========================
# 4. plot: active vs inactive L1 mappability, faceted by track
# =========================
p <- ggplot(df_long, aes(x = active, y = mean_mappability, fill = active)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
  # ggpubr::stat_compare_means(
  #   method = "wilcox.test",
  #   label  = "p.format",
  #   vjust  = 1
  # ) +
  ggpubr::stat_compare_means(
    method = "t.test",
    label  = "p.format",
    vjust  = 1
  )+
  facet_wrap(~track, scales = "free_y") +
  ggpubr::theme_pubr(base_size = 9.6) +
  labs(
    x    = "Transcriptionally active L1PA2",
    y    = "Mean mappability (multi-read)",
    fill = "Active L1PA2"
  )

print(p)

# =========================
# 5. 输出结果表
# =========================
out_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/hg38_mappability"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  file.path(out_dir, "fulllength_L1PA2_mappability0.5Fil.pdf"),
  plot = p,
  width = 14, height = 8, units = "cm", limitsize = FALSE
)

out_tsv <- file.path(out_dir, "fulllength_L1PA2_mean_mappability0.5Fil.tsv")
data.table::fwrite(df_long, out_tsv, sep = "\t")
message("Done. Wrote: ", out_tsv)


# check genetic or inter generic
# all_activeL1_antisense_PMM <- all_activeL1 %>%
#   dplyr::filter(str_detect(sampleID, "PMM")) %>%
#   dplyr::filter(str_detect(group, "_antisense_")) %>% 
#   dplyr::select(
#     te_id = Geneid, repName, repClass, repFamily,
#     seqnames = genoName, start = genoStart, end = genoEnd,
#     strand
#   ) %>% 
#   unique() %>% 
#   as_granges() 
# 
# library(GenomicFeatures)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# 
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# genes_gr <- genes(txdb)
# 
# all_activeL1_antisense_PMM$genic <-
#   countOverlaps(all_activeL1_antisense_PMM, genes_gr) > 0


