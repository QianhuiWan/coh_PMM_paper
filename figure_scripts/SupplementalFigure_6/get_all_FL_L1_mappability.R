library(tidyverse)
library(plyranges)
library(rtracklayer)
library(data.table)

# 1. read in activated and inactivated L1s ========
activeL1_sense <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/activationScore/filtered_L1/sense_L1s.tsv")
activeL1_antisense <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/activationScore/filtered_L1/antisense_L1s.tsv")
all_activeL1 <- rbind(activeL1_sense, activeL1_antisense)
all_activeL1$active <- rep("Yes", nrow(all_activeL1))

map <- c(
  "CD19-GSM2171212" = "BMPC4",
  "CD19-GSM2171213" = "BMPC5",
  "CD19-GSM2171214" = "BMPC6",
  "CD19-GSM2171215" = "BMPC7")
all_activeL1$sampleID <- str_replace_all(all_activeL1$sampleID, map)

all_activeL1_PMM <- all_activeL1 %>%
  dplyr::filter(str_detect(sampleID, "PMM")) %>%
  dplyr::select(
    te_id = Geneid, repName, repClass, repFamily,
    seqnames = genoName, start = genoStart, end = genoEnd,
    strand) %>%
  unique() %>% 
  as_granges() #%>% # n=625
  # plyranges::filter(repName %in% c("L1PA2")) 
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

repeat_GR_FL_L1 <-
  repeat_bed_df %>%
  dplyr::mutate(seqnames = genoName,
    start = genoStart + 1, end = genoEnd) %>%
  dplyr::select(-c(genoName:genoEnd)) %>%
  unique() %>%
  as_granges() %>%
  plyranges::filter(repFamily == "L1") %>% 
  plyranges::filter(width(.) > 5900)

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

calc_mean_mappability <- function(bw_file, gr, 
                                  score_col = "mean_mappability") {
  bw_data <- rtracklayer::import.bw(bw_file, which = gr)
  hits <- findOverlaps(gr, bw_data)
  
  agg <- tapply(
    bw_data$score[subjectHits(hits)],
    queryHits(hits),
    mean,
    na.rm = TRUE)
  
  result <- rep(NA_real_, length(gr))
  result[as.integer(names(agg))] <- as.numeric(agg)
  
  gr_df <- as.data.frame(gr)
  gr_df[[score_col]] <- result
  
  gr_df
}

message("Computing mean umap mappability for active L1s ...")
umap_active_PMM <- calc_mean_mappability(umap_bw, all_activeL1_PMM)

message("Computing mean umap mappability for all TE and full length L1s ...")
# umap_allTE  <- calc_mean_mappability(umap_bw, allTE)
umap_FL_L1  <- calc_mean_mappability(umap_bw, repeat_GR_FL_L1)

message("Computing mean bismap mappability for active L1s ...")
bismap_active_PMM <- calc_mean_mappability(bismap_bw, all_activeL1_PMM)

message("Computing mean bismap mappability for all TE and full length L1s ...")
# bismap_allTE  <- calc_mean_mappability(bismap_bw, allTE)
bismap_FL_L1  <- calc_mean_mappability(bismap_bw, repeat_GR_FL_L1)

# combine into tidy long table
df_active <- umap_active_PMM %>% 
  dplyr::select(te_id, umap_mean = mean_mappability) %>% 
  dplyr::left_join(bismap_active_PMM, by = "te_id") %>% 
  dplyr::select(te_id, umap_mean, bismap_mean = mean_mappability) #%>% 
  # dplyr::mutate(active = rep("Yes", nrow(.)))
  
df_FL_L1 <- umap_FL_L1 %>% 
  dplyr::select(te_id, umap_mean = mean_mappability) %>% 
  dplyr::left_join(bismap_FL_L1, by = "te_id") %>% 
  dplyr::select(te_id, umap_mean, bismap_mean = mean_mappability, repName) %>% 
  dplyr::mutate(repName_less = ifelse(
    repName %in% c("L1HS", "L1PA2", "L1PA3", "L1PA4", 
                   "L1PA5", "L1PA6", "L1PA7", "L1PA8"), 
    repName, "others"))


# pivot to long for plotting (one row per L1 × track)
df_long_active <- df_active %>%
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

df_long <- df_FL_L1 %>%
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

# =========================
# 4. plot: active vs inactive L1 mappability, faceted by track
# =========================

p_active_L1 <- ggplot(df_long_active, aes(x = track, 
                                  y = mean_mappability
)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey") +
  ggpubr::theme_pubr(base_size = 9.6) +
  labs(
    x    = "All transcriptional active L1s",
    y    = "Mean mappability (multi-read)"
  )

print(p_active_L1)


p_allFL_L1 <- ggplot(df_long, aes(x = repName_less, 
                                  y = mean_mappability
                                  )) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey") +
  # geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
  # ggpubr::stat_compare_means(
  #   method = "wilcox.test",
  #   label  = "p.format",
  #   vjust  = 1
  # ) +
  # ggpubr::stat_compare_means(
  #   method = "t.test",
  #   label  = "p.format",
  #   vjust  = 1
  # )+
  facet_wrap(~track, scales = "free_y", nrow = 2) +
  ggpubr::theme_pubr(base_size = 9.6) +
  labs(
    x    = "All full length L1s",
    y    = "Mean mappability (multi-read)"
  )

print(p_allFL_L1)

p_allFL_L12 <- df_long %>% 
  dplyr::filter(str_detect(track, "umap")) %>% 
  ggplot(aes(x = track, 
             y = mean_mappability)) +
  geom_violin(
    fill = "grey85",
    color = "grey40",
    width = 0.8,
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.12,
    fill = "white",
    color = "black",
    outlier.shape = NA
  ) +
  ggpubr::theme_pubr(base_size = 9.6) +
  labs(
    x    = "All full length L1s",
    y    = "Mean mappability (multi-read)"
  )

print(p_allFL_L12)

# =========================
# 5. get results
# =========================
out_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/hg38_mappability"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  file.path(out_dir, "all_active_L1_mappability.pdf"),
  plot = p_active_L1,
  width = 14, height = 8, units = "cm", limitsize = FALSE
)

ggsave(
  file.path(out_dir, "all_fulllength_L1_mappability.pdf"),
  plot = p_allFL_L1,
  width = 14, height = 12, units = "cm", limitsize = FALSE
)

ggsave(
  file.path(out_dir, "all_fulllength_L1_mappability_violin.pdf"),
  plot = p_allFL_L12,
  width = 6, height = 8, units = "cm", limitsize = FALSE
)

out_tsv <- file.path(out_dir, "fulllength_L1s_mean_mappability.tsv")
data.table::fwrite(df_long, out_tsv, sep = "\t")
message("Done. Wrote: ", out_tsv)


