
# R3: only get the KLF binding at the >4500nt L1 from L1base V2

# args <- commandArgs(trailingOnly = TRUE)
# cell_type <- args[1]
cellType <- "K562"

outdir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/L1_TFs"
# n_bins <- c(1000, 500, 50, 2000, 500)

library(tidyverse)
library(GenomicRanges)
library(plyranges)
library(IRanges)
library(purrr)

library(ggplot2)
library(dplyr)

# step1: add TF (including KLFs) binding coordinates to L1 overlapped TFs
# ~/githubRepo/coh_PMM/RNAseq_scripts/L1_ORF1_5UTR_predict/use_encode_chipSeq/detect_L1bindTFs_R2.R
L1_bindingTFs_cellType <- read_tsv(paste0(outdir, "/all_overlaps_L1_binding_TFs_", 
                                          cellType, ".tsv"))
## get the GR object
L1_bindingTFs_cellType_GR <- L1_bindingTFs_cellType %>% 
  separate(TFbindingPos_coord, into = c("TF_chr", "TF_range"), sep = ":") %>%
  separate(TF_range, into = c("TF_start", "TF_end"), sep = "-") %>%
  mutate(
    TF_start = as.numeric(TF_start),
    TF_end = as.numeric(TF_end)) %>%
  dplyr::select(seqnames = TF_chr, start = TF_start, end = TF_end, 
                TF_id, TF_name, name:peakPosition, 
                TE_chr = seqnames, TE_start = start, TE_end = end, 
                TE_width = width, TE_strand = strand, te_id:age) %>% 
  as_granges()

## filter KLF only 
L1_bindingTFs_cellType_KLF_GR <- L1_bindingTFs_cellType_GR %>% 
  plyranges::filter(str_detect(TF_name, "KLF"))

# step2: add ORF annotation to full length L1s
## note: full length L1 may contain several TE annotations from rmsk
FL_L1_anno <- read_csv(paste0(outdir,"/hg38_Ens84_38_L1_above4500nt_annotation.csv"))


## then need to get coord for L1 parts:
# ---- Load L1 annotation TSV: change to GR coord for 5p and 3p of L1 ----
FL_L1_anno_df <- FL_L1_anno %>%
  # dplyr::mutate(Chr = paste0("chr", Chr)) %>% 
  mutate(
    strand_chr = case_when(
      Strand == -1 ~ "-",
      Strand == 1 ~ "+",
      TRUE ~ "*"
    ),
    L1_5prime = if_else(strand_chr == "+", Start, End),
    L1_3prime = if_else(strand_chr == "+", End, Start),
    L1_5prime =L1_5prime + 1 # Convert BED-style start to 1-based inclusive
  )

# ---- Function to convert relative positions to genomic ----
relative_to_genomic <- function(l1_start, strand, rel_start, rel_end) {
  if (strand == "+") {
    genomic_start <- l1_start + rel_start
    genomic_end   <- l1_start + rel_end
  } else if (strand == "-") {
    genomic_start <- l1_start - rel_end
    genomic_end   <- l1_start - rel_start
  } else {
    stop("Strand must be '+' or '-'")
  }
  list(
    start = min(genomic_start, genomic_end),
    end = max(genomic_start, genomic_end)
  )
}


# ---- Helper to convert any feature ----
make_feature_gr <- function(df, rel_start_col, rel_end_col, region_name) {
  pmap_dfr(df, function(...) {
    row <- list(...)
    coords <- relative_to_genomic(
      l1_start = row$Start, # original coord, not bed format coord
      strand = row$strand_chr,
      rel_start = row[[rel_start_col]],
      rel_end = row[[rel_end_col]]
    )
    tibble(
      seqnames = paste0("chr", row$Chr),
      start = coords$start+1,
      end = coords$end,
      strand = row$strand_chr,
      ID = row$ID,
      region = region_name
    )
  }) %>% 
    as_granges()
}

# ---- Create GRanges for all features ----

# ORF1
orf1_gr <- make_feature_gr(df=FL_L1_anno_df, 
                           rel_start_col = "ORF1 Start", 
                           rel_end_col = "ORF1 End", region_name = "ORF1")

# ORF2
# orf2_gr <- make_feature_gr(FL_L1_anno_df,  rel_start_col = "ORF2 Start", 
#                            rel_end_col = "ORF2 End", region_name = "ORF2")

# update full length L1 annotation
FL_L1_anno_GR_df <- FL_L1_anno_df %>% 
  dplyr::select(ID, Chr:End, Strand, `Intactness Score`, 
                # `ORF1 Start`, `ORF1 End`, `ORF2 Start`, `ORF2 End`,
                L1_5prime,L1_3prime, strand_chr) %>% 
  left_join(as.data.frame(orf1_gr), by="ID") %>% 
  dplyr::mutate(ORF1_start = start, ORF1_end = end, ORF1 = region) %>% 
  dplyr::select(-(seqnames:region))
  # left_join(as.data.frame(orf2_gr), by="ID") %>% 
  # dplyr::mutate(ORF2_start = start, ORF2_end = end, ORF2 = region) %>% 
  # dplyr::select(-(seqnames:region))

FL_L1_anno_GR <- FL_L1_anno_GR_df %>% 
  dplyr::mutate(seqnames = paste0("chr", Chr), 
                start = L1_5prime, 
                end = L1_3prime,
                strand = strand_chr) %>% 
  dplyr::select(ID, seqnames:strand) %>% 
  as_granges()

# 5'UTR: from L1_5prime to ORF1 start - 1 for + strand
utr5_gr <- FL_L1_anno_GR_df %>% 
  dplyr::mutate(seqnames = paste0("chr", Chr), 
                start = ifelse(strand_chr=="+", pmin(L1_5prime, ORF1_start - 1), pmin(ORF1_end + 1, L1_3prime)), 
                end = ifelse(strand_chr=="+", pmax(L1_5prime, ORF1_start - 1), pmax(ORF1_end + 1, L1_3prime)),
                strand = strand_chr,
                region = rep("5UTR", nrow(.))) %>% 
  dplyr::select(ID, seqnames:region) %>% 
  as_granges()

# Linker sequence to 3UTR : ORF1 end + 1 to the end of L1 for + strand
linkerToEnd_gr <- FL_L1_anno_GR_df %>% 
  dplyr::mutate(seqnames = paste0("chr", Chr), 
                start = ifelse(strand_chr=="+", pmin(ORF1_end + 1, L1_3prime), pmin(L1_5prime, ORF1_start - 1)), 
                end = ifelse(strand_chr=="+", pmax(ORF1_end + 1, L1_3prime), pmax(L1_5prime, ORF1_start - 1)),
                strand = strand_chr,
                region = rep("afterORF1", nrow(.))) %>% 
  dplyr::select(ID, seqnames:region) %>% 
  as_granges()

# ---- Combine all features into one GRanges ----
FL_L1_features_gr <- c(utr5_gr, orf1_gr, linkerToEnd_gr)

# ---- Done! ----
# You now have: l1_features_gr (GRanges with all features + region labels)
# Optional: save as RDS or tsv
# saveRDS(l1_features_gr, "L1_feature_GRanges.rds")
write_tsv(as.data.frame(FL_L1_features_gr), 
          file = paste0(outdir, "/L1_above4500nt_features_gr_df.tsv"))







# step3: Now we can bin the MACS2 values by regions ############################
## Load RepeatMasker BED file
repeat_bed_path <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_repeatMasker_hg38_addTEid.tsv"
repeat_ranges <- read_tsv(repeat_bed_path, 
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
                          )) %>%
  dplyr::mutate(seqnames = genoName, 
                start = genoStart + 1, end = genoEnd) %>%
  dplyr::select(-c(genoName:genoEnd)) %>% 
  as_granges()

FL_L1_bindingTF_GR <- 
  find_overlaps(FL_L1_anno_GR,L1_bindingTFs_cellType_KLF_GR) #425



# ---- Helper: Bin and summarize TF signal per feature ----
bin_feature_signal <- function(feature_gr, tf_gr, nbins, region_name, tf_name) {
  message("Binning ", region_name, " ", tf_name, " MACS2 into ",
          nbins, " bins...")
  
  # Tile into bins
  feature_binned <- tile(feature_gr, n = nbins) %>% unlist()
  
  # Assign bin metadata
  feature_binned$bin_id <- rep(seq_len(nbins), times = length(feature_gr))
  feature_binned$ID <- rep(feature_gr$ID, each = nbins)
  feature_binned$region <- region_name
  
  # Overlap with TF signal
  result <- feature_binned %>%
    join_overlap_left(tf_gr) %>%
    as_tibble() %>%  # 👈 Convert GRanges to data frame (with all metadata)
    group_by(ID, region, bin_id) %>%
    summarise(avg_signal = mean(singnalValue, na.rm = TRUE), .groups = "drop") %>% 
    mutate(TF = tf_name)
  return(result)
}

# ---- Parameters ----
# features_list <- list(utr5_gr, orf1_gr, linkerToEnd_gr)
feature_names <- c("5UTR", "ORF1", "afterORF1")
bin_counts <- c(500, 250, 1000)  # <-- your custom bin counts

utr5_gr_fil <- utr5_gr[width(utr5_gr)>=500]
orf1_gr_fil <- orf1_gr[width(orf1_gr)>=250]
linkerToEnd_gr_fil <- linkerToEnd_gr[width(linkerToEnd_gr)>=1000]
features_list <- list(utr5_gr_fil, orf1_gr_fil, linkerToEnd_gr_fil)

# TF list
tf_list <- split(L1_bindingTFs_cellType_KLF_GR, L1_bindingTFs_cellType_KLF_GR$TF_name)

# features needed for function
feature_df <- tibble(
  gr = features_list,
  nbins = bin_counts,
  region_name = feature_names
)


# main loop, each TF feature bin ###############################################
bin_results <- map_dfr(names(tf_list), function(tf_name) {
  tf_gr <- tf_list[[tf_name]]
  
  pmap_dfr(feature_df, function(gr, nbins, region_name) {
    bin_feature_signal(gr, tf_gr, nbins, region_name, tf_name)
  })
})



# Define standard chromosome names
# # Collect all your feature GRs and the TF peak GR
# L1_bindingTFs_cellType_KLF_GR <- L1_bindingTFs_cellType_GR
# all_grs <- list(utr5_gr, orf1_gr, linker_gr, orf2_gr, utr3_gr, 
#                 L1_bindingTFs_cellType_KLF_GR)
# # Find chromosomes they all share
# shared_chrs <- purrr::reduce(map(all_grs, seqlevels), intersect)
# # Optional: Print for debugging
# message("Shared chromosomes: ", paste(shared_chrs, collapse = ", "))
# 
# # Define a custom cleaning function that removes chrY (or anything not shared)
# clean_gr_shared <- function(gr, keep_set=shared_chrs) {
#   shared <- intersect(seqlevels(gr), keep_set)
#   gr %>%
#     keepSeqlevels(shared, pruning.mode = "coarse") %>%
#     sortSeqlevels()
# }
# 
# # Clean TF peak GRanges
# TF_gr <- clean_gr(L1_bindingTFs_cellType_KLF_GR)
# # Clean all feature GRanges
# features_list <- map(features_list, clean_gr)


# ---- Normalize bin ID (0–1 scale) for plotting ----
bin_results <- bin_results %>%
  group_by(region) %>%
  mutate(norm_bin = bin_id / max(bin_id)) %>%
  ungroup()

# ---- Save output ----
write_tsv(bin_results, "binned_TF_signal_custom_bins.tsv")
message("Binning complete: saved as 'binned_TF_signal_custom_bins.tsv'")

# ---- Plot mean signal across normalized bins ----
# avg_profile <- bin_results %>%
#   dplyr::filter(TF %in% c("KLF6")) %>% 
#   group_by(region, norm_bin) %>%
#   summarize(mean_signal = mean(avg_signal, na.rm = TRUE), .groups = "drop")

bin_results %>% 
  dplyr::filter(TF %in% c("KLF6")) %>% 
  dplyr::mutate(region = factor(region, levels = c("5UTR", "ORF1", "afterORF1"))) %>% 
  ggplot(aes(x = bin_id, y = avg_signal, color = region)) +
  geom_point()+
  geom_smooth()+
  # geom_line(size = 1.2) +
  facet_wrap(~region, scales = "free") +
  theme_minimal(base_size = 14) +
  labs(
    # title = "Average TF Signal Across L1 Features (normalized bins)",
    title = "Average TF Signal Across L1 Features",
    # x = "Normalized Position (0–1)",
    y = "Average MACS2 Signal",
    color = "Region"
  ) +
  theme(legend.position = "right")







# Step 1: Normalization bin_id → bin_prop ################################
bin_results_plot <- bin_results %>%
  group_by(region, ID) %>%
  mutate(bin_prop = bin_id / max(bin_id)) %>%
  ungroup()

df_clean <- bin_results_plot %>%
  filter(TF == "KLF16") %>%
  # filter(TF == "KLF6") %>%
  mutate(avg_signal = ifelse(is.na(avg_signal) | is.nan(avg_signal) | is.infinite(avg_signal), 0, avg_signal))

# normalizing function for each KLF
get_density_by_region <- function(df_region, region_name,
                                  adjust = 3,
                                  n = 1000,
                                  normalize = FALSE) {
  df_region <- df_region %>%
    mutate(avg_signal = ifelse(
      is.na(avg_signal) | is.nan(avg_signal) | is.infinite(avg_signal),
      0,
      avg_signal
    ))
  
  total_signal <- sum(df_region$avg_signal)
  x_range <- range(df_region$bin_id)
  x_seq <- seq(x_range[1], x_range[2], length.out = n)
  
  # no signal → flat 0 curve
  if (total_signal == 0) {
    return(data.frame(x = x_seq, y = 0, region = region_name))
  }
  
  # weighted density
  dens <- density(
    x = df_region$bin_id,
    weights = df_region$avg_signal / total_signal,
    from = x_range[1],
    to = x_range[2],
    n = n,
    adjust = adjust
  )
  
  y_scaled <- dens$y * total_signal
  
  if (normalize) {
    y_scaled <- y_scaled / max(y_scaled)
    # y_scaled <- y_scaled / max(avg_signal)
  }
  
  data.frame(x = dens$x, y = y_scaled, region = region_name)
}


density_all_regions <- df_clean %>%
  group_split(region) %>%
  map_dfr(~ get_density_by_region(.x, unique(.x$region), normalize = FALSE)) %>% 
  dplyr::mutate(region = factor(region, levels = c("5UTR", "ORF1", "afterORF1")))


ggplot(density_all_regions, aes(x = x, y = y)) +
  geom_line(color = "#1f77b4", size = 1.2) +
  facet_wrap(~region, scales = "free_x") +
  labs(
    x = "Bin index (region-specific)",
    y = "Normalized TF binding signal",
    # title = "KLF6 binding across L1 sub-regions"
    title = "KLF16 binding across L1 sub-regions"
  ) +
  theme_pubr(base_size = 12)
ggsave(paste0(outdir, "/L1_above4500nt_KLF6_binding.pdf"),
       device = "pdf", width = 15, height = 10, units = "cm")














# ---- Save Output ----
output_file <- paste0(outdir, "/FL_L1_binned_signal_", cellType, ".tsv")
write_tsv(bin_results, output_file)
message("✅ Binned signal data saved to: ", output_file)



