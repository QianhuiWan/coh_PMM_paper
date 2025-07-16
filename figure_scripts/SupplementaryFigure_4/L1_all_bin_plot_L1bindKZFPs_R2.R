
# R2: only get the KZFP binding at the >4500nt L1 from L1base V2 +
# Use MACS scores from chip-exo bed files
# + only focus on L1PA2 binding top20 KZFPs

# args <- commandArgs(trailingOnly = TRUE)
# KZFP_name <- args[1]

hg38_bed_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/chipExo/GSE78099/KZFP_exo_hg38"
indir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/L1_TFs"
outdir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs"

library(tidyverse)
library(readr)
library(GenomicRanges)
library(plyranges)
library(rtracklayer)
library(IRanges)
library(purrr)
library(ggplot2)
library(ggpubr)
library(dplyr)

################################################################################
# step1: add KZFP binding coordinates to L1 overlapped KZFPs
################################################################################
# ~/githubRepo/coh_PMM/RNAseq_scripts/KZFP_metaData/get_L1binding_KZFPs_R4_updatedActivationScore.R
L1binding_KZFPs <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/all_overlaps_L1_binding_KZFPs_PMM_R4.tsv")

## get the GR object for KZFPs
L1binding_KZFPs_GR <- L1binding_KZFPs %>% 
  separate(KZFPbindingPos_coord, into = c("KZFP_chr", "KZFP_range"), sep = ":") %>%
  separate(KZFP_range, into = c("KZFP_start", "KZFP_end"), sep = "-") %>%
  mutate(
    KZFP_start = as.numeric(KZFP_start),
    KZFP_end = as.numeric(KZFP_end)) %>%
  dplyr::select(seqnames = KZFP_chr, start = KZFP_start, end = KZFP_end, 
                KZFP_id, KZFP_name, name:score, 
                TE_chr = seqnames, TE_start = start, TE_end = end, 
                TE_width = width, TE_strand = strand, te_id:age) %>% 
  as_granges()

## filter top20 KZFP only 
L1PA2binding_KZFPs_GR <- L1binding_KZFPs_GR %>% 
  plyranges::filter(repName == "L1PA2")
top20_KZFP_byL1PA2bindingSites <- L1PA2binding_KZFPs_GR$KZFP_name %>% 
  table() %>% sort(decreasing = TRUE) %>% head(20) %>% names()

L1binding_KZFPs_top20_GR <- L1PA2binding_KZFPs_GR %>% 
  plyranges::filter(KZFP_name %in% top20_KZFP_byL1PA2bindingSites)

## check how many KLF binding to L1
table(L1binding_KZFPs_top20_GR$KZFP_name) %>% sort()
### KLF1, KLF10, KLF13, KLF16 and KLF6

# Split into list for downstream binning --
# top20L1binding_by_KZFP <- 
#   split(L1binding_KZFPs_top20_GR, L1binding_KZFPs_top20_GR$KZFP_name)

################################################################################
# step2: add ORF annotation to full length L1s
################################################################################
## note: full length L1 may contain several TE annotations from rmsk
FL_L1_anno <- read_csv(paste0(indir,"/hg38_Ens84_38_L1_above4500nt_annotation.csv"))

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

# write_tsv(as.data.frame(FL_L1_features_gr),
#           file = paste0(indir, "/L1_above4500nt_features_gr_df.tsv"))
# print("L1_above4500nt_features_gr_df.tsv saved already")

################################################################################
# step3: Now we can bin the MACS2 values (-log10(p-value)) by regions 
################################################################################
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

# find_overlaps(FL_L1_anno_GR, L1binding_KZFPs_GR)
# FL_L1_bindingKZFP_GR <-
#   find_overlaps(FL_L1_anno_GR, L1binding_KZFPs_top20_GR) #273/288

# filter for PMM activated L1s only ############################################


# ---- Helper: Bin and summarize KZFP signal per feature ----
bin_feature_signal <- function(feature_gr, KZFP_gr, nbins, region_name, KZFP_name) {
  message("Binning ", region_name, " ", KZFP_name, " MACS2 into ",
          nbins, " bins...")
  
  # Tile into bins
  feature_binned <- tile(feature_gr, n = nbins) %>% unlist()
  
  # Assign bin metadata
  feature_binned$bin_id <- rep(seq_len(nbins), times = length(feature_gr))
  feature_binned$ID <- rep(feature_gr$ID, each = nbins)
  feature_binned$region <- region_name
  
  # Overlap with KZFP signal
  result <- feature_binned %>%
    join_overlap_left(KZFP_gr) %>%
    as_tibble() %>%  # Convert GRanges to data frame (with all metadata)
    group_by(ID, region, bin_id) %>%
    summarise(avg_score = mean(score, na.rm = TRUE), .groups = "drop") %>% 
    mutate(KZFP = KZFP_name)
  return(result)
}

# ---- Parameters ----
# features_list <- list(utr5_gr, orf1_gr, linkerToEnd_gr)
feature_names <- c("5UTR", "ORF1", "afterORF1")
bin_counts <- c(50, 25, 100)  # <-- our custom bin counts

utr5_gr_fil <- utr5_gr[width(utr5_gr)>=50]
orf1_gr_fil <- orf1_gr[width(orf1_gr)>=25]
linkerToEnd_gr_fil <- linkerToEnd_gr[width(linkerToEnd_gr)>=100]
features_list <- list(utr5_gr_fil, orf1_gr_fil, linkerToEnd_gr_fil)

# KZFP list
KZFP_list <- split(L1binding_KZFPs_top20_GR, 
                   L1binding_KZFPs_top20_GR$KZFP_name)
# features needed for function
feature_df <- tibble(
  gr = features_list,
  nbins = bin_counts,
  region_name = feature_names
)

# main loop, each KZFP feature bin ###############################################
# bin_results <- map_dfr(names(KZFP_list), function(KZFP_name) {
#   KZFP_gr <- KZFP_list[[KZFP_name]]
# 
#   pmap_dfr(feature_df, function(gr, nbins, region_name) {
#     bin_feature_signal(gr, KZFP_gr, nbins, region_name, KZFP_name)
#   })
# })

# ---- Save output ----
# write_tsv(bin_results, paste0(outdir, "/FL_L1PA2_binned_KZFP_score_293T_PMM.tsv"))
message("Binning complete: 
        saved as 'FL_L1PA2_binned_KZFP_score_293T_PMM.tsv'")

# ---- Plot mean signal across normalized bins ----
bin_results <- read_tsv(paste0(outdir, "/FL_L1PA2_binned_KZFP_score_293T_PMM.tsv"),
                        col_types = cols(
                          ID = col_integer(),
                          region = col_character(),
                          bin_id = col_integer(),
                          avg_score = col_double(),
                          KZFP = col_character()))

p1 <- bin_results %>% 
  dplyr::filter(KZFP %in% c("ZNF141", "ZNF382")) %>%
  dplyr::mutate(region = 
                  factor(region, levels = c("5UTR", "ORF1", "afterORF1"))) %>% 
  ggplot(aes(x = bin_id, y = avg_score, color = KZFP)) +
  geom_point()+
  geom_smooth()+
  # geom_line(size = 1.2) +
  facet_wrap(~region, scales = "free_x", 
             labeller = labeller(region = c(`5UTR` = "5'UTR", 
                                            ORF1 =  "ORF1",
                                            afterORF1 = "3'ORF1"))) +
  labs(
    # these score values are -log(pval) for each bp, binding signal vs input background
    title = "Average KZFP Score Across L1 Features",
    x = "Bin index",
    y = "Average KZFP score",
    color = "Region"
  ) +
  theme_pubr(base_size = 12)

ggsave(paste0(outdir, "/L1PA2_above4500nt_KZFP_binding_score_293T_PMM.pdf"),
       plot = p1,
       device = "pdf", width = 15, height = 10, units = "cm")

# Normalization bin_id → bin_prop ################################
df_clean <- bin_results %>%
  mutate(avg_score = ifelse(
    is.na(avg_score) | is.nan(avg_score) | is.infinite(avg_score), 0, avg_score))

## normalizing function for each KZFP
get_density_by_region <- function(df_region, region_name,
                                  adjust = 3,
                                  n = 1000,
                                  normalize = FALSE) {
  
  # total_signal <- sum(df_region$avg_signal)
  total_score <- sum(df_region$avg_score)
  x_range <- range(df_region$bin_id)
  x_seq <- seq(x_range[1], x_range[2], length.out = n)
  
  # no signal → flat 0 curve
  # if (total_signal == 0) {
  if (total_score == 0) {
    return(data.frame(x = x_seq, y = 0, region = region_name))
  }
  
  # weighted density
  dens <- density(
    x = df_region$bin_id,
    # weights = df_region$avg_signal / total_signal,
    weights = df_region$avg_score / total_score,
    from = x_range[1],
    to = x_range[2],
    n = n,
    adjust = adjust
  )
  
  # y_scaled <- dens$y * total_signal
  y_scaled <- dens$y * total_score
  
  if (normalize) {
    y_scaled <- y_scaled / max(y_scaled)
    # y_scaled <- y_scaled / max(avg_signal)
  }
  
  data.frame(x = dens$x, y = y_scaled, region = region_name)
}

# plot all L1 binding KZFP, PMM
density_all_KZFPs <- df_clean %>%
  # filter(KZFP %in% c("ZNF93", "ZNF649", "ZNF680", 
  #                    "ZNF141", "ZNF382", "ZNF425")) %>%
  group_by(KZFP, region) %>%
  group_split() %>%
  map_dfr(~ {
    df_density <- get_density_by_region(.x, unique(.x$region), normalize = FALSE)
    df_density$KZFP <- unique(.x$KZFP)
    df_density
  }) %>% 
  dplyr::mutate(region = factor(region, levels = c("5UTR", "ORF1", "afterORF1")))

p2 <- ggplot(density_all_KZFPs, aes(x = x, y = y, color = KZFP)) +
  geom_line()+
  facet_wrap(~region, scales = "free_x", 
             labeller = labeller(region = c(`5UTR` = "5'UTR", 
                                            ORF1 =  "ORF1",
                                            afterORF1 = "3'ORF1"))) +
  labs(
    x = "Bin index",
    y = "Weighted KZFP binding density",
    # title = "KLF6 binding across L1 sub-regions"
    title = "KZFP binding"
  ) +
  theme_pubr(base_size = 12)
ggsave(paste0(outdir, "/L1PA2_above4500nt_KZFP_binding_scoreDensity_293T_PMM.pdf"),
       plot = p2,
       device = "pdf", width = 15, height = 12, units = "cm")

# plot top L1 binding KZFP, PMM
density_top_KZFPs <- density_all_KZFPs %>%
  filter(KZFP %in% c("ZNF28", "ZNF680",
                     "ZNF141", "ZNF382", "ZNF425")) 

p3 <- ggplot(density_top_KZFPs, aes(x = x, y = y, color = KZFP)) +
  geom_line()+
  facet_wrap(~region, scales = "free_x", 
             labeller = labeller(region = c(`5UTR` = "5'UTR", 
                                            ORF1 =  "ORF1",
                                            afterORF1 = "3'ORF1"))) +
  labs(
    x = "Bin index",
    y = "Weighted KZFP binding density",
    title = "KZFP binding"
  ) +
  theme_pubr(base_size = 12)

ggsave(paste0(outdir, "/L1PA2_above4500nt_keyKZFP_binding_scoreDensity_293T_PMM.pdf"),
       plot = p3,
       device = "pdf", width = 15, height = 10, units = "cm")
