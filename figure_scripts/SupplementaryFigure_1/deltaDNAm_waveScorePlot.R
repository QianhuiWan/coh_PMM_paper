
library(tidyverse)
library(plyranges)
library(ggplot2)
library(hexbin)
outDir <- '/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1'


# load hg18 liftover hg38 wave bedgraph file ###################################
HEPG2_Wave_hg38GR <- read_bed_graph(paste0(outDir, "/waveScore_hg38/HEPG2_Wave_hg38.bedgraph"))

# delta DNAm data ##############################################################
## load DNAm for all PMM and controls
meth_snps_fil_matGR <- readRDS(paste0(outDir, "/meth_snps_fil_matGR_withXY.rds"))

## calculate delta DNAm: PMM-restB
# Mean of columns ending with "rest"
mean_rest <- meth_snps_fil_matGR %>%
  select(ends_with("rest")) %>%
  mcols() %>%
  as.data.frame() %>%
  rowMeans(na.rm = TRUE)

# Mean of columns starting with "PMM"
mean_pmm <- meth_snps_fil_matGR %>%
  select(starts_with("PMM")) %>%
  mcols() %>%
  as.data.frame() %>%
  rowMeans(na.rm = TRUE)

# Add the results back to the GRanges object
meth_snps_fil_matGR_deltaDNAm <- 
  meth_snps_fil_matGR[,0] %>%
  mutate(mean_rest = mean_rest, 
         mean_pmm = mean_pmm,
         deltaDNAm = mean_pmm - mean_rest) %>% 
  filter(!is.nan(deltaDNAm))


# plot delta DNAm against wave score
meth_snps_fil_matGR_deltaDNAm_waveS <- 
  find_overlaps(meth_snps_fil_matGR_deltaDNAm, HEPG2_Wave_hg38GR)


# plot, use ggplot
# plot <- meth_snps_fil_matGR_deltaDNAm_waveS %>%
#   as.data.frame() %>%
#   ggplot(aes(x = score, y = deltaDNAm)) +
#   geom_point(alpha = 0.2, size = 0.5) +  # Transparent and smaller points
#   geom_smooth(method = "loess", se = FALSE, color = "blue") +  # Add smooth trendline
#   labs(x = "Wave Score (HEPG2)", y = "Delta DNAm", 
#        title = "Trend of Delta DNAm by Wave Score") +
#   theme_minimal()
# 
# ggsave(paste0(outDir, "/deltaDNAm_againstHEPG2_Wave_hg38.pdf"),
#        plot = plot,
#        device = "pdf", width = 22, height = 15, units = "cm")

# use hexbin
# library(hexbin)
# plot <- meth_snps_fil_matGR_deltaDNAm_waveS %>%
#   as.data.frame() %>%
#   ggplot(aes(x = score, y = deltaDNAm)) +
#   geom_hex(bins = 30) +  # Create hexagonal bins
#   geom_smooth(method = "loess", se = FALSE, color = "blue") +  # Add smooth trendline
#   labs(x = "Wave Score", y = "Delta DNAm", title = "Hexbin Trend of Delta DNAm by Wave Score") +
#   theme_minimal()
# 
# ggsave(paste0(outDir, "/deltaDNAm_againstHEPG2_Wave_hg38_hexbin.pdf"),
#        plot = plot,
#        device = "pdf", width = 22, height = 15, units = "cm")



# plot mean + SD
# plot <- meth_snps_fil_matGR_deltaDNAm_waveS %>%
#   as.data.frame() %>%
#   dplyr::mutate(wave_score = factor(score, levels = unique(sort(score)))) %>%
#   dplyr::group_by(wave_score) %>%
#   dplyr::summarize(
#     mean_deltaDNAm = mean(deltaDNAm, na.rm = TRUE),
#     sd_deltaDNAm = sd(deltaDNAm, na.rm = TRUE)
#   ) %>%
#   ggplot(aes(x = wave_score, y = mean_deltaDNAm)) +
#   geom_point(color = "red", size = 2) +  # Plot aggregated mean points
#   geom_errorbar(aes(ymin = mean_deltaDNAm - sd_deltaDNAm, 
#                     ymax = mean_deltaDNAm + sd_deltaDNAm), 
#                 width = 0.2, color = "black") +  # Add error bars for SD
#   geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "blue") +  # Add trendline
#   labs(x = "Wave Score (HEPG2)", y = "Delta DNAm", 
#        title = "Trend of Delta DNAm by Wave Score") +
#   theme_minimal()
# 
# ggsave(paste0(outDir, "/deltaDNAm_againstHEPG2_Wave_hg38_meanAndSD.pdf"),
#        plot = plot,
#        device = "pdf", width = 22, height = 15, units = "cm")



# plot in 10 bins for x axis
# plot_in10Bins <- 
#   meth_snps_fil_matGR_deltaDNAm_waveS %>%
#   as.data.frame() %>%
#   dplyr::mutate(
#     wave_bins = cut(score, breaks = 10)  # Bin score into 10 equal intervals
#   ) %>%
#   dplyr::filter(!is.na(deltaDNAm)) %>%   # Remove NA values for deltaDNAm
#   ggplot(aes(x = wave_bins, y = deltaDNAm)) +
#   geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +  # Reduce outlier size and transparency
#   geom_hline(yintercept = 0, linetype = "dashed", color = "blue", size = 0.8) +  # Add dashed line at y=0
#   labs(x = "Score Bins", y = "deltaDNAm") +               # Label axes
#   theme_minimal(base_size = 10)                           # Minimal theme for better readability
# 
# 
# ggsave(paste0(outDir, "/deltaDNAm_againstHEPG2_Wave_hg38_x10Bins.pdf"),
#        plot = plot_in10Bins, 
#        device = "pdf", width = 22, height = 15, units = "cm")


# plot in x bins with width of 1
plot_width1_xBins <- 
  meth_snps_fil_matGR_deltaDNAm_waveS %>%
  as.data.frame() %>%
  ggplot(aes(x = cut_width(score, width = 1, boundary = 0), y = deltaDNAm)) +
  geom_boxplot(
    # outlier.size = 0.5, outlier.alpha = 0.5, 
    outlier.shape = NA) +  # Adjust outlier size and transparency
  stat_summary(fun = "mean", geom = "point", color = "red", size = 1.5) +  # Add mean points
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue", size = 0.8) +  # Add dashed line at y=0
  labs(
    x = "Wave score Bins (HEPG2)", 
    y = "delta DNAm\n(PMM-restingB)", 
    title = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))  # Rotate x-axis labels

ggsave(paste0(outDir, "/deltaDNAm_againstHEPG2_Wave_hg38_width1_xBins_rmNA.pdf"),
       plot = plot_width1_xBins, 
       device = "pdf", width = 35, height = 15, units = "cm")


# plot in x bins with width of 0.5
# plot_width0_5_xBins <- 
#   meth_snps_fil_matGR_deltaDNAm_waveS %>%
#   as.data.frame() %>%
#   ggplot(aes(x = cut_width(score, width = 0.5, boundary = 0), y = deltaDNAm)) +
#   geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +  # Adjust outlier size and transparency
#   stat_summary(fun = "mean", geom = "point", color = "red", size = 1.5) +  # Add mean points
#   geom_hline(yintercept = 0, linetype = "dashed", color = "blue", size = 0.8) +  # Add dashed line at y=0
#   labs(
#     x = "Wave score Bins (HEPG2)", 
#     y = "delta DNAm\n(PMM-restingB)", 
#     title = ""
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))  # Rotate x-axis labels
# 
# 
# ggsave(paste0(outDir, "/deltaDNAm_againstHEPG2_Wave_hg38_width0_5_xBins.pdf"),
#        plot = plot_width0_5_xBins, 
#        device = "pdf", width = 35, height = 15, units = "cm")


# plot in x bins with width of 0.1
# plot_width0_1_xBins <- 
#   meth_snps_fil_matGR_deltaDNAm_waveS %>%
#   as.data.frame() %>%
#   ggplot(aes(x = cut_width(score, width = 0.1, boundary = 0), y = deltaDNAm)) +
#   geom_boxplot(outlier.size = 0.5, 
#                outlier.alpha = 0.5, 
#                outlier.shape = NA) +  # Adjust outlier size and transparency
#   stat_summary(fun = "mean", geom = "point", color = "red", size = 1.5) +  # Add mean points
#   geom_hline(yintercept = 0, linetype = "dashed", color = "blue", size = 0.8) +  # Add dashed line at y=0
#   labs(
#     x = "Wave score Bins (HEPG2)", 
#     y = "delta DNAm\n(PMM-restingB)", 
#     title = ""
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))  # Rotate x-axis labels
# 
# 
# ggsave(paste0(outDir, "/deltaDNAm_againstHEPG2_Wave_hg38_width0_1_xBins.pdf"),
#        plot = plot_width0_1_xBins, 
#        device = "pdf", width = 35, height = 15, units = "cm")

