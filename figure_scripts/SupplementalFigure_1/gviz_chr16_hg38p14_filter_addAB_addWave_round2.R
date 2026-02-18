
# filter out 3 B_EBV samples; NA values are grey in heatmap
# if set `window = -1` for DNAm dataTracks, it will make this script very slow
# png to pdf, res=1200 dpi
# add ABsignal tracks
# R2: change A/B from hist to mountain plots; mountain plot not suitable 
# R2 not suitable -> failed

.libPaths()
options(scipen = 999)

library(tidyverse)
library(dplyr)
library(magrittr)
library(plyranges)
library(Gviz)
library(rtracklayer)
library(RColorBrewer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
seq_info <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

# Gviz plot for chr16
# get for bw files for PMM samples
genome <- "hg38"
chrom <- "chr16"
start <- 1
end <- 90338345

#### Ideogram Track:
iTrack <- IdeogramTrack(genome = genome, chromosome = chrom,  name = "chr6",
                        from = start, to = end, showTitle = FALSE,
                        fontsize = 9, fontcolor = "black")

#### Genome Axis Track:
gTrack <- GenomeAxisTrack(col = "black", cex = 1, name = "Genome Axis",
                          fontcolor = "black", add53 = TRUE, lwd = 0.5,
                          fontsize = 9, showTitle = TRUE, col.title = "black")

#### Create DNAm data tracks
outDir <- '/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1/'
meth_snps_fil_mat <- readRDS(paste0(outDir, "/meth_snps_fil_matGR_withXY.rds"))
meth_sub <- meth_snps_fil_mat[seqnames(meth_snps_fil_mat) %in% "chr16",
                              !colnames(mcols(meth_snps_fil_mat)) %in%
                                c("B1_EBV", "B2_EBV", "B3_EBV")]
# Replace NA with 0 (optional)
# mcols(meth_sub)[] <- lapply(mcols(meth_sub), function(x) {
#   ifelse(is.na(x), 0, x)
# })
# Keep rows where at least one sample (metadata column) is NOT NA
# keep_rows <- rowSums(!is.na(as.data.frame(mcols(meth_sub)))) > 0
# meth_sub <- meth_sub[keep_rows]

# Your ranked indices (e.g., from some prior computation)
ranked_indices <- c("B1_rest", "B2_rest", "B3_rest",
                    paste0("PMM", c(14, 18, 17, 13, 3,
                                    12, 16, 15, 11, 4, 6,
                                    9, 1, 2, 7)))
# Reorder the input vectors
mcols(meth_sub) <- mcols(meth_sub)[, ranked_indices]

# Create DataTrack for each BigWig file
data_tracks <- DataTrack(range = meth_sub,
                         genome = "hg38", chromosome = chrom,
                         start = start, end = end,
                         name = "DNAm", #DNA Methylation
            col.title = "black",
            # missingAsZero=TRUE,
            # lwd.border.title=1.5,
            ylim= c(0, 1),
            yTicksAt = c(0, 0.2, 0.4, 0.6, 0.8, 1),
            # groups= colnames(mcols(meth_snps_fil_mat[seqnames(meth_snps_fil_mat)%in%"chr16",])),
            # title.width = 1,
            background.title = "white", # Ensure the title background is visible
            col.axis = "black",
            # cex.axis= 0.33,
            # rotation.title = 0,
            # cex.title = 0.45,
            fontsize = 9,
            # type = "smooth", #window = -1,
            type = c("heatmap"),
            # gradient = c("white", "blue"),  # Color gradient for scores
            # type = c("a","p"),
            # span = 1/3,
            # window = -1,
            windowSize = 1000000,
            showSampleNames = TRUE,
            col.sampleNames = "black",
            # cex.sampleNames = 0.6,
            showTitle = TRUE)  # Ensure titles are shown

# Set NA color (grey or transparent)
displayPars(data_tracks)$na.color <- "#FFFFFF00" # for transparent

# Load GTF or GFF file using GenomicFeatures
gtf_data <- import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/gencode_v46_chr_patch_hapl_scaff_annotation.gtf")
mcols(gtf_data)$phase[mcols(gtf_data)$type == "stop_codon"] <- NA
export(gtf_data, "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/gencode_v46_chr_patch_hapl_scaff_annotation_phaseModified.gtf")

txdb <- makeTxDbFromGFF("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/gencode_v46_chr_patch_hapl_scaff_annotation_phaseModified.gtf", 
                        format = "gtf")

# Extract gene coordinates
genes <- genes(txdb)

new_seqinfo <- Seqinfo(seqnames = seqnames(seqinfo(BSgenome.Hsapiens.UCSC.hg38))[1:24],
                       seqlengths = seqlengths(seqinfo(BSgenome.Hsapiens.UCSC.hg38))[1:24],
                       isCircular = isCircular(seqinfo(BSgenome.Hsapiens.UCSC.hg38))[1:24],
                       genome = genome(seqinfo(BSgenome.Hsapiens.UCSC.hg38))[1:24])

valid_seqlevels <- intersect(seqlevels(genes) , seqlevels(new_seqinfo))

# Prune GRanges object to keep only valid seqlevels
genes_GR <- keepSeqlevels(genes, valid_seqlevels, pruning.mode = "coarse")

# Replace Seqinfo with the new Seqinfo object
seqinfo(genes_GR) <- new_seqinfo
strand(genes_GR) <- "*"

lcpm_df <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_star_allo_round4/counts_tx_round4/lcpm_TMM.tsv")
# deGenes <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_star_allo_round4/counts_tx_round4/topTable_res.tsv")

data <- import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/gencode_v46_chr_patch_hapl_scaff_annotation.gtf")
data <- as.data.frame(data, stringsAsFactors = FALSE)
data <- data[data$type == "transcript",]

lcpm_df_wide <-  lcpm_df %>%
  # dplyr::mutate(cpm = 2^lcpm) %>%
  pivot_wider(names_from = sample, values_from = lcpm,
              id_cols = gene_name)
lcpm_df_wide[-1] <- lapply(lcpm_df_wide[-1], function(x) ifelse(x < 0, 0, x))

library(dplyr)
genes_GR <- genes_GR %>% as.data.frame() %>%
  left_join(unique(data[,c("gene_id", "gene_name")]), by = "gene_id") %>%
  left_join(lcpm_df_wide, by = "gene_name") %>%
  mutate(
    # mean values: BMPC, PMM
    BMPC_lcpm_mean = rowMeans(dplyr::select(., BMPC1:BMPC3), na.rm = TRUE),
    PMM_lcpm_mean = rowMeans(dplyr::select(., PMM11:PMM7), na.rm = TRUE),
    lcpm_mean = rowMeans(dplyr::select(., BMPC1:PMM7), na.rm = TRUE),
    # sd values: BMPC, PMM
    BMPC_lcpm_sd = apply(dplyr::select(., BMPC1:BMPC3), 1, function(x) sd(x, na.rm = TRUE)),
    PMM_lcpm_sd = apply(dplyr::select(., PMM11:PMM7), 1, function(x) sd(x, na.rm = TRUE)),
    lcpm_sd = apply(dplyr::select(., BMPC1:PMM7), 1, function(x) sd(x, na.rm = TRUE))) %>% 
  replace_na(list(lcpm_mean = 0)) %>%
  replace_na(list(lcpm_sd = 0)) %>%
  replace_na(list(PMM_lcpm_sd = 0)) %>%
  replace_na(list(PMM_lcpm_mean = 0)) %>%
  dplyr::mutate(sd_mean = PMM_lcpm_sd/(PMM_lcpm_mean)) %>%
  replace_na(list(sd_mean = 0)) %>%
  as_granges()

mean_track <- DataTrack(
  range = genes_GR,
  data = genes_GR$PMM_lcpm_mean,
  genome = "hg38",  # Change to appropriate genome version
  chromosome = chrom, start = start, end = end,
  name = "Mean lcpm\n (PMM)", #Mean lcpm\n (PMM)
  # type = "smooth",
  span = 1/50,
  window = -1,
  # window = 10000,
  windowSize = 1000000,
  # ylim= c(0, 15), 
  ylim= c(-1, 8),
  yTicksAt = c(0, 2, 4, 6, 8),
  type = c("smooth","h"),
  # transformation = function(x) { x[x < 0] <- 0; x },
  col.title = "black",
  background.title = "white",
  col.axis="black",
  col = "blue")

variation_track <- DataTrack(
  range = genes_GR,
  data = genes_GR$sd_mean,
  chromosome = chrom, start = start, end = end,
  genome = "hg38",  # Change to appropriate genome version
  name = "SD/Mean lcpm\n (PMM)", #SD/Mean lcpm\n (PMM)
  type = c("smooth", "h"),
  span = 1/60,
  window = -1,
  windowSize = 1000000,
  # transformation = function(x) { x[x < 0] <- 0; x },
  # ylim= c(0, 5),
  ylim = c(-0.1, 1.5),
  yTicksAt = c(0, 0.5, 1),
  col.title = "black",
  background.title = "white",
  col.axis = "black",
  cex.axis = 0.65,
  col = "blue")

# ABana: S1R3: use 100000 bins (100kb bins) and use **PMM samples** only
absignal_allChr_GR_R3 <- read_tsv(paste0(outDir,
                                         "/absignal_100kbBins_PMMonly_allChr.tsv")) %>%
  as_granges()

# ABana: S1R4: use 100000 bins (100kb bins) and use **resting B cell samples** only
absignal_allChr_GR_R4 <- read_tsv(paste0(outDir,
                                         "/absignal_100kbBins_restBonly_allChr.tsv")) %>%
  as_granges()

# get PMD bed for all PMM sample
pmd_PMMs <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/pmd_allPMM_samples/allPMM_samples_pmd.bed")
pmd_i_GR <- import(pmd_PMMs) %>%
  plyranges::filter(seqnames %in%  c(paste0("chr", c(1:22, "X", "Y")))) %>%
  plyranges::mutate(pmd_id = paste0("pmd_", 1:length(.)))

#### PMD track
PMDTrack <- AnnotationTrack(range=pmd_i_GR, gensome = genome,
                            chromosome = chrom,
                            start = start, end = end,
                            name = "PMDs", #PMDs
                            id = "PMD",
                            fontcolor.feature = "black",
                            background.title = "transparent",
                            col.title = "black",
                            #groupAnnotation=NULL,
                            shape = "box",
                            showId = FALSE,
                            showFeatureId = FALSE,
                            stacking = "dense")


# Get the eigenvector values
pc1 <- absignal_allChr_GR_R3$pc

# Build two-group matrix: one for positive, one for negative
pc_matrix <- rbind(
  A = ifelse(pc1 > 0, pc1, 0),
  B = ifelse(pc1 < 0, -pc1, 0)  # Flip B so it's positive (for plotting)
)

mountain_track <- DataTrack(
  range = gr,
  data = pc_matrix,
  groups = c("A", "B"),
  type = "mountain",
  genome = "hg38",
  chromosome = chrom,
  start = start,
  end = end,
  ylim = c(0, 2),  # Adjust depending on scale of PC1
  name = "Compartment Score",
  col.mountain = c("firebrick3", "dodgerblue3"),  # A = red, B = blue
  fill.mountain = c("firebrick3", "dodgerblue3"),
  col.axis = "black",
  background.title = "white",
  col.title = "black",
  cex.axis = 0.6,
  cex.title = 0.65)

absignal_PMM_track <- DataTrack(
  range = absignal_allChr_GR_R3,
  data = absignal_allChr_GR_R3$pc,
  groups = c("open", "closed"),
  chromosome = chrom, start = start, end = end,
  genome = "hg38",  # Change to appropriate genome version
  name = "Eigenvector\n(PMM)", #Eigenvector\n(PMM)
  # type = c("hist"),
  type = "mountain",
  # span = 1/60,
  window = -1,
  windowSize = 1000000,
  transformation = function(x) { (-1) * x },
  # ylim= c(0, 5),
  ylim = c(-2, 2),
  yTicksAt = c(-1, 0, 1),
  col.title = "black",
  background.title = "white",
  col.axis = "black",
  cex.axis = 0.65,
  col.mountain = c("firebrick3", "dodgerblue3"),  # A = red, B = blue
  fill.mountain = c("firebrick3", "dodgerblue3"),
  # col.histogram = "grey50",
  # fill.histogram = "grey50"
  )

absignal_restB_track <- DataTrack(
  range = absignal_allChr_GR_R4,
  data = absignal_allChr_GR_R4$pc,
  chromosome = chrom, start = start, end = end,
  genome = "hg38",  # Change to appropriate genome version
  name = "Eigenvector\n(Resting B)", #Eigenvector\n(Resting B)
  type = c("histogram"),
  window = -1,
  windowSize = 1000000,
  transformation = function(x) { (-1) * x },
  ylim = c(-2, 2),
  yTicksAt = c(-1, 0, 1),
  col.title = "black",
  background.title = "white",
  col.axis = "black",
  cex.axis = 0.65,
  col.histogram = ifelse(absignal_allChr_GR_R4$pc < 0, "grey20", "grey50"),
  fill.histogram = ifelse(absignal_allChr_GR_R4$pc < 0, "grey20", "grey50"))

# Wave score
HEPG2_Wave_hg38GR <- read_bed_graph(paste0(outDir, "/waveScore_hg38/HEPG2_Wave_hg38.bedgraph"))

HEPG2_Wave_track <- DataTrack(
  range = HEPG2_Wave_hg38GR,
  data = HEPG2_Wave_hg38GR$score,
  chromosome = chrom, start = start, end = end,
  genome = "hg38",  # Change to appropriate genome version
  name = "Wave score\n(HEPG2)", #Wave score\n(HEPG2)
  type = c("hist"),
  window = -1,
  windowSize = 1000000,
  # transformation = function(x) { (-1) * x },
  ylim = c(-10, 100),
  yTicks = c(0, 40, 80),
  col.title = "black",
  background.title = "white",
  col.axis = "black",
  cex.axis = 0.65,
  # col = ifelse(HEPG2_Wave_hg38GR$score < 0, "deeppink4", "grey50"),
  # col = "grey20",
  col = ifelse(HEPG2_Wave_hg38GR$score < 0, "grey20", "grey50"),
  fill = ifelse(HEPG2_Wave_hg38GR$score < 0, "grey20", "grey50"))

# Combine tracks
# tracks <- c(iTrack, gTrack, mean_track,
#             variation_track, data_tracks)

tracks <- c(iTrack, gTrack, mean_track,
            variation_track,
            # data_tracks,
            PMDTrack,
            absignal_PMM_track, absignal_restB_track,
            HEPG2_Wave_track)

#### Set up relative size for each track:
sizes <- c(1.2, 2, 3, 3,
           # 4,
           1.2,
           1.5, 1.5,
           1.5)  # Adjust size for data tracks

# Plot the tracks
dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1"

tiff(filename = paste0(dir, "/gviz_chr16_addAB_addWave_addPMD_rmNames.tiff"),
     width = 12, height = 9, units = "cm", res = 1200)

plotTracks(tracks, chromosome = chrom, from = start, to = end,
           showTitle = FALSE,  # Ensure track titles are shown
           # lwd.border.title=1.5,
           # title.width = 1,
           sizes = sizes
           )
dev.off()


pdf(file = paste0(dir, "/gviz_chr16_addAB_addWave_addPMD_rmNames.pdf"),
    width = 12 / 2.54, height = 9 / 2.54)
plotTracks(tracks, chromosome = chrom, from = start, to = end,
           showTitle = FALSE,  # Ensure track titles are shown
           # lwd.border.title=1.5,
           title.width = 1,
           sizes = sizes)
dev.off()


# cairo_pdf(paste0(dir, "/gviz_chr16_addHeatmap_addAB_addWave_addPMD_rmNames.pdf"),
#           width = 10 / 2.54, height = 8 / 2.54)
# plotTracks(tracks, chromosome = chrom, from = start, to = end,
#            showTitle = TRUE,  # Ensure track titles are shown
#            # lwd.border.title=1.5,
#            title.width = 1,
#            sizes = sizes)
# dev.off()

# save as svg
# svg(paste0(dir, "/gviz_chr16_addHeatmap_addAB_addWave_addPMD_rmNames.svg"),
#     width = 12 / 2.54, height = 9 / 2.54)  # Adjust width and height as needed
# plotTracks(tracks, chromosome = chrom, from = start, to = end,
#            showTitle = TRUE,  # Ensure track titles are shown
#            title.width = 1,
#            sizes = sizes)
# dev.off()

# png to pdf
library(png)
library(grid)
library(magick)
# List all PNG files you want to combine
# grDevices::pdf(file = paste0(dir, "/gviz_chr16_addHeatmap_addAB_addWave_addPMD_TIFF.pdf"),
#                width = 12, height = 10)
# png_file <- paste0(dir, "/gviz_chr16_addHeatmap_addAB_addWave_addPMD.tiff")
# img <- readPNG(png_file)        # Read the PNG image
# readTIFF
# grid.raster(img) 
# grDevices::dev.off()

# img <- image_read(paste0(dir, "/gviz_chr16_addHeatmap_addAB_addWave_addPMD.tiff"))
# image_write(img, paste0(dir, "/gviz_chr16_addHeatmap_addAB_addWave_addPMD_TIFF.pdf"),
#             format = "pdf")

