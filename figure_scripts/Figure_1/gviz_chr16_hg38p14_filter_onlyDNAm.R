
# filter out 3 B_EBV samples; NA values are grey in heatmap
# if set `window = -1` for DNAm dataTracks, it will make this script very slow
# png to pdf, res=1200
# add ABsignal tracks

.libPaths()
options(scipen=999)

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
print(seqlengths(seq_info)[chrom])


#### Create DNAm data tracks
outDir <- '/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1/'
meth_snps_fil_mat <- readRDS(paste0(outDir, "/meth_snps_fil_matGR_withXY.rds"))
meth_sub <- meth_snps_fil_mat[seqnames(meth_snps_fil_mat)%in%"chr16", 
                              !colnames(mcols(meth_snps_fil_mat)) %in% 
                                c("B1_EBV", "B2_EBV", "B3_EBV")]
bw_files <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/eRRBS_hg38_p14_bowtie1_filter/*/03*/*.bw")
meth_BMPC <- lapply(bw_files, import)
names(meth_BMPC) <- c("BMPC1", "BMPC2", "BMPC3")
# meth_BMPC: named list of GRanges, names are sample names (BMPC1/BMPC2/BMPC3)
meth_BMPC_renamed <- Map(function(gr, nm) {
  mcols(gr)$score <- mcols(gr)$score / 100
  colnames(mcols(gr)) <- nm
  gr
}, meth_BMPC, names(meth_BMPC))

meth_BMPC_gr <- Reduce(c, meth_BMPC_renamed)

meth_sub <- meth_sub %>% 
  plyranges::join_overlap_left(meth_BMPC_gr)

#Replace NA with 0 (optional)
# mcols(meth_sub)[] <- lapply(mcols(meth_sub), function(x) {
#   ifelse(is.na(x), 0, x)
# })
# Keep rows where at least one sample (metadata column) is NOT NA
# keep_rows <- rowSums(!is.na(as.data.frame(mcols(meth_sub)))) > 0
# meth_sub <- meth_sub[keep_rows]

# Your ranked indices (e.g., from some prior computation)
ranked_indices <- c("B1_rest", "B2_rest", "B3_rest", 
                    "BMPC1", "BMPC2", "BMPC3",
                    paste0("PMM", c(14, 18, 17, 13, 3, 
                                    12, 16, 15, 11, 4, 6, 
                                    9, 1, 2, 7)))
# Reorder the input vectors
mcols(meth_sub) <- mcols(meth_sub)[, ranked_indices]

# Create DataTrack for each BigWig file
data_tracks <- 
  DataTrack(range = meth_sub,
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
            background.title = "white",  # Ensure the title background is visible
            col.axis="black",
            # cex.axis= 0.33,
            cex.title = 9.6 / 12,
            # rotation.title = 0,
            # fontsize=7,
            # type = "smooth", #window = -1,
            type = c("heatmap"),
            # gradient = c("white", "blue"),  # Color gradient for scores
            # type = c("a","p"),
            # span=1/3,
            # window = -1,
            windowSize = 1000000,
            # windowSize = 10000000,
            showSampleNames = TRUE,
            fontfamily = "arial",
            fontface=1,
            fontface.title = 1,
            # showSampleNames = FALSE,
            col.sampleNames="black",
            cex.sampleNames = 9.6 / 12,
            showTitle = TRUE) # Ensure titles are shown

# Set NA color (grey or transparent)
displayPars(data_tracks)$na.color <- "#FFFFFF00" # for transparent


# tracks <- c(gTrack, data_tracks)

# #### Set up relative size for each track:
# sizes <- c(1.2, 4)  # Adjust size for data tracks


# Plot the tracks
dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1"

tiff(filename = paste0(dir, "/gviz_chr16_Heatmap_v4.tiff"),
     width = 12, height = 7.7, units = "cm", res = 1200)

plotTracks(data_tracks, chromosome = chrom, from = start, to = end,
           showTitle = FALSE # Ensure track titles are shown
           # lwd.border.title=1.5,
           # title.width = 1
           # sizes = sizes
           )
dev.off()


# pdf(file = paste0(dir, "/gviz_chr16_Heatmap_v7.pdf"),
#     width = 12 / 2.54, height = 7.7/ 2.54)
# plotTracks(data_tracks, chromosome = chrom, from = start, to = end,
#            showTitle = FALSE # Ensure track titles are shown
#            # lwd.border.title=1.5,
#            # title.width = 1
#            # sizes = sizes
# )
# dev.off()

