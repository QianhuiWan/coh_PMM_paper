library(tidyverse)
library(dplyr)
library(magrittr)
library(plyranges)
library(Gviz)
library(rtracklayer)
library(RColorBrewer)

# Gviz plot for PMM region
# load data
repeat_GR <- import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/hg38_rmsk.bed", format = "BED") %>%
  plyranges::filter(seqnames %in% c(paste0("chr", 1:22)))

cpgIsland_GR <- import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/cpgIsland/UCSC_CpGislands_unmask.bed") %>%
  plyranges::filter(seqnames %in% c(paste0("chr", 1:22)))

# get for bw files for PMM samples
genome <- "hg38"
chrom <- "chr5"
start <- 25377500
end <- 25392500

# add extra space to plot (since plot range should be longer than this DMR)
minbase <- start 
maxbase <- end

#### Ideogram Track:
iTrack <- IdeogramTrack(genome = genome, chromosome = chrom,  # Fix typo in chromosome name
                        from = minbase, to = maxbase,
                        fontsize = 9.6, fontcolor = "black", showTitle = TRUE)

#### Genome Axis Track:
gTrack <- GenomeAxisTrack(col = "black", 
                          cex = 0.5,
                          fontcolor = "black", add53 = TRUE, lwd = 0.5, 
                          # fontsize = 9.6, 
                          # labelPos = "below", 
                          name = "", exponent=6,
                          showTitle = TRUE, col.title = "black")


#### Repeat Track:
repeatTrack <- AnnotationTrack(range = repeat_GR, genome = genome, chromosome = chrom,
                               start = start, end = end,
                               name = "TE", title.width = 1,
                               id = repeat_GR$name,
                               fontcolor.feature = "black",
                               col.title = "black",
                               background.title = "white",  # Use white for title background
                               showTitle = TRUE,
                               showId = FALSE,
                               showFeatureId = FALSE,  # Ensure IDs are shown
                               cex.title = 0.5,
                               fontsize=9.6,
                               rotation.title = 0,
                               stacking = "dense")

#### CpG island track:
CpGislandTrack <- AnnotationTrack(range = cpgIsland_GR, genome = genome, chromosome = chrom,
                                  start = start, end = end,
                                  name = "CGI", title.width = 1,
                                  id = cpgIsland_GR$name,
                                  fontcolor.feature = "black",
                                  col.title = "black",
                                  background.title = "white",  # Use white for title background
                                  showTitle = TRUE,
                                  showFeatureId = FALSE,
                                  rotation.title = 0,
                                  cex.title = 0.5,
                                  fontsize=9.6,
                                  showId = TRUE,  # Show IDs for CpG island features
                                  stacking = "dense")

#### Create RNAseq data tracks:
file_paths <- list.files(path = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe", pattern = "*rpkm.bw", full.names = TRUE, recursive = TRUE)
# file_paths <- list.files(path = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe", pattern = "*cpm.bw", full.names = TRUE, recursive = TRUE)

# Your ranked indices (e.g., from some prior computation)
ranked_indices <- c("BMPC1", "BMPC2", "BMPC3", 
                    paste0("PMM", c(15, 4, 7, 14, 16, 3, 6, 11)))

# Reorder the input vectors
file_paths_meta <- file_paths %>% 
  as.data.frame() %>% 
  `colnames<-`("file_paths") %>% 
  dplyr::mutate(sample_sa = str_replace(basename(file_paths), 
                                        "_rpkm.bw", "")) %>% 
  dplyr::mutate(sample = str_replace_all(sample_sa, "_..+", "")) %>% 
  dplyr::mutate(names = str_replace(sample_sa, "_", "\n")) %>% 
  dplyr::arrange(match(sample, ranked_indices))

# data track in order:
data_tracks <- lapply(1:length(file_paths_meta$file_paths), function(i) {
  DataTrack(
    range = file_paths_meta$file_paths[i], genome = "hg38", 
    chromosome = chrom,
    name = file_paths_meta$names[i], 
    col.title = "black",
    background.title = "white",
    col.axis = "black",
    cex.title = 0.45,
    fontsize = 9,
    ylim = c(-100, 500),
    yTicksAt = c(0, 350),
    # rotation.title = 0,
    # ylim = c(-0.2, 1.2),
    # yTicksAt = c(0, 1),
    window = -1,
    type = "histogram",
    col.histogram = "darkblue",
    # fill.histogram = "lightblue",
    fill.histogram = "darkblue",
    showTitle = TRUE
  )
})


# lapply(data_tracks, function(track) {
#   displayPars(track) <- list(col = "darkblue", fill = "lightblue", showTitle = TRUE, type = "hist", window = -1)
# })

# Combine tracks
tracks <- c(iTrack, gTrack, repeatTrack, CpGislandTrack, data_tracks)

#### Set up relative size for each track:
sizes <- c(1.2, 3, 1.2, 1.2, rep(2, length(data_tracks)))  # Adjust size for data tracks

# Plot the tracks
dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/"

pdf(file = paste0(dir, "/gviz_antisense_chr5_example_rpkm.pdf"),
    width = 10 / 2.54, height = 13.5 / 2.54)
plotTracks(tracks,
           from = minbase - 1, to = maxbase + 1,
           showTitle = TRUE,  # Ensure track titles are shown
           # lwd.border.title=1.5, 
           title.width = 1,
           sizes = sizes)
dev.off()











