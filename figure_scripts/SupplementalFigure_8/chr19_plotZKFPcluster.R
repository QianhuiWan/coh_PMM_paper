

# library(karyoploteR)
# markers <- data.frame(chr=rep("chr19", 10), pos=(1:10*10e6), 
#                       labels=paste0("Gene", 1:10))
# kp <- plotKaryotype(chromosomes="chr19")
# kpAddBaseNumbers(kp)
# kpPlotMarkers(kp, chr=markers$chr, x=markers$pos, labels=markers$labels)


# Load required libraries
library(karyoploteR)
library(tidyverse)
library(GenomicRanges)


# ------------------------------
# 1. Load your gene coordinates
# ------------------------------

kzfp_genes <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/kzfps_hg38_all.tsv")
kzfp_genes_GR <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/kzfps_hg38_all_geneCoordinates.tsv") %>% 
  as_granges()

# ------------------------------
# 2. Load your cluster GRanges
# ------------------------------

KZFP_cluster_GR <- 
  read_tsv(paste0(KZFP_dir, "/KZFP_clusters_coord_hg38.tsv")) %>% 
  as_granges()

# ------------------------------
# 3. Plotting
# ------------------------------

# Make sure seqnames are in UCSC format ("chr19")
library(GenomeInfoDb)
seqlevelsStyle(kzfp_genes_GR) <- "UCSC"
seqlevelsStyle(KZFP_cluster_GR) <- "UCSC"
# Filter clusters to only those on chr19
KZFP_cluster_chr19 <- KZFP_cluster_GR[seqnames(KZFP_cluster_GR) == "chr19"] %>% 
  GenomicRanges::reduce(ignore.strand=TRUE) %>% 
  plyranges::mutate(clusterID = paste0("C", 1:length(.)))

kzfp_genes_GR_chr19 <- kzfp_genes_GR %>% plyranges::filter(seqnames == "chr19")%>% 
  GenomicRanges::sort(ignore.strand=TRUE)
  
kzfp_genes_GR_chr19_C6 <- kzfp_genes_GR %>% 
  filter_by_overlaps(KZFP_cluster_chr19[KZFP_cluster_chr19$clusterID=="C6",])%>% 
  GenomicRanges::sort(ignore.strand=TRUE)



pdf(paste0(KZFP_dir, "/kzfp_chr19_cluster6.pdf"), width = 12, height = 4)
# Initialize plot
kp <- plotKaryotype(chromosomes = "chr19")
kpAddBaseNumbers(kp)

# Plot clusters as blocks below chromosome
kpPlotRegions(
  kp,
  data = KZFP_cluster_chr19,
  col = "#FFC30066",     # transparent yellow fill
  border = "#FFC30066",    
  r0 = 0.08,
  r1 = -0.09,
  # r0 = 0.1, r1 = 0.3,
  # r0 = 0.01, r1 = 0.1,
  data.panel = 1        # default
)

# Label cluster regions
# Label cluster regions INSIDE the yellow blocks
for (i in seq_along(KZFP_cluster_chr19)) {
  cluster <- KZFP_cluster_chr19[i]
  
  if (length(cluster) > 0 && !is.na(start(cluster))) {
    kpText(
      kp,
      chr = as.character(seqnames(cluster)),
      x = (start(cluster) + end(cluster)) / 2,   # center of block
      y = 0,                                   # middle between r0 = 0.1 and r1 = 0.3
      labels = paste0("C", i),
      cex = 0.8,
      col = "black",                             # good contrast with yellow
      font = 2                                   # optional: bold
    )
  }
}


# ---- Spread gene labels horizontally across a defined range ----

# Total number of genes
n_genes <- length(kzfp_genes_GR_chr19_C6)

# Define the horizontal span for label placement (e.g., the cluster region)
x_min <- min(start(kzfp_genes_GR_chr19))
x_max <- max(end(kzfp_genes_GR_chr19_C6))
label_xs <- seq(x_min, x_max, length.out = n_genes)

# Set fixed y position for labels
label_y <- 0.25
# Define colors by strand
strand_colors <- ifelse(strand(kzfp_genes_GR_chr19_C6) == "+","#d62728", "#1f77b4")

# Plot gene names with lines to true start positions
for (i in seq_along(kzfp_genes_GR_chr19_C6)) {
  gene <- kzfp_genes_GR_chr19_C6[i]
  label <- mcols(gene)$gene_name
  
  # Draw text label (vertical)
  kpText(
    kp,
    chr = as.character(seqnames(gene)),
    x = label_xs[i],    # evenly spaced across region
    y = label_y,
    labels = label,
    srt = 90,
    cex = 1,
    pos = 3,
    adj = 2,   # left-aligned
    col = strand_colors[i], 
    data.panel = 1
  )
  
  # Draw connector line from gene start to label
  kpSegments(
    kp,
    chr = as.character(seqnames(gene)),
    x0 = start(gene),
    x1 = label_xs[i],
    y0 = 0.08,       # top of yellow block
    y1 = label_y,
    col = "#999999",
    lwd = 1,
    data.panel = 1
  )
}

dev.off()








