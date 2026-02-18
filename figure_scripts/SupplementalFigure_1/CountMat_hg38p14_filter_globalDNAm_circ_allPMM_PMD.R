library(tidyverse)
library(magrittr)
library(GenomicRanges)
library(plyranges)
library(circlize)
library(ComplexHeatmap) # Might be needed for color ramp, but circlize handles colors well

# Output directory
outDir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1"
if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load Data
# ------------------------------------------------------------------------------
message("Loading methylation matrix...")
# Path from implementation plan
meth_mat_path <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1/meth_snps_cov_fil_matGR_withXY.rds"
meth_snps_cov_fil_mat <- readRDS(meth_mat_path)

message("Loading PMD track...")
# Path from implementation plan
pmd_path <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/pmd_allPMM_samples/allPMM_samples_pmd.bed"
# Import BED: seqnames, start, end
# Filter for 1-22 only (User request: remove sex chromosomes)
# Trim to ensure ranges don't exceed chrom sizes (Fix warning)
# Need chrom sizes for trim(), we load them next step anyway.
# We can just filter seqnames here, trim later or ignore strictly if well formed.
# Better to trim after loading chrom sizes.
pmd_gr <- rtracklayer::import(pmd_path) %>%
    # Filter for 1-22 only (User request: remove sex chromosomes)
    plyranges::filter(seqnames %in% paste0("chr", 1:22))

# Reduce small gaps to make PMDs clearer (User request)
# Merge PMDs within 500kb (0.5MB) of each other
pmd_gr <- GenomicRanges::reduce(pmd_gr, min.gapwidth = 5e5)
# Filter out very small PMDs (<1MB) to reduce noise
pmd_gr <- pmd_gr[width(pmd_gr) > 1e6]

# ------------------------------------------------------------------------------
# 2. Process Methylation Data for Heatmap
# ------------------------------------------------------------------------------
# Extract mcols (methylation values) and GRanges info
m <- mcols(meth_snps_cov_fil_mat)
gr <- granges(meth_snps_cov_fil_mat)

# Select PMM samples
# --- Option 1: User defined order (Inner -> Outer) ---
# manual_order_nums <- c(11, 14, 18, 17, 13, 3, 12, 16, 15, 4, 6, 9, 1, 2, 7)
# manual_order <- paste0("PMM", manual_order_nums)
# # Intersect with actual columns to be safe
# available_cols <- grep("^PMM", colnames(m), value = TRUE)
# PMM_cols <- manual_order[manual_order %in% available_cols]
# message("User ordered samples: ", paste(PMM_cols, collapse = ", "))

# --- Option 2: Sort by Mean Methylation (High -> Low) ---
PMM_cols <- grep("^PMM", colnames(m), value = TRUE)
# Sort PMM samples by overall mean methylation (High -> Low)
mean_meth <- colMeans(as.matrix(m[, PMM_cols]), na.rm = TRUE)
PMM_cols <- names(sort(mean_meth, decreasing = TRUE))
message("Sorted samples (High -> Low): ", paste(PMM_cols, collapse = ", "))

message("Selected samples: ", paste(PMM_cols, collapse = ", "))

# Subset matrix
m_pmm <- as.data.frame(m[, PMM_cols])

# Combine with genomic coordinates
df_meth <- as.data.frame(gr) %>%
    dplyr::select(chr = seqnames, start, end) %>%
    bind_cols(m_pmm)

# Convert to numeric 0-1 (assuming original is 0-100 based on previous scripts, checking check:
# Ref script line 87: ifelse(PMM > 1, PMM/100, PMM).
# So likely 0-100. Let's convert to 0-1.
df_meth[, PMM_cols] <- lapply(df_meth[, PMM_cols], function(x) ifelse(x > 1, x / 100, x))

# Binning
# Circos heatmap works best with reduced data points.
# Let's bin the genome into windows (e.g. 5MB or 1MB)
message("Binning genome...")
# Create tiling window
# Read standard chrom sizes
chrom_sizes_df <- read.table("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/filtered_hg38_p14.chrom.sizes",
    col.names = c("chr", "size")
)
# Filter for 1-22
valid_chroms <- paste0("chr", 1:22)
chrom_sizes_df <- chrom_sizes_df %>% dplyr::filter(chr %in% valid_chroms)

# Now we can trim pmd_gr correctly if needed, or just rely on valid_chroms
seqlevels(pmd_gr, pruning.mode = "coarse") <- valid_chroms

# Create named vector
chrom_sizes <- setNames(chrom_sizes_df$size, chrom_sizes_df$chr)

# 1MB windows for higher resolution (User request)
window_size <- 1e6
tiles <- tileGenome(chrom_sizes, tilewidth = window_size, cut.last.tile.in.chrom = TRUE)

# Calculate average methylation in each bin
# Convert methylation df back to GRanges for overlap
meth_gr <- makeGRangesFromDataFrame(df_meth, keep.extra.columns = TRUE)

# Use binnedAverage (or similar manual aggregation)
# binnedAverage is slow for many columns. Let's use findOverlaps + data.table/dplyr aggregation
hits <- findOverlaps(tiles, meth_gr)

# Aggregate
# This might be memory intensive. Alternative:
# Iterate over samples?
# Or simple aggregate:
iso <- as.data.frame(mcols(meth_gr))
# iso$tile_id <- subjectHits(hits) # This assumes many meth sites map to one tile

# Calculate mean per tile
# Note: indices in 'queryHits' correspond to 'tiles'
# We need to group by queryHits(hits)
agg_res <- iso[subjectHits(hits), ] %>%
    dplyr::mutate(tile_idx = queryHits(hits)) %>%
    dplyr::group_by(tile_idx) %>%
    dplyr::summarise(across(all_of(PMM_cols), \(x) mean(x, na.rm = TRUE)))

# Map back to tiles
mcols(tiles) <- NULL # clear
mcols(tiles)[agg_res$tile_idx, PMM_cols] <- agg_res[, PMM_cols]

# Filter tiles to only those overlapping PMDs (User request - Reverted to show ALL data)
# tiles are GRanges, pmd_gr is GRanges
# tiles_pmd <- subsetByOverlaps(tiles, pmd_gr)

# Prepare data frame for circos
# Format: chr, start, end, value1, value2...
df_plot <- as.data.frame(tiles) %>%
    dplyr::select(chr = seqnames, start, end, all_of(PMM_cols)) %>%
    na.omit()

# ------------------------------------------------------------------------------
# 3. Circos Plot
# ------------------------------------------------------------------------------
out_pdf <- file.path(outDir, "circos_PMM_heatmap_PMD.pdf")
pdf(out_pdf, width = 8, height = 8)

circos.clear()
circos.par(start.degree = 90, gap.after = c(rep(1, length(unique(df_plot$chr)) - 1), 10))

# Initialize with genomic ranges (no ideogram bands)
# Ensure chrom_sizes order is correct
# Initialize with Ideogram (Standard hg38)
# This adds the outermost track with cytobands and labels
# Initialize with Ideogram (Standard hg38)
# plotType = c("ideogram", "axis") -> Suppresses default "chr1" labels
# plotType = c("ideogram") -> Suppresses default "chr1" labels AND the axis/ticks
# User request: "remove the ticks outside of ideogram, since the mb labels overlapped"
circos.initializeWithIdeogram(
    species = "hg38",
    chromosome.index = valid_chroms,
    plotType = c("ideogram")
)

# Add custom labels (Numbers only) to the Ideogram track (Track 1)
circos.track(track.index = 1, panel.fun = function(x, y) {
    chr_name <- CELL_META$sector.index
    # Show only number
    label <- gsub("chr", "", chr_name)
    # Plot label slightly above the ideogram
    # cex = 1.2 -> Larger than default, axis labels (0.9) are slightly smaller than this (User request)
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), label,
        facing = "inside", niceFacing = TRUE, adj = c(0.5, 0), cex = 1.2
    )
}, bg.border = NA)

# Track 1: PMD Regions (Outer)
# ------------------------------------------------------------------------------
# We can use circos.genomicRegion or just highlight them
# Let's verify PMD format: chr, start, end...
df_pmd <- as.data.frame(pmd_gr)[, c("seqnames", "start", "end")]
colnames(df_pmd)[1] <- "chr"

circos.genomicTrack(df_pmd,
    ylim = c(0, 1), track.height = 0.05,
    bg.border = NA, panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = "forestgreen", border = NA, ytop = 1, ybottom = 0)
    }
)
# Legend/Text for track? Or just visual
# Manually add track name?
# circos.text(...) inside panel?

# Track 2: PMM Heatmap (Inner)
# ------------------------------------------------------------------------------
# Color ramp: Blue (0) -> White (0.5) -> Red (1) (Standard methylation)
# Or Yellow -> Blue
# Color ramp: White (0) -> SkyBlue (0.5) -> Navy (1)
# Multi-stop gradient increases contrast in the middle range
col_fun <- colorRamp2(c(0, 0.5, 1), c("white", "deepskyblue", "navy"))

circos.genomicHeatmap(df_plot,
    col = col_fun, side = "inside", border = NA,
    heatmap_height = 0.7, connection_height = NULL
)

# Label the rows (Samples)
# Add text labels at the end of the last chromosome sector (chr22)
# Since gap.after is 10, there is space.
# We need to disable clipping or track overflow warning
# Label the rows (Samples) using circos.yaxis for correctness
# circos.yaxis draws labels at defined 'at' positions.
# For circos.genomicHeatmap, the Y-range is typically c(0, n_columns).
# So we need coordinates 0.5, 1.5, ..., N-0.5
circos.yaxis(
    at = 1:length(PMM_cols),
    labels = PMM_cols,
    sector.index = valid_chroms[length(valid_chroms)],
    track.index = get.current.track.index(),
    side = "right",
    labels.cex = 0.6,
    tick = FALSE # No distinct ticks needed, just labels
)

# Add Legend (using ComplexHeatmap logic separate from circos, strictly circos doesn't support legend easily inside)
# Usually we verify the visual output first.

# Add Legend
lgd_meth <- Legend(
    at = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"),
    col_fun = col_fun,
    title_position = "topcenter",
    title = "DNAm"
)

# Sample Legend (Inner to Outer)
# Reverse the order for legend so "Top" corresponds to "Outer" (Top of heatmap)?
# Default: Inner is bottom, Outer is top.
# Let's list from Outer (Top) to Inner (Bottom) to matches visually?
# Or just state "Inner -> Outer".
# Let's label "Samples\n(Outer -> Inner)" and reverse the list so the top item is the outer ring.
lgd_samples <- Legend(
    labels = rev(PMM_cols),
    title = "Samples\n(Outer -> Inner)",
    title_position = "topcenter",
    type = "points", pch = NA, background = NA # Pure text
)

# Pack legends and draw legends
# pd <- packLegend(lgd_meth, lgd_samples, direction = "vertical") # Vertical packing to stack them
pd <- packLegend(lgd_meth, direction = "vertical")
# Position at Top-Right
draw(pd, x = unit(0.95, "npc"), y = unit(0.95, "npc"))

dev.off()
message("Circos plot saved to: ", out_pdf)
