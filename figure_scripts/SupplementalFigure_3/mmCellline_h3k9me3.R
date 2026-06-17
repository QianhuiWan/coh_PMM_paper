#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(plyranges)
  library(rtracklayer)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# =========================
# 用户需要改的参数
# =========================
# promoter_bed <- "L1_activated_promoters.bed"   # 你的 activated L1 promoter 区域（BED）
activeL1_sense <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/activationScore/filtered_L1/sense_L1s.tsv")
activeL1_antisense <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe_uniqueMap/activationScore/filtered_L1/antisense_L1s.tsv")
all_activeL1 <- rbind(activeL1_sense, activeL1_antisense)
all_activeL1$active <- rep("Yes", nrow(all_activeL1))

map <- c(
  "CD19-GSM2171212" = "BMPC4",
  "CD19-GSM2171213" = "BMPC5",
  "CD19-GSM2171214" = "BMPC6",
  "CD19-GSM2171215" = "BMPC7"
)
all_activeL1$sampleID <- str_replace_all(all_activeL1$sampleID, map)

all_activeL1_PMM_pro <- all_activeL1 %>% 
  dplyr::filter(str_detect(sampleID, "PMM")) %>% 
  dplyr::select(te_id = Geneid, repName, repClass, repFamily, 
                seqnames = genoName, start = genoStart, end = genoEnd,
                strand) %>% 
  unique() %>% 
  as_granges() %>% 
  promoters(upstream = 1000, downstream = 1000)
  
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
repeat_GR_inact_L1_pro <- 
  repeat_bed_df %>%
  dplyr::mutate(seqnames = genoName, 
                start = genoStart + 1, end = genoEnd) %>%
  dplyr::select(-c(genoName:genoEnd)) %>% 
  unique() %>% 
  as_granges() %>% 
  plyranges::filter(repFamily=="L1") %>% 
  plyranges::filter(width(.) > 5900) %>% 
  plyranges::filter(!(te_id %in% all_activeL1_PMM_pro$te_id)) %>% 
  promoters(upstream = 1000, downstream = 1000)

bw_dir       <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/track_hub_settings/histone_hubDir/histone"     
meta    <- tibble( sample = c("MM1S", "KMS11"), 
                        group = c("MM_cell_line", "MM_cell_line"))      # 样本信息：sample, group（至少两列）

# bigWig 文件名如何对应 sample：
# 默认规则：bigWig 文件 basename 去掉 .bw/.bigWig 后 = sample
# 例：PMM01.bw -> sample = "PMM01"
bw_pattern <- "\\.(bw|bigWig)$"

# 你想用哪种 promoter summary：
# "wmean" = 加权平均(推荐，考虑每段 width)
# "mean"  = 简单平均(不考虑 width)
# "sum"   = score * width 累加（类似 area）
# "max"   = 最大值
# summary_method <- "wmean"
summary_method <- "mean"

# 遇到某 promoter 在某 bigWig 没信号时，用 0 还是 NA：
fill_no_signal <- 0  # 改成 NA 也行
# 是否对结果做 log1p（画图更稳）
do_log1p <- FALSE
# do_log1p <- TRUE

# =========================
# 找到 bigWig 文件
# =========================
bw_files <- list.files(bw_dir, pattern = bw_pattern, full.names = TRUE)
if (length(bw_files) == 0) stop("No bigWig files found in: ", bw_dir)

bw_tbl <- data.frame(
  bw_path = bw_files,
  bw_base = basename(bw_files),
  sample  = str_replace(basename(bw_files), bw_pattern, ""),
  stringsAsFactors = FALSE
) %>% 
  dplyr::mutate( sample  = ifelse(str_detect(sample, "KMS"), "KMS11", "MM1S"))

# 保留 metadata 里存在的样本
bw_tbl <- bw_tbl %>% semi_join(meta, by = "sample")
if (nrow(bw_tbl) == 0) stop("No bigWig samples match your metadata sample column.")

message("Found ", nrow(bw_tbl), " bigWig files matching metadata.")


# =========================
# 核心函数：单个 bigWig 在 promoters 上算 summary
# =========================
calc_bw_on_regions <- function(bw_path, regions, method = c("wmean","mean","sum","max"),
                               fill_no_signal = 0) {
  method <- match.arg(method)
  
  # 只导入 regions 覆盖到的片段，速度会快很多
  bw <- BigWigFile(bw_path)
  
  # import 会返回 bigWig 里有信号的区段（GRanges，score列）
  sig <- import(bw, which = regions)
  
  # 默认全都没信号
  out <- rep(fill_no_signal, length(regions))
  names(out) <- mcols(regions)$te_id
  
  if (length(sig) == 0) return(out)
  
  # 计算每个 promoter 的 summary
  hits <- findOverlaps(regions, sig, ignore.strand = TRUE)
  if (length(hits) == 0) return(out)
  
  qh <- queryHits(hits)
  sh <- subjectHits(hits)
  
  # promoter 与 signal 的交集区段
  inter <- pintersect(ranges(regions[qh]), ranges(sig[sh]))
  inter_w <- width(inter)
  sc <- mcols(sig)$score[sh]
  
  dt <- data.table(
    promoter_idx = qh,
    w = inter_w,
    score = sc
  )
  
  if (method == "max") {
    val <- dt[, .(v = max(score, na.rm = TRUE)), by = promoter_idx]
  } else if (method == "sum") {
    # 类似 area：score * width 累加
    val <- dt[, .(v = sum(score * w, na.rm = TRUE)), by = promoter_idx]
  } else if (method == "mean") {
    val <- dt[, .(v = mean(score, na.rm = TRUE)), by = promoter_idx]
  } else if (method == "wmean") {
    # 加权平均：sum(score*width)/sum(width)
    val <- dt[, .(v = sum(score * w, na.rm = TRUE) / sum(w, na.rm = TRUE)), by = promoter_idx]
  }
  
  out[val$promoter_idx] <- val$v
  names(out) <- mcols(regions)$te_id
  out
}

# =========================
# 批量计算：每个样本一列（promoter × sample matrix）
# =========================
message("Computing signal summaries (method = ", summary_method, ") ...")

mat_list_activeL1 <- lapply(seq_len(nrow(bw_tbl)), function(i) {
  s <- bw_tbl$sample[i]
  p <- bw_tbl$bw_path[i]
  message("  - ", s, " : ", basename(p))
  v <- calc_bw_on_regions(p, all_activeL1_PMM_pro, method = summary_method, fill_no_signal = fill_no_signal)
  data.frame(promoter = names(v), sample = s, value = as.numeric(v), stringsAsFactors = FALSE)
})

df_long_activeL1 <- bind_rows(mat_list_activeL1) %>%
  left_join(meta, by = "sample")
df_long_activeL1$active <- rep("Yes", nrow(df_long_activeL1))

mat_list_inact_L1 <- lapply(seq_len(nrow(bw_tbl)), function(i) {
  s <- bw_tbl$sample[i]
  p <- bw_tbl$bw_path[i]
  message("  - ", s, " : ", basename(p))
  v <- calc_bw_on_regions(p, repeat_GR_inact_L1_pro, method = summary_method, fill_no_signal = fill_no_signal)
  data.frame(promoter = names(v), sample = s, value = as.numeric(v), stringsAsFactors = FALSE)
})

df_long_inact_L1 <- bind_rows(mat_list_inact_L1) %>%
  left_join(meta, by = "sample")
df_long_inact_L1$active <- rep("No", nrow(df_long_inact_L1))

df_long <- rbind(unique(df_long_activeL1), unique(df_long_inact_L1))

if (do_log1p) {
  df_long <- df_long %>% mutate(value = log1p(value))
}

# =========================
# 画图：group 之间 boxplot + jitter
# =========================
p <- ggplot(df_long, aes(x = sample, y = value, fill = active)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_boxplot() +
  # geom_jitter(width = 0.2, size = 0.8, alpha = 0.6) +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    label = "p.signif", vjust = 1
  ) +
  # facet_wrap(~sample, scales = "free_x") +
  ggpubr::theme_pubr(base_size = 9.6)+
  ylim(c(0,8))+
  labs(
    x = "",
    y = ifelse(do_log1p, "log1p(promoter H3K9me3 signal)", 
               paste0("L1 promoter H3K9me3 signal (", summary_method, ")")), 
    fill = "Transcriptional active L1s"
  )

print(p)
ggsave(paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/chipSeq_histone_MMcellline",
              "/l1_promoter_h3k9me3.pdf"),
       plot = p,
       width = 10, height = 8, units = "cm", limitsize = FALSE)

# =========================
# 输出结果表
# =========================
out_tsv <- paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/chipSeq_histone_MMcellline/",
                  "fulllength_L1promoter_signal_", summary_method, ifelse(do_log1p, "_log1p", ""), ".tsv")
fwrite(df_long, out_tsv, sep = "\t")
message("Done. Wrote: ", out_tsv)

