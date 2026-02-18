

# in this round3, change to methrix R packge
library(magrittr)
library(methrix)
library(stringr)
library(plyranges)
library(rtracklayer)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MafDb.1Kgenomes.phase3.GRCh38)

# read bed.graph files
# files <- list.files("/Users/qianhuiwan/gitHubRepo/coh_PMM/data/testRes_DNAm/reformatedCpGmeth", pattern = "formatted_perc.txt", full.names = TRUE) 
files <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/*/03.DNMtoolsPosMethylation/*_CpG.meth.gz")

print(basename(files))

outDir <- '/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1'

# get sample id and group
sampleID <- str_replace_all(basename(files), "_CpG.meth.gz", "")
group <- str_replace_all(sampleID, "[:digit:]", "")
# group <- ifelse(group == "B_rest", 0, ifelse(group == "PMM", 1, 2))

# meth_snps_fil_mat <- readRDS(paste0(outDir, "/meth_snps_fil_matGR_withXY.rds"))
meth_snps_fil_mat <- readRDS("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1/meth_snps_cov_fil_matGR_withXY.rds")

# interact with repeat, CpG island and PMD #####################################
repeat_GR <- rtracklayer::import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/hg38_rmsk.bed", format = "BED") %>%
  plyranges::filter(seqnames %in% c(paste0("chr", c(1:22, "X", "Y"))))

cpgIsland_GR <- rtracklayer::import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/cpgIsland/UCSC_CpGislands_unmask.bed") %>%
  plyranges::filter(seqnames %in%  c(paste0("chr", c(1:22, "X", "Y")))) %>% 
  trim()
plyranges::write_bed(cpgIsland_GR, "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/cpgIsland/UCSC_CpGislands_unmask_fil.bed")

cpgIsland_masked_GR <- rtracklayer::import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/cpgIsland/UCSC_CpGislands.bed") %>%
  plyranges::filter(seqnames %in%  c(paste0("chr", c(1:22, "X", "Y"))))

# for each sample
pmdDNAm_list <- list()
notinPMD_DNAm_list <- list()

for (i in 1:ncol(mcols(meth_snps_fil_mat))) {
  # general info needed
  print(names(mcols(meth_snps_fil_mat))[i])
  ## get the sample DNAm
  cpg_i_GR <- meth_snps_fil_mat[, i]
  
  ## get PMD bed for samples
  pmd_PMMs <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/pmd_allPMM_samples/allPMM_samples_pmd.bed")
  pmd_i_GR <- import(pmd_PMMs) %>%
    plyranges::filter(seqnames %in%  c(paste0("chr", c(1:22, "X", "Y")))) %>%
    plyranges::mutate(pmd_id = paste0("pmd_", 1:length(.)))
  
  # get all PMD DNAm
  pmdDNAm_list[[i]] <- find_overlaps(cpg_i_GR, pmd_i_GR[,0]) %>% 
    plyranges::mutate(beta = mcols(.)[[1]], 
                      samples = rep(names(mcols(meth_snps_fil_mat))[i], 
                                    length(.))) %>% plyranges::select(-1)
  # repeat CpG island not in PMDs
  notinPMD_DNAm_list[[i]] <- filter_by_non_overlaps(cpg_i_GR, pmd_i_GR[,0]) %>% 
    plyranges::mutate(beta = mcols(.)[[1]], 
                      samples = rep(names(mcols(meth_snps_fil_mat))[i], 
                                    length(.))) %>% plyranges::select(-1)
}

# boxplot ######################################################################

# PMD/non-PMD DNAm diff between groups addP ####################################
names(pmdDNAm_list) <- names(mcols(meth_snps_fil_mat))
pmdDNAm_dfl <- lapply(pmdDNAm_list, as.data.frame)
pmdDNAm_df <- do.call(rbind, pmdDNAm_dfl)

pmd_DNAm_plot <- 
  pmdDNAm_df %>%
  # Create the 'group' column by removing digits from 'samples'
  dplyr::mutate(group=str_replace_all(samples, "[:digit:]", "")) %>% 
  na.omit() %>% 
  ggplot(aes(
    x=reorder(group, -beta),
    # x=factor(samples, levels=levels(sample_levels$samples_reordered)), 
    y=beta, fill=group))+
  # geom_density_ridges() +
  stat_boxplot(geom = "errorbar",
               width = 0.15, linetype = 2, # Line type
               lwd = 0.5)+
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE)+
  labs(title = "", x = "", y="DNA methylation")+
  theme_pubr()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(paste0(outDir, "/pmd_DNAm.pdf"),
       plot = pmd_DNAm_plot, 
       device = "pdf", width = 12, height = 10, units = "cm")

# non-PMD DNAm
names(notinPMD_DNAm_list) <- names(mcols(meth_snps_fil_mat))
notinPMD_DNAm_dfl <- lapply(notinPMD_DNAm_list, as.data.frame)
notinPMD_DNAm_df <- do.call(rbind, notinPMD_DNAm_dfl)

notinPMD_DNAm_plot <- 
  notinPMD_DNAm_df %>%
  # Create the 'group' column by removing digits from 'samples'
  dplyr::mutate(group=str_replace_all(samples, "[:digit:]", "")) %>% 
  na.omit() %>% 
  ggplot(aes(
    x=reorder(group, -beta),
    # x=factor(samples, levels=levels(sample_levels$samples_reordered)), 
    y=beta, fill=group))+
  # geom_density_ridges() +
  stat_boxplot(geom = "errorbar",
               width = 0.15, linetype = 2, # Line type
               lwd = 0.5)+
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE)+
  labs(title = "", x = "", y="DNA methylation")+
  theme_pubr()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(paste0(outDir, "/non_pmd_DNAm.pdf"),
       plot = pmd_DNAm_plot, 
       device = "pdf", width = 12, height = 10, units = "cm")



