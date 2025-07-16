
# args    <- commandArgs(trailingOnly = TRUE)

#infile  <- "/data/Project/huyukai/Project_WGS/Job/normal_frags-0329/20220321/sample.info.txt"
# infile  <- args[1]
#outfile <- '/data/Project/huyukai/Project_WGS/Job_new/normal_frags-rerun-0429/chosenmed_ratio'
# outDir <- args[2]
# bampath <- '/data/Project/huyukai/Project_WGS/Job_new/normal_frags-rerun-0429/chosenmed_ratio'
# bampath <- args[3]

# set the library path:
# .libPaths("/home/qwan/miniconda3/envs/coh/lib/R/library")
# library(BiocManager)
# BiocManager::install("methrix", "plyranges", "stringr", "rtracklayer", 
#                      "ggplot2", "ggridges", "ggpubr",
#                      "BSgenome.Hsapiens.UCSC.hg38", 
#                      "MafDb.1Kgenomes.phase3.GRCh38")

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
  plyranges::filter(seqnames %in%  c(paste0("chr", c(1:22, "X", "Y"))))

cpgIsland_masked_GR <- rtracklayer::import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/cpgIsland/UCSC_CpGislands.bed") %>%
  plyranges::filter(seqnames %in%  c(paste0("chr", c(1:22, "X", "Y"))))

# for each sample
repeatDNAm_list <- list()
cpgIslandDNAm_list <- list()
cpgIsland_maskedDNAm_list <- list()
pmdDNAm_list <- list()

repeatCpGislandDNAm_list <- list()
repeatPMDdnam_list <- list()
cpgIslandPMDdnam_list <- list()
repeatCpGislandPMDdnam_list <- list()
repeat_cpgIsland_notinPMD_DNAm_list <- list()

for (i in 1:ncol(mcols(meth_snps_fil_mat))) {
  # general info needed
  print(names(mcols(meth_snps_fil_mat))[i])
  ## get the sample DNAm
  cpg_i_GR <- meth_snps_fil_mat[, i]
  
  ## get PMD bed for this matched sample
  # pmd_list <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/*/07.PMD/*_pmd.bed")
  # pattern = paste0(names(mcols(meth_snps_fil_mat))[i], "_pmd.bed")
  # pmd_i <- pmd_list[str_detect(pmd_list, pattern = pattern)]
  # pmd_i_GR <- import(pmd_i) %>% 
  #   plyranges::filter(seqnames %in%  c(paste0("chr", c(1:22, "X", "Y"))))
  
  pmd_PMMs <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/pmd_allPMM_samples/allPMM_samples_pmd.bed")
  pmd_i_GR <- import(pmd_PMMs) %>%
    plyranges::filter(seqnames %in%  c(paste0("chr", c(1:22, "X", "Y")))) %>%
    plyranges::mutate(pmd_id = paste0("pmd_", 1:length(.)))
  
  # get all repeat DNAm
  repeatDNAm_list[[i]] <- find_overlaps(cpg_i_GR, repeat_GR) %>% 
    plyranges::mutate(beta = mcols(.)[[1]], 
                      samples = rep(names(mcols(meth_snps_fil_mat))[i], 
                                    length(.))) %>% plyranges::select(-1)

  # get all CpGisland DNAm
  cpgIslandDNAm_list[[i]] <- find_overlaps(cpg_i_GR, cpgIsland_GR) %>% 
    plyranges::mutate(beta = mcols(.)[[1]], 
                      samples = rep(names(mcols(meth_snps_fil_mat))[i], 
                                    length(.))) %>% plyranges::select(-1)
  
  cpgIsland_maskedDNAm_list[[i]] <- find_overlaps(cpg_i_GR, cpgIsland_masked_GR) %>% 
    plyranges::mutate(beta = mcols(.)[[1]], 
                      samples = rep(names(mcols(meth_snps_fil_mat))[i], 
                                    length(.))) %>% plyranges::select(-1)
  
  # get all PMD DNAm
  pmdDNAm_list[[i]] <- find_overlaps(cpg_i_GR, pmd_i_GR) %>% 
    plyranges::mutate(beta = mcols(.)[[1]], 
                      samples = rep(names(mcols(meth_snps_fil_mat))[i], 
                                    length(.))) %>% plyranges::select(-1)
  
  # get repeat CpGisland DNAm
  repeatCpGislandDNAm_list[[i]] <- find_overlaps(repeatDNAm_list[[i]], 
                                                 cpgIslandDNAm_list[[i]][, c(-2,-3)])
  
  # get all repeat DNAm in PMDs
  repeatPMDdnam_list[[i]] <- find_overlaps(repeatDNAm_list[[i]], pmd_i_GR)
  
  # get all CpGisland DNAm in PMDs
  cpgIslandPMDdnam_list[[i]] <- find_overlaps(cpgIslandDNAm_list[[i]], pmd_i_GR)
  
  # get all CpG island DNAm in repeat in PMDs
  repeatCpGislandPMDdnam_list[[i]] <- find_overlaps(repeatPMDdnam_list[[i]],
                                                    cpgIslandPMDdnam_list[[i]][, c(-2,-3,-4,-5)])
  # repeat CpG island not in PMDs
  repeat_cpgIsland_notinPMD_DNAm_list[[i]] <- filter_by_non_overlaps(repeatCpGislandDNAm_list[[i]], 
                                                                     pmd_i_GR[,0])
}


# boxplot
## masked CpG island
names(cpgIsland_maskedDNAm_list) <- names(mcols(meth_snps_fil_mat))
cpgIslandDNAm_dfl <- lapply(cpgIsland_maskedDNAm_list, as.data.frame)
cpgIslandDNAm_df <- do.call(rbind, cpgIslandDNAm_dfl)
# cpgIslandDNAm_df <- bind_rows(cpgIslandDNAm_dfl, .id = "basename")
# nonrepeat_cpgIslandDNAm_list <- lapply(cpgIslandDNAm_list,
#                           function(x){filter_by_non_overlaps(x, repeat_GR)})
# nonrepeat_cpgIslandDNAm_dfl <- lapply(nonrepeat_cpgIslandDNAm_list, as.data.frame)
# nonrepeat_cpgIslandDNAm_df <- do.call(rbind, nonrepeat_cpgIslandDNAm_dfl)
library(tidytext)

sample_levels <- cpgIslandDNAm_df %>% 
  dplyr::mutate(group=str_replace_all(samples, "[:digit:]", "")) %>%
  dplyr::filter(group != "B_EBV") %>% 
  # na.omit() %>% 
  group_by(group) %>%
  arrange(beta, .by_group = TRUE) %>%
  mutate(samples_reordered = factor(samples, levels = unique(samples)))
  
masked_cpgIslandDNAm_plot <- 
  cpgIslandDNAm_df %>% 
  # Create the 'group' column by removing digits from 'samples'
  dplyr::mutate(group=str_replace_all(samples, "[:digit:]", "")) %>% 
  dplyr::filter(group != "B_EBV") %>% 
  # na.omit() %>% 
  ggplot(aes(x=factor(samples, levels=levels(sample_levels$samples_reordered)), 
             y=beta, fill=group))+
  # geom_density_ridges() +
  stat_boxplot(geom = "errorbar",
               width = 0.15, linetype = 2, # Line type
               lwd = 0.5)+
  geom_boxplot(outlier.shape = NA)+
  labs(title = "CGI", x = "", y="DNA methylation")+
  theme_pubr()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(paste0(outDir, "/masked_cpgIslandDNAm.pdf"),
       plot = masked_cpgIslandDNAm_plot, 
       device = "pdf", width = 12, height = 10, units = "cm")

## repeat CpG island in PMD
names(repeatCpGislandPMDdnam_list) <- names(mcols(meth_snps_fil_mat))
repeatCpGislandPMDdnam_dfl <- lapply(repeatCpGislandPMDdnam_list, as.data.frame)
repeatCpGislandPMDdnam_df <- do.call(rbind, repeatCpGislandPMDdnam_dfl)

repeatCpGisland_PMD_DNAm_plot <- 
  repeatCpGislandPMDdnam_df %>%
  # Create the 'group' column by removing digits from 'samples'
  dplyr::mutate(group=str_replace_all(samples, "[:digit:]", "")) %>% 
  dplyr::filter(group != "B_EBV") %>% 
  na.omit() %>% 
  ggplot(aes(
    x=reorder(samples, -beta),
    # x=factor(samples, levels=levels(sample_levels$samples_reordered)), 
    y=beta, fill=group))+
  # geom_density_ridges() +
  stat_boxplot(geom = "errorbar",
               width = 0.15, linetype = 2, # Line type
               lwd = 0.5)+
  geom_boxplot(outlier.shape = NA)+
  labs(title = "Repeat CGI in PMDs", x = "", y="DNA methylation")+
  theme_pubr()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(paste0(outDir, "/repeatCpGislandPMDdnam.pdf"),
       plot = repeatCpGisland_PMD_DNAm_plot, 
       device = "pdf", width = 12, height = 10, units = "cm")

## repeat CpG island not in PMD
names(repeat_cpgIsland_notinPMD_DNAm_list) <- names(mcols(meth_snps_fil_mat))
repeat_cpgIsland_notinPMD_dfl <- lapply(repeat_cpgIsland_notinPMD_DNAm_list, as.data.frame)
repeat_cpgIsland_notinPMD_df <- do.call(rbind, repeat_cpgIsland_notinPMD_dfl)

repeatCpGisland_notInPMD_DNAm_plot <- 
  repeat_cpgIsland_notinPMD_df %>%
  # Create the 'group' column by removing digits from 'samples'
  dplyr::mutate(group=str_replace_all(samples, "[:digit:]", "")) %>% 
  dplyr::filter(group != "B_EBV") %>% 
  na.omit() %>% 
  ggplot(aes(
    x=reorder(samples, -beta),
    # x=factor(samples, levels=levels(sample_levels$samples_reordered)), 
    y=beta, fill=group))+
  # geom_density_ridges() +
  stat_boxplot(geom = "errorbar",
               width = 0.15, linetype = 2, # Line type
               lwd = 0.5)+
  geom_boxplot(outlier.shape = NA)+
  labs(title = "Repeat CGI not in PMDs", x = "", y="DNA methylation")+
  theme_pubr()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(paste0(outDir, "/repeatCpGislandNotInPMDdnam.pdf"),
       plot = repeatCpGisland_notInPMD_DNAm_plot, 
       device = "pdf", width = 12, height = 10, units = "cm")




# non-repeat CpG-island in PMDs
names(cpgIslandPMDdnam_list) <- names(mcols(meth_snps_fil_mat))

## unmasked CGI
nonrepeat_cpgIsland_inPMD_DNAm_list <- 
  lapply(cpgIslandPMDdnam_list, 
         function(x){filter_by_non_overlaps(x, repeat_GR)})
nonrepeat_cpgIslandinPMD_DNAm_dfl <- lapply(nonrepeat_cpgIsland_inPMD_DNAm_list, as.data.frame)
nonrepeat_cpgIslandinPMD_DNAm_df <- do.call(rbind, nonrepeat_cpgIslandinPMD_DNAm_dfl)
## masked CGI
nonrepeat_cpgIsland_inPMD_DNAm_list <- lapply(
  cpgIsland_maskedDNAm_list,
  function(x) {
    x_filtered <- filter_by_overlaps(x, pmd_i_GR)
    x_nonrepeat <- filter_by_non_overlaps(x_filtered, repeat_GR)
    return(x_nonrepeat)
  }
)

nonrepeat_cpgIslandinPMD_DNAm_dfl <- lapply(nonrepeat_cpgIsland_inPMD_DNAm_list, as.data.frame)
nonrepeat_cpgIslandinPMD_DNAm_df <- do.call(rbind, nonrepeat_cpgIslandinPMD_DNAm_dfl)

sample_levels <- c("B3_rest", "B2_rest", "B1_rest", 
                   "PMM11", "PMM14", "PMM1", "PMM13", "PMM3", 
                   "PMM9", "PMM15", "PMM4", "PMM12", "PMM2", "PMM16", "PMM18",
                   "PMM17", "PMM6", "PMM7")

nonrepeat_cpgIslandinPMD_DNAm_plot <- 
  nonrepeat_cpgIslandinPMD_DNAm_df %>%
  # Create the 'group' column by removing digits from 'samples'
  dplyr::mutate(group=str_replace_all(samples, "[:digit:]", "")) %>% 
  dplyr::filter(group != "B_EBV") %>% 
  na.omit() %>% 
  dplyr::mutate(samples = factor(samples, levels = unique(samples))) %>% 
  ggplot(aes(
    # x=reorder(samples, -beta),
    x=factor(samples, levels=sample_levels),
    y=beta, fill=group))+
  # geom_density_ridges() +
  stat_boxplot(geom = "errorbar",
               width = 0.15, linetype = 2, # Line type
               lwd = 0.5)+
  geom_boxplot(outlier.shape = NA)+
  labs(title = "Non-repeat masked CGI in PMDs", x = "", y="DNA methylation")+
  theme_pubr()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(paste0(outDir, "/nonrepeat_masked_cpgIslandinPMD_DNAm.pdf"),
       plot = nonrepeat_cpgIslandinPMD_DNAm_plot, 
       device = "pdf", width = 12, height = 10, units = "cm")



# non-repeat CpG-island not in PMDs
## unmasked CGI
nonrepeat_cpgIsland_notinPMD_DNAm_list <- lapply(
  cpgIslandDNAm_list,
  function(x) {
    x_nonPMD <- filter_by_non_overlaps(x, pmd_i_GR)
    x_nonrepeat <- filter_by_non_overlaps(x_nonPMD, repeat_GR)
    return(x_nonrepeat)
  }
)
nonrepeat_cpgIsland_notinPMD_DNAm_dfl <- lapply(nonrepeat_cpgIsland_notinPMD_DNAm_list, as.data.frame)
nonrepeat_cpgIsland_notinPMD_DNAm_df<- do.call(rbind, nonrepeat_cpgIsland_notinPMD_DNAm_dfl)

sample_levels <- c("B3_rest", "B1_rest", "B2_rest",
                   "PMM11", "PMM1", "PMM15", "PMM14", "PMM4", "PMM13", 
                   "PMM9", "PMM3",  "PMM2", "PMM16", "PMM18",
                   "PMM12", "PMM17", "PMM6", "PMM7")

## masked CGI
nonrepeat_cpgIsland_notinPMD_DNAm_list <- lapply(
  cpgIsland_maskedDNAm_list,
  function(x) {
    x_nonPMD <- filter_by_non_overlaps(x, pmd_i_GR)
    x_nonrepeat <- filter_by_non_overlaps(x_nonPMD, repeat_GR)
    return(x_nonrepeat)
  }
)
nonrepeat_cpgIsland_notinPMD_DNAm_dfl <- lapply(nonrepeat_cpgIsland_notinPMD_DNAm_list, as.data.frame)
nonrepeat_cpgIsland_notinPMD_DNAm_df <- do.call(rbind, nonrepeat_cpgIsland_notinPMD_DNAm_dfl)

sample_levels <- c("B3_rest", "B1_rest", "B2_rest", 
                   "PMM11", "PMM12",  "PMM1", "PMM15", "PMM14", "PMM4", "PMM13", 
                   "PMM9", "PMM2", "PMM3", "PMM16", "PMM6", "PMM18",
                    "PMM17", "PMM7")

nonrepeat_cpgIsland_notinPMD_DNAm_plot <- 
  nonrepeat_cpgIsland_notinPMD_DNAm_df %>%
  # Create the 'group' column by removing digits from 'samples'
  dplyr::mutate(group=str_replace_all(samples, "[:digit:]", "")) %>% 
  dplyr::filter(group != "B_EBV") %>% 
  na.omit() %>% 
  ggplot(aes(
    # x=reorder(samples, -beta),
    x=factor(samples, levels=sample_levels),
    y=beta, fill=group))+
  # geom_density_ridges() +
  stat_boxplot(geom = "errorbar",
               width = 0.15, linetype = 2, # Line type
               lwd = 0.5)+
  geom_boxplot(outlier.shape = NA)+
  labs(title = "Non-repeat masked CGI not in PMDs", x = "", y="DNA methylation")+
  theme_pubr()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(paste0(outDir, "/nonRepeatCpGisland_masked_NotInPMDdnam.pdf"),
       plot = nonrepeat_cpgIsland_notinPMD_DNAm_plot, 
       device = "pdf", width = 12, height = 10, units = "cm")



