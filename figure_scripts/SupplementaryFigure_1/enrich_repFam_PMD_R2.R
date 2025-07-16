
options(scipen=999)

# R2: modify enrichment: only use CGI number in PMD for contingency table
# only consider PMD, DNAm BMPC-PMM >0.2


# args    <- commandArgs(trailingOnly = TRUE)

#infile  <- "/data/Project/huyukai/Project_WGS/Job/normal_frags-0329/20220321/sample.info.txt"
# infile  <- args[1]
#outfile <- '/data/Project/huyukai/Project_WGS/Job_new/normal_frags-rerun-0429/chosenmed_ratio'
# outDir <- args[2]
#bampath <- '/data/Project/huyukai/Project_WGS/Job_new/normal_frags-rerun-0429/chosenmed_ratio'
# bampath <- args[3]

# set the library path:
# .libPaths("/home/qwan/miniconda3/envs/coh/lib/R/library")
library(BiocManager)
# BiocManager::install("methrix", "plyranges", "stringr", "rtracklayer", 
#                      "ggplot2", "ggridges", "ggpubr",
#                      "BSgenome.Hsapiens.UCSC.hg38", 
#                      "MafDb.1Kgenomes.phase3.GRCh38")

# in this round3, change to methrix R packge
library(tidyverse)
library(magrittr)
library(methrix)
library(stringr)
library(plyranges)
library(rtracklayer)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(BSgenome.Hsapiens.UCSC.hg38)
# library(MafDb.1Kgenomes.phase3.GRCh38)
library(extrafont)
# Run once to import all system fonts
# font_import(prompt = FALSE)
# Then register the fonts in your session
# loadfonts()
# fonts()[grepl("Arial", fonts(), ignore.case = TRUE)] # no Arial font

outDir <- '/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/pmd_res'

meth_snps_filGR <- readRDS("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1/meth_snps_cov_fil_matGR_withXY.rds")
PMM_metadata <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/PMM_sampleSheet.tsv")

random_regions_GR <- readRDS("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal/res_round4/random_regionsGR_v.rds")

# interact with repeat, CpG island and PMD #####################################
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
repeat_GR <- 
  repeat_bed_df %>%
  dplyr::mutate(seqnames = genoName, 
                start = genoStart + 1, end = genoEnd) %>%
  dplyr::select(-c(genoName:genoEnd)) %>% 
  as_granges()

cpgIsland_GR <- rtracklayer::import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/cpgIsland/UCSC_CpGislands_unmask.bed") %>%
  plyranges::filter(seqnames %in%  c(paste0("chr", c(1:22, "X", "Y")))) %>% 
  plyranges::mutate(cpgIsland_id = paste("cpg", 1:length(.), sep = "_"))

# for each sample
repeatDNAm_list <- list()
cpgIslandDNAm_list <- list()
pmdDNAm_list <- list()

repeatCpGislandDNAm_list <- list()
repeatPMDdnam_list <- list()
cpgIslandPMDdnam_list <- list()
repeatCpGislandPMDdnam_list <- list()
repeat_cpgIsland_notinPMD_DNAm_list <- list()

for (i in 1:ncol(mcols(meth_snps_filGR))) {
  
  # general info needed
  print(names(mcols(meth_snps_filGR))[i])
  ## get the sample DNAm
  cpg_i_GR <- meth_snps_filGR[, i]
  
  # get PMD bed for this matched sample
  # pmd_list <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/*/07.PMD/*_pmd.bed")
  # pattern = paste0(names(mcols(meth_snps_filGR))[i], "_pmd.bed")
  # pmd_i <- pmd_list[str_detect(pmd_list, pattern = pattern)]
  # pmd_i_GR <- import(pmd_i) %>%
  #   plyranges::filter(seqnames %in%  c(paste0("chr", c(1:22, "X", "Y"))))
  # get PMD bed for all PMM sample
  pmd_PMMs <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/pmd_allPMM_samples/allPMM_samples_pmd.bed")
  pmd_i_GR <- import(pmd_PMMs) %>%
    plyranges::filter(seqnames %in%  c(paste0("chr", c(1:22, "X", "Y")))) %>%
    plyranges::mutate(pmd_id = paste0("pmd_", 1:length(.)))
  
  # get all repeat DNAm
  repeatDNAm_list[[i]] <- find_overlaps(cpg_i_GR, repeat_GR) %>% 
    plyranges::mutate(beta = mcols(.)[[1]], 
                      samples = rep(names(mcols(meth_snps_filGR))[i], 
                                    length(.))) %>% plyranges::select(-1)

  # get all CpGisland DNAm
  cpgIslandDNAm_list[[i]] <- find_overlaps(cpg_i_GR, cpgIsland_GR) %>% 
    plyranges::mutate(beta = mcols(.)[[1]], 
                      samples = rep(names(mcols(meth_snps_filGR))[i], 
                                    length(.))) %>% 
    plyranges::select(-1)
  
  # get all PMD DNAm
  pmdDNAm_list[[i]] <- find_overlaps(cpg_i_GR, pmd_i_GR) %>% 
    plyranges::mutate(beta = mcols(.)[[1]], 
                      samples = rep(names(mcols(meth_snps_filGR))[i], 
                                    length(.))) %>% plyranges::select(-1)
  
  # get repeat CpGisland DNAm
  repeatCpGislandDNAm_list[[i]] <- find_overlaps(repeatDNAm_list[[i]], 
                                                 cpgIslandDNAm_list[[i]][, c(-2,-3,-4)])
  
  # get all repeat DNAm in PMDs
  repeatPMDdnam_list[[i]] <- find_overlaps(repeatDNAm_list[[i]], pmd_i_GR)
  
  # get all CpGisland DNAm in PMDs
  cpgIslandPMDdnam_list[[i]] <- find_overlaps(cpgIslandDNAm_list[[i]], pmd_i_GR)
  
  # get all CpG island DNAm in repeat in PMDs
  repeatCpGislandPMDdnam_list[[i]] <- find_overlaps(repeatPMDdnam_list[[i]],
                                                    cpgIslandPMDdnam_list[[i]][, c(-3,-4,-5,-6)])
  # repeat CpG island not in PMDs
  repeat_cpgIsland_notinPMD_DNAm_list[[i]] <- filter_by_non_overlaps(repeatCpGislandDNAm_list[[i]], 
                                                                     pmd_i_GR[,0])
}


## all repeat CGI ##############################################################
### get the number of all repeat (i.e. repeat with CpGs)
names(repeatCpGislandDNAm_list) <- names(mcols(meth_snps_filGR))
repeatFam_num_list <- 
  lapply(repeatCpGislandDNAm_list,
         function(x){
           repFamG <- x %>% mcols() %>% 
             as.data.frame() %>% 
             dplyr::filter(!is.na(beta)) %>% 
             dplyr::select(te_id, repFamily) %>% 
             unique()
           repFamG_num <- repFamG$repFamily %>% table() %>% 
             sort(decreasing = T) %>% as.data.frame()
           colnames(repFamG_num) <-  c("repFamily", "numTE_repFamWithCGI")
           repFamG_num <- repFamG_num %>% 
             dplyr::mutate(total_numTEwithCGI = sum(numTE_repFamWithCGI)) 
           return(repFamG_num)
           })
### Add list names as a new column using lapply
repeatFam_num_dfl <- lapply(names(repeatFam_num_list), function(name) {
  df <- repeatFam_num_list[[name]]          # Extract the data frame
  df$sample <- name                  # Add a new column with the list name
  return(df)
})
repeatFam_num_df <- do.call(rbind, repeatFam_num_dfl)
repeatFam_num_df_fil <- repeatFam_num_df %>% 
  dplyr::select(-sample) %>% 
  unique()


## repeat CpG island in PMD  ###################################################
names(repeatCpGislandPMDdnam_list) <- names(mcols(meth_snps_filGR))

repeatCpGislandPMD_Fam_val_list <- 
  lapply(repeatCpGislandPMDdnam_list,
         function(x){
           repFamG_unique <- x %>% mcols() %>% 
             as.data.frame() %>% 
             dplyr::filter(!is.na(beta)) %>% 
             dplyr::select(te_id, repFamily) %>% 
             unique()
           
           rep_meanBeta <- x %>% mcols() %>% 
             as.data.frame() %>% 
             dplyr::filter(!is.na(beta)) %>% 
             dplyr::group_by(samples, te_id) %>% 
             dplyr::summarise(te_beta = mean(beta)) %>% 
             dplyr::left_join(repFamG_unique, by = "te_id") %>% 
             as.data.frame()
           
           return(rep_meanBeta)
         })
repeatCpGislandPMD_Fam_val_df <- do.call(rbind, repeatCpGislandPMD_Fam_val_list)

teid_repFamilyMatch <- mcols(repeat_GR[,c("te_id", "repFamily")]) %>% 
  as.data.frame() %>% unique()

### select repeat with repeatCpGislandPMD DNAm BMPC-PMM > 0.2 ###################
### VMR
mean_PMM <- repeatCpGislandPMD_Fam_val_df %>% 
  dplyr::filter(str_detect(samples, "PMM")) %>%
  dplyr::group_by(te_id) %>% 
  dplyr::summarise(
    mean_beta_PMM = mean(te_beta, na.rm = TRUE),  # Calculate mean
    sd_beta_PMM = sd(te_beta, na.rm = TRUE)       # Calculate sd
  ) %>% 
  dplyr::mutate(sd_mean_PMM = sd_beta_PMM/mean_beta_PMM)

mean_B <- repeatCpGislandPMD_Fam_val_df %>% 
  dplyr::filter(str_detect(samples, "_rest")) %>%
  dplyr::group_by(te_id) %>% 
  dplyr::summarise(
    mean_beta_B = mean(te_beta, na.rm = TRUE),  # Calculate mean
    sd_beta_B = sd(te_beta, na.rm = TRUE)       # Calculate sd
  ) %>% 
  dplyr::mutate(sd_mean_B = sd_beta_B/mean_beta_B)

rep_select <- cbind(mean_PMM, mean_B[-1]) %>% 
  dplyr::mutate(DNAm_diff = mean_beta_B-mean_beta_PMM) %>% 
  dplyr::filter(DNAm_diff > 0.2) %>%
  dplyr::left_join(teid_repFamilyMatch, "te_id") # n=6735

### VMR in PMD
repeatCpGislandPMD_Fam_val_fil_list <- 
  lapply(repeatCpGislandPMD_Fam_val_list,
         function(x){
           pmd_CGI_num <- x$repFamily %>% table() %>% as.data.frame()
           colnames(pmd_CGI_num) <-  c("repFamily", "numTE_repFam_CGI_PMD")
           pmd_CGI_num <- pmd_CGI_num %>% 
             dplyr::mutate(total_numTEwithCGI_PMD = sum(numTE_repFam_CGI_PMD))
             
           repFamG <- x %>% 
             dplyr::filter(te_id %in% rep_select$te_id)
           sample <- repFamG$samples %>% unique()
           
           repFamG_num <- repFamG$repFamily %>% table() %>% 
             sort(decreasing = T) %>% as.data.frame()
           colnames(repFamG_num) <-  c("repFamily", "numTE_repFam_VMR_PMD")
           repFamG_num <- repFamG_num %>% 
             dplyr::mutate(total_numTE_VMR_PMD = sum(numTE_repFam_VMR_PMD)) %>% 
             dplyr::left_join(pmd_CGI_num, by = "repFamily")
           repFamG_num$sample <- rep(sample, nrow(repFamG_num))
           return(repFamG_num)
         })

repeatCpGislandPMD_Fam_val_fil_df <- 
  do.call(rbind, repeatCpGislandPMD_Fam_val_fil_list)

repeatCpGislandPMD_Fam_val_fil_df_fil <- repeatCpGislandPMD_Fam_val_fil_df %>% 
  remove_rownames() %>% 
  dplyr::select(-sample) %>% 
  unique()

# enrichment analysis
## for repeat subgroups with CpG island in PMDs
# function #####################################################################

# function to apply Fisher's exact test for each row
apply_fisher_test <- function(row) {
  # Extract values
  x <- row$numTE_repFam_VMR_PMD
  n <- row$numTE_repFam_CGI_PMD 
  k <- row$total_numTE_VMR_PMD 
  N <- row$total_numTEwithCGI_PMD 
  # Create contingency table
  contingency_table <- matrix(c(x, # repFam + VMR
                                n - x, # repFam + non-VMR
                                k - x, # otherRepFam + VMR
                                N - k - (n - x)), # otherRepFam + non-VMR
                              nrow = 2, byrow = TRUE)
  # Compute proportions
  # prop_repFam: Proportion of repFam PMD
  # prop_otherRepFam: Proportion of otherRepFam PMD
  # prop_repFam <- x/n
  # prop_otherRepFam <- (k - x)/(N - n)
  # Compute enrichment score
  # enrichment_score <- prop_repFam / prop_otherRepFam
  
  # Calculate proportions
  prop_repFam <- x / k
  prop_expected_repFam <- n / N
  enrichment <- prop_repFam / prop_expected_repFam
  
  # Expected count for repFam in VMR
  expected_count <- k * (n / N)
  
  # Perform Fisher's exact test
  fisher_result <- fisher.test(contingency_table, alternative = "greater")
  # Return p-value and odds ratio
  # print(fisher_result)
  # Return results
  return(c(
    p_value = fisher_result$p.value,
    odds_ratio = fisher_result$estimate,
    enrichment_score = enrichment,
    expected_count = expected_count
  ))
}

te_enrichment <- repeatFam_num_df_fil %>% 
  # dplyr::left_join(repeatCGI_Fam_val_fil_df_fil, by = c("repFamily")) %>% 
  dplyr::left_join(repeatCpGislandPMD_Fam_val_fil_df_fil, by = c("repFamily")) %>%
  mutate(across(where(is.factor), as.character)) %>%  
  mutate_all(~replace(., is.na(.), 0)) %>%  
  dplyr::filter(numTE_repFam_VMR_PMD >=1) %>%
  dplyr::filter(numTE_repFamWithCGI - numTE_repFam_VMR_PMD >=1)

# fisher exact test
results <- sapply(1:nrow(te_enrichment), function(i) apply_fisher_test(te_enrichment[i, ]))
## Transpose the results matrix and set column names
results <- t(results) %>% as.data.frame()
colnames(results) <- c("p_value", "odds_ratio", "enrichment_score", "expected_count")
te_enrichment$p_value <- results$p_value
te_enrichment$odds_ratio <- results$odds_ratio
te_enrichment$enrichment_score <- results$enrichment_score
te_enrichment$expected_count <- results$expected_count

# fdr
te_enrichment$adj_p_value <- p.adjust(te_enrichment$p_value, method = "fdr")

# get retrotransposons
retroTE_fam <- repeat_GR[repeat_GR$repClass%in% 
                           c("LINE", "SINE", "LTR", "Retroposon"),]$repFamily %>% 
  table() %>% names() %>% 
  .[!str_detect(., "\\?")]


p <- te_enrichment %>% 
  dplyr::filter(numTE_repFam_VMR_PMD >1) %>%
  # dplyr::filter(adj_p_value<0.05) %>%
  dplyr::filter(repFamily %in% retroTE_fam) %>% 
  # dplyr::filter(str_detect(sample, "PMM")) %>% 
  # dplyr::mutate(group=str_replace_all(sample, "[:digit:]", "")) %>% 
  ggplot(aes(x=reorder(repFamily,-(enrichment_score)), 
             y=enrichment_score, 
             colour = adj_p_value, 
             size = numTE_repFam_VMR_PMD))+
  geom_point()+
  geom_hline(yintercept = 1, linetype = "dashed",linewidth = 0.6, col="red")+
  # facet_wrap(~group, nrow = 3)+
  labs(x="Retrotransposon family", y = "Observed/Expected", 
       colour = "FDR", 
       size = "Number of TE with CGI\nDNAm difference in PMDs")+
  guides(
    colour = guide_colourbar(position = "right"),
    size   = guide_legend(position = "top")
  ) +
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.9))

ggsave(paste0(outDir, "/repeatFam_PMD_enrichment_retroTE_R2_fdr0.05.pdf"),
       plot = p, 
       device = "pdf", width = 20, height = 12.5, units = "cm")

