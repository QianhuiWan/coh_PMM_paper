
outDir <- '/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1'

library(tidyverse)
library(GenomicRanges)
library(plyranges)
library(rtracklayer)

PMM_metadata <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/PMM_sampleSheet.tsv")

lcpm_PMM <- readr::read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/lcpm_TMM_PMM_STARcounts.tsv")

# S1: get DE genes in IFN, NFKB and p53 pathway MM
keyDEgene_id_symbol <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR1/key_pathway_DEgenes_MM.tsv")


# S2: get promoter regions 
# Define promoter regions around the TSS
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
# get all key DEgene GRanges
keyDE_genes <- getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version", 
                                   "external_gene_name", "description",
                                   "chromosome_name",
                                   "start_position", "end_position", "strand"),
                    # filters = "ensembl_gene_id_version",
                    filters = "ensembl_gene_id",
                    # values = unique(keyDEgene_id_symbol$gene_id_version),  
                    values = unique(keyDEgene_id_symbol$gene_id),  
                    mart = ensembl)

# Convert to GRanges
keyDE_ranges <- GRanges(
  seqnames = keyDE_genes$chromosome_name,
  ranges = IRanges(start = keyDE_genes$start_position, end = keyDE_genes$end_position),
  strand = ifelse(keyDE_genes$strand == 1, "+", "-"),
  gene_name = keyDE_genes$external_gene_name,
  gene_id = keyDE_genes$ensembl_gene_id_version,
  description = keyDE_genes$description
)
promoter_ranges <- promoters(keyDE_ranges, upstream = 2000, downstream = 2000)

# S3: get CpG overlapped promoter regions 
promoter_ranges <- promoter_ranges %>% 
  as.data.frame() %>% 
  plyranges::mutate(seqnames = paste0("chr", seqnames)) %>% 
  as_granges() 
## load CpG islands (unmasked)
cpgIsland_GR <- rtracklayer::import("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/cpgIsland/UCSC_CpGislands_unmask.bed") %>%
  plyranges::filter(seqnames %in%  c(paste0("chr", c(1:22, "X", "Y")))) %>% 
  plyranges::mutate(cpgIsland_id = paste("cpg", 1:length(.), sep = "_"))
# get DMR overlapped KZFP promoter overlapped CpG islands
dmrs <- readRDS(file="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1/pmm_bsseq_tstat_R2_dmrs.rds")
dmrs_GR <- dmrs %>% mutate(seqnames = chr) %>% as_granges()

dmr_keyDEpromoter <- find_overlaps(dmrs_GR, promoter_ranges) #1
dmr_keyDEpromoter$DMR_id <- paste0(dmr_keyDEpromoter$gene_name, paste0("_DMR", 1:length(dmr_keyDEpromoter)))

# find_overlaps(dmr_keyDEpromoter, cpgIsland_GR)
CpGisland_keyDEpromoter <- find_overlaps(cpgIsland_GR, promoter_ranges) #1
CpGisland_keyDEpromoter$CGI_id <- paste0(CpGisland_keyDEpromoter$gene_name, paste0("_CGI", 1:length(CpGisland_keyDEpromoter)))

# S4: save to bed file 

# S5: CpG methylation against gene expression
## get gene expression
PMM_metadata <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/PMM_sampleSheet.tsv")
PMM_metadata$group <- str_replace_all(PMM_metadata$group, "BMPC", "control")
PMM_ss_L1PAnum <- PMM_metadata %>% 
  dplyr::select(sampleName, L1PA2, L1PA3) %>% 
  dplyr::filter(str_detect(sampleName, "PMM"))

lcpm_PMM <- readr::read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_fastp_2passStar_sameAsMMRFpipe/gene_counts/lcpm_TMM_PMM_STARcounts.tsv")
keyDE_lcpm_PMM <- lcpm_PMM %>% 
  dplyr::mutate(gene_id = str_replace_all(gene_name, ".+_", ""),
                gene_name = str_replace_all(gene_name, "_.+", "")) %>% 
  dplyr::filter(gene_name %in% keyDE_ranges$gene_name) %>% 
  mutate(lcpm = ifelse(lcpm < 0, 0, lcpm)) %>%
  dplyr::filter(sample %in%  PMM_ss_L1PAnum$sampleName) %>% 
  left_join(PMM_ss_L1PAnum, by = c("sample"="sampleName")) %>% 
  dplyr::mutate(omics = rep("gene expression", nrow(.))) %>% 
  `colnames<-`(c("gene_name", "sample", "value", "ID", 
                 "L1PA2", "L1PA3", "omics")) %>% 
  dplyr::select(ID, gene_name, sample, value, L1PA2:omics)


## get mean DNAm of CEBPA_DMR1 also a CpG island
meth_snps_filGR <- readRDS("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/DNAm/DNAm_hg38_p14_abismal_round3_filter/res_round1/meth_snps_cov_fil_matGR_withXY.rds")
nms <- colnames(values(meth_snps_filGR))
cols <- lapply(nms, function(.) { col <- ensym(.); quo(mean(!!col)) })
# the column names for the summarize op
names(cols) <- paste0("mean_", nms)


CpGisland_keyDEpromoter_df <- 
  # dmr_keyDEpromoter %>% 
  CpGisland_keyDEpromoter %>% 
  find_overlaps(meth_snps_filGR) %>% 
  plyranges::group_by(CGI_id) %>% 
  plyranges::summarise(!!! cols) %>% 
  as.data.frame() %>% 
  `colnames<-`(str_replace_all(colnames(.), "mean_", "")) %>% 
  as.data.frame() %>% 
  dplyr::select(-starts_with("B")) %>% 
  pivot_longer(cols = c(starts_with("PMM")),
               names_to = "sample",
               values_to = "CGI_DNAm") %>%
  dplyr::filter(sample %in%  PMM_ss_L1PAnum$sampleName) %>% 
  left_join(PMM_ss_L1PAnum, by = c("sample"="sampleName")) %>% 
  dplyr::mutate(omics = rep("DNA methylation", nrow(.))) %>% 
  dplyr::mutate(gene_name = str_remove_all(CGI_id, "_..+")) %>% 
  `colnames<-`(c("ID", "sample", "value", 
                 "L1PA2", "L1PA3", "omics","gene_name")) %>% 
  dplyr::select(ID, gene_name, sample, value, L1PA2:omics)

# sig: APP_CGI39, SECTM1_CGI33
p1 <- 
  CpGisland_keyDEpromoter_df %>% 
  # rbind(keyDE_lcpm_PMM) %>% 
  dplyr::filter(ID %in% c("APP_CGI39", "SECTM1_CGI33", "CEBPA_CGI34")) %>% 
  ggplot(aes(x = L1PA2, y = value)) + 
  geom_point() + geom_smooth(method = "lm")+ 
  facet_wrap(~ID, scales = "free_x")+
  ggpubr::stat_cor(method = "pearson", size=3, label.y = 1.2) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3)+
  labs(x = "number of L1PA2", y="DNA methylation")+
  theme_set(ggpubr::theme_pubr(base_size = 12))

 
ggsave(paste0(outDir, "/keyDE_DNAm_L1PA2_cor_v2.pdf"),
       plot = p1, 
       device = "pdf", 
       width = 18, height = 10,
       # width = 22, height = 35, 
       units = "cm")  


# sig: APP_CGI39, SECTM1_CGI33
p2 <- 
  keyDE_lcpm_PMM %>% 
  # rbind(keyDE_lcpm_PMM) %>% 
  dplyr::filter(gene_name %in% c("APP", "SECTM1", "CEBPA")) %>% 
  ggplot(aes(x = L1PA2, y = value)) + 
  geom_point() + geom_smooth(method = "lm")+ 
  facet_wrap(~gene_name, scales = "free_x")+
  ggpubr::stat_cor(method = "pearson", size=3, label.y = 10) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3)+
  labs(x = "number of L1PA2", y="Gene expression")+
  theme_set(ggpubr::theme_pubr(base_size = 12))

ggsave(paste0(outDir, "/keyDE_geneExpre_L1PA2_cor_v2.pdf"),
       plot = p2, 
       device = "pdf", 
       width = 18, height = 10,
       # width = 22, height = 35, 
       units = "cm")  

