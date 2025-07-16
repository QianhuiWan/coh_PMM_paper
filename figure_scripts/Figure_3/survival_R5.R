
# only use BM primary samples
.libPaths()
library(tidyverse)
library(TCGAbiolinks)
library(edgeR)
library(survival)
library(survminer)
options(scipen=999)

out_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess"

all_data <- readRDS(paste0(out_dir, "/all_data_geneCounts_MMRF.rds"))
clinical_data <- all_data@colData %>% as.data.frame()

clinical_data_BM <- clinical_data %>% as.data.frame() %>% 
  dplyr::filter(site_of_resection_or_biopsy == "Bone marrow") %>% 
  dplyr::filter(specimen_type == "Bone Marrow NOS")

clinical_data_BMprimary <- clinical_data %>% as.data.frame() %>% 
  dplyr::filter(site_of_resection_or_biopsy == "Bone marrow") %>% 
  dplyr::filter(specimen_type == "Bone Marrow NOS") %>% 
  dplyr::filter(tumor_descriptor =="Primary") #764

# survival data
# OS.time: OS.time: Duration from the starting point to the event or censoring.
# OS: Binary indicator of whether the event occurred (1) or was censored (0).
# get ZNF141 group

ZNF141_info <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR1/ZNF141_info.tsv")

patientsIn4Groups <- 
  read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_redo_res/all_activeL1PA2_patientsIn4Groups_764.tsv")


# get MMRF patient ID
MMRF_ids_raw <- read_csv("/net/nfs-irwrsrchnas01/labs/dschones/Seq/MMRF-dbgap/SraRunTableRNABoneMarrow.csv", col_names = FALSE)
MMRF_ids <- MMRF_ids_raw[, c("X1", "X10", "X39")]
colnames(MMRF_ids) <- c("SRR_id", "patient_barCode", "patient_id")
MMRF_ids$barcode <- str_remove_all(MMRF_ids$patient_barCode, "_T..+")


clinical_data_addIDs <- clinical_data_BMprimary %>% as.data.frame() %>% 
  # dplyr::filter(specimen_type == "Bone Marrow NOS") %>% 
  left_join(MMRF_ids, by = "barcode") %>% 
  left_join(ZNF141_info, by = c("SRR_id"="sampleID")) %>% 
  left_join(patientsIn4Groups, by = c("SRR_id"="sampleID")) %>% 
  dplyr::mutate(OS = ifelse(vital_status=="Alive", 0, 1)) %>% 
  arrange(L1PA2) %>% 
  # arrange(TEtranscript_num) %>% 
  dplyr::mutate(
    sampleID = factor(SRR_id, levels = SRR_id),
    sample_index = as.numeric(factor(sampleID)),  # Create a numeric index
    quantile_group = cut(
      as.numeric(factor(sampleID)), 
      breaks = quantile(as.numeric(factor(sampleID)), 
                        probs = seq(0, 1, 0.25), na.rm = TRUE), 
      include.lowest = TRUE, 
      labels = c("Q1", "Q2", "Q3", "Q4")
    )) %>% 
  mutate(L1_quantile_group = cut(
    L1PA2,
    breaks = quantile(L1PA2, probs = seq(0, 1, 0.25), na.rm = TRUE),
    include.lowest = TRUE,
    labels = c("Q1", "Q2", "Q3", "Q4")
  ))

# fit surv
fit <- survfit(Surv(days_to_last_known_disease_status, OS) ~ L1_quantile_group, 
               data = clinical_data_addIDs)
# Plot Kaplan-Meier curve
p1 <- ggsurvplot(fit, data = clinical_data_addIDs, conf.int = TRUE, pval = TRUE,
                 # fun = "pct", 
                 # test.for.trend=TRUE,
                 legend.labs = c("Q1", "Q2", "Q3", "Q4"),
                 xlab = "Time (days)",  
                 surv.median.line = "hv")
ggsave(paste0(out_dir, "/L1PA2_quantile_survivalRate_R5.pdf"),
       plot = p1$plot,
       device = "pdf", width = 10, height = 8, units = "cm")


fit2 <- survfit(Surv(days_to_last_known_disease_status, OS) ~ ZNF141group, 
               data = clinical_data_addIDs)
p2 <- ggsurvplot(fit2, data = clinical_data_addIDs, 
                 legend.labs = c("High_ZNF141", "Medium_ZNF141", "Zero_ZNF141"),
                 conf.int = TRUE, pval = TRUE)
# ggsave(paste0(out_dir, "/ZNF141_3groups_survivalRate.pdf"), 
#        plot = p2$plot, 
#        device = "pdf", width = 15, height = 10, units = "cm")

# sink()  # Turn off logging


# get TCGA count mat, match barcode to clinial data
raw_counts <- read_tsv(paste0(out_dir, "/raw_counts_geneCounts_MMRF_BMprimary.tsv")) %>% 
  tibble::column_to_rownames(var = "gene_id")

raw_counts_df_BMprimary <- raw_counts
## create DGElist
### get group info.
if (unique(clinical_data_addIDs$barcode == colnames(raw_counts_df_BMprimary))) {
  groups <- clinical_data_addIDs$L1_quantile_group
  groups[is.na(groups)] <- "nogroup"
}

dge_rawCounts <- DGEList(counts=raw_counts_df_BMprimary, 
                         # genes=raw_counts_df_BMprimary[,1:3], 
                         #remove.zeros = TRUE,
                         group = groups)


## normalize counts with TMM method
dge_rawCounts <- calcNormFactors(dge_rawCounts, method = "TMM")
normalized_lcpm <- cpm(dge_rawCounts, normalized.lib.sizes = TRUE, log = TRUE)

boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF141"),] ~groups)
boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF382"),] ~groups)
boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF649_"),] ~groups)
boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF84_"),] ~groups)
boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF93"),] ~groups)

boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF141"),] ~clinical_data_addIDs$quantile_group)
boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF382"),] ~clinical_data_addIDs$quantile_group)
boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF649_"),] ~clinical_data_addIDs$quantile_group)
boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF84_"),] ~clinical_data_addIDs$quantile_group)
boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF93"),] ~clinical_data_addIDs$quantile_group)

boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF765_"),] ~groups)
boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF17_"),] ~groups)
boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF433_"),] ~groups)

boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF433_"),] ~clinical_data_addIDs$iss_stage)
boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF382_"),] ~clinical_data_addIDs$iss_stage)
boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF141_"),] ~clinical_data_addIDs$iss_stage)
boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "MPHOSPH8_"),] ~clinical_data_addIDs$iss_stage)

boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF84_"),] ~clinical_data_addIDs$iss_stage)
boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF17_"),] ~clinical_data_addIDs$iss_stage)


boxplot(clinical_data_addIDs$L1PA2 ~clinical_data_addIDs$iss_stage)

clinical_data_addIDs %>% dplyr::select(sampleID, iss_stage, L1_quantile_group, L1PA2) %>% 
  dplyr::mutate(group = paste0(iss_stage ,"_", L1_quantile_group)) %>% 
  ggplot(aes(x = group, fill = L1_quantile_group)) + geom_bar() +
  ggpubr::theme_pubr()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


plot(as.vector(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF649_"),]) ~clinical_data_addIDs$ZNF141)
plot(as.vector(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF84_"),]) ~clinical_data_addIDs$ZNF141)
plot(as.vector(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF382_"),]) ~clinical_data_addIDs$ZNF141)
plot(as.vector(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF93_"),]) ~clinical_data_addIDs$ZNF141)
plot(as.vector(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF17_"),]) ~clinical_data_addIDs$ZNF141)


## get pub plot for ZNF17
# better to add the 10 KZFPs without gene symbol
kzfps <- read_tsv(file = paste0("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/kzfps_hg38_allHaveGeneSymbol.tsv")) %>% 
  dplyr::filter(!(ensembl_gene_id %in% c("ENSG00000285124", 
                                         "ENSG00000292386", "ENSG00000274195",
                                         "ENSG00000285286"))) 

for (ZNF in unique(kzfps$hgnc_symbol)[370]) {
  ZNF = unique(kzfps$hgnc_symbol)[370]
  gene_symbol <- str_replace_all(rownames(normalized_lcpm), "_..+", "")
  if (any(gene_symbol %in% ZNF) == FALSE) {
    print(paste0(ZNF, " not detected"))
  }else{
    ZNF_info <- 
      # normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF382_"),] %>% 
      # normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF93_"),] %>%
      # normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF17_"),] %>%
      # normalized_lcpm[str_detect(rownames(normalized_lcpm), paste0(ZNF, "_")),] %>%
      normalized_lcpm[gene_symbol %in% ZNF,] %>%
      as.data.frame() %>% 
      `colnames<-`("ZNF_lcpm") %>% 
      rownames_to_column(var = "barcode") %>% 
      # dplyr::mutate(ZNF382group = ifelse(lcpm>0, "haveZNF382", "NoZNF382"))
      arrange(ZNF_lcpm) 
  
    breaks <- unique(quantile(ZNF_info$ZNF_lcpm, probs = seq(0, 1, length.out = 5), na.rm = TRUE))
    # breaks <- unique(quantile(ZNF_info$ZNF_lcpm, probs = seq(0, 1, 0.25), na.rm = TRUE))
    ZNF_info <- ZNF_info%>% 
      mutate(ZNF_quantile_group = cut(
        ZNF_lcpm,
        # breaks = quantile(ZNF_lcpm, probs = seq(0, 1, 0.25), na.rm = TRUE),
        breaks = breaks,
        include.lowest = TRUE, 
        labels = c("Q1", "Q2", "Q3", "Q4")
      ))
    
  clinical_data_addIDs_addZNF <- 
    clinical_data_addIDs %>% 
    dplyr::left_join(ZNF_info, by = "barcode")
    
  fit3 <- survfit(Surv(days_to_last_known_disease_status, OS) ~ ZNF_quantile_group, 
                  data = clinical_data_addIDs_addZNF)
  p3 <- ggsurvplot(fit3, data = clinical_data_addIDs_addZNF, 
                   legend.labs = c("ZNF17_Q1", "ZNF17_Q2", "ZNF17_Q3", "ZNF17_Q4"),
                   # palette = c("purple", "lightcyan", "darkgreen", "red"),
                   palette = c("#957DAD", "#56B4E9", "#009E73", "#E57373"),
                   xlab = "Time (days)",  
                   surv.median.line = "hv",
                   conf.int = TRUE, pval = TRUE)
  ggsave(paste0(out_dir, "/", ZNF, "_4groups_survivalRate.pdf"),
         plot = p3$plot,
         device = "pdf", width = 15, height = 10, units = "cm")
  }
}


# check E2F high samples
e2fHigh_samples_MM <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR1/e2fHigh_samples_MM.tsv")

table(colnames(normalized_lcpm)==clinical_data_addIDs$barcode)
colnames(normalized_lcpm) <- as.character(clinical_data_addIDs$sampleID)

normalized_lcpm_e2fHighMM <- normalized_lcpm[,colnames(normalized_lcpm)%in%colnames(e2fHigh_samples_MM)]
clinical_data_addIDs_e2fHighMM <- clinical_data_addIDs[clinical_data_addIDs$sampleID %in% colnames(e2fHigh_samples_MM),]
table(colnames(normalized_lcpm_e2fHighMM)==clinical_data_addIDs_e2fHighMM$sampleID)



boxplot(normalized_lcpm_e2fHighMM[str_detect(rownames(normalized_lcpm_e2fHighMM), "ZNF433_"),] ~clinical_data_addIDs_e2fHighMM$iss_stage)
boxplot(normalized_lcpm_e2fHighMM[str_detect(rownames(normalized_lcpm_e2fHighMM), "ZNF382_"),] ~clinical_data_addIDs_e2fHighMM$iss_stage)
boxplot(normalized_lcpm_e2fHighMM[str_detect(rownames(normalized_lcpm_e2fHighMM), "ZNF141_"),] ~clinical_data_addIDs_e2fHighMM$iss_stage)
boxplot(normalized_lcpm_e2fHighMM[str_detect(rownames(normalized_lcpm_e2fHighMM), "MPHOSPH8_"),] ~clinical_data_addIDs_e2fHighMM$iss_stage)
boxplot(normalized_lcpm_e2fHighMM[str_detect(rownames(normalized_lcpm_e2fHighMM), "ZNF17_"),] ~clinical_data_addIDs_e2fHighMM$iss_stage)

boxplot(normalized_lcpm_e2fHighMM[str_detect(rownames(normalized_lcpm_e2fHighMM), "ZNF84_"),] ~clinical_data_addIDs_e2fHighMM$iss_stage)



boxplot(clinical_data_addIDs_e2fHighMM$L1PA2 ~clinical_data_addIDs_e2fHighMM$iss_stage)

plot(normalized_lcpm_e2fHighMM[str_detect(rownames(normalized_lcpm_e2fHighMM), "ZNF141_"),] ~clinical_data_addIDs_e2fHighMM$L1PA2)
# Fit a linear regression model
lm_fit <- lm(normalized_lcpm_e2fHighMM[str_detect(rownames(normalized_lcpm_e2fHighMM), "ZNF141_"),] ~clinical_data_addIDs_e2fHighMM$L1PA2)
summary(lm_fit)
# Add the regression line
abline(lm_fit, col = "red", lwd = 2)

plot(normalized_lcpm_e2fHighMM[str_detect(rownames(normalized_lcpm_e2fHighMM), "ZNF382_"),] ~clinical_data_addIDs_e2fHighMM$L1PA2)
lm_fit <- lm(normalized_lcpm_e2fHighMM[str_detect(rownames(normalized_lcpm_e2fHighMM), "ZNF382_"),] ~clinical_data_addIDs_e2fHighMM$L1PA2)
summary(lm_fit)
# Add the regression line
abline(lm_fit, col = "red", lwd = 2)

plot(normalized_lcpm_e2fHighMM[str_detect(rownames(normalized_lcpm_e2fHighMM), "ZNF17_"),] ~clinical_data_addIDs_e2fHighMM$L1PA2)
# Fit a linear regression model
lm_fit <- lm(normalized_lcpm_e2fHighMM[str_detect(rownames(normalized_lcpm_e2fHighMM), "ZNF17_"),] ~clinical_data_addIDs_e2fHighMM$L1PA2)
summary(lm_fit)
# Add the regression line
abline(lm_fit, col = "red", lwd = 2)



