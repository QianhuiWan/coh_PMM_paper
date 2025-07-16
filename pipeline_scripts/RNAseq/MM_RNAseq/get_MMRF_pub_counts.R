
.libPaths()
library(tidyverse)
library(TCGAbiolinks)
library(survival)
library(survminer)
options(scipen=999)

out_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess"

# Query for all MMRF projects with RNA-Seq data
projects <- TCGAbiolinks::getGDCprojects()
mmrf_projects <- projects[grepl("MMRF", projects$project_id), "project_id"] #n=1

# sink("TCGA_countMat_download.log")
tumor_normal_projects <- list()
tumor_normal_project_expQuery <- list()
# Loop through each project to check and download both tumor and normal samples
for (project in mmrf_projects) {
  # counts query
  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    # data.type = c("Aligned Reads", "Unaligned Reads"),
    data.type = "Gene Expression Quantification"
    # workflow.type = "STAR - Counts",
    # sample.type = c("Primary Tumor", "Solid Tissue Normal")
    # sample.type = c("Solid Tissue Normal")
  )
  
  # Check for available results before downloading
  results <- getResults(query)
  if (!is.null(results) && nrow(results) > 0) {
    # Download the data to the specific project directory
    GDCdownload(query = query, directory = out_dir, files.per.chunk=4)
    # Add project to the list of completed downloads
    tumor_normal_projects <- c(tumor_normal_projects, project)
    tumor_normal_project_expQuery[[project]] <- query
  } else {
    message("No results found for project: ", project)
  }
}


# Display the projects with both tumor and normal samples
tumor_normal_projects
tumor_normal_project_expQuery


# clinical_data <- getResults(tumor_normal_project_expQuery$`MMRF-COMMPASS`)
all_data <- GDCprepare(tumor_normal_project_expQuery$`MMRF-COMMPASS`, 
                       directory = out_dir)

saveRDS(all_data, file = paste0(out_dir, "/all_data_geneCounts_MMRF.rds"))

clinical_data <- all_data@colData %>% as.data.frame()
write_tsv(clinical_data, file = paste0(out_dir, "/clinical_data_all_MMRF.tsv"))

clinical_data_BM <- clinical_data %>% as.data.frame() %>% 
  dplyr::filter(site_of_resection_or_biopsy == "Bone marrow") %>% 
  dplyr::filter(specimen_type == "Bone Marrow NOS")
write_tsv(clinical_data_BM, file = paste0(out_dir, "/clinical_data_BM.tsv"))

clinical_data_BMprimary <- clinical_data %>% as.data.frame() %>% 
  dplyr::filter(site_of_resection_or_biopsy == "Bone marrow") %>% 
  dplyr::filter(specimen_type == "Bone Marrow NOS") %>% 
  dplyr::filter(tumor_descriptor =="Primary")
write_tsv(clinical_data_BMprimary, file = paste0(out_dir, "/clinical_data_BMprimary.tsv")) #764


# survival data
# OS.time: OS.time: Duration from the starting point to the event or censoring.
# OS: Binary indicator of whether the event occurred (1) or was censored (0).
# get ZNF141 group

ZNF141_info <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s4_countGene_R2/res_updateGroupR1/ZNF141_info.tsv")
# patientsIn4Groups <-
# read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/all_activeL1s_patientsIn4Groups.tsv")

# patientsIn4Groups <-
#   read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/all_activeL1PA2_patientsIn4Groups.tsv")

patientsIn4Groups <-
  read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_R2res/all_activeL1PA2_patientsIn4Groups.tsv")


# patientsIn4Groups <-
#   read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s3_countTx_R5/all_activeL1s_patientsIn4Groups.tsv")

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
  dplyr::mutate(OS = ifelse(vital_status=="Alive", 0, 1),
                OStime = days_to_death)


dup_patientIDs <- clinical_data_addIDs[duplicated(clinical_data_addIDs$patient),]$patient

dup_patientIDs_SRRids <- clinical_data_addIDs[duplicated(clinical_data_addIDs$patient),]$SRR_id

clinical_data_addIDs_fil <- clinical_data_addIDs %>% 
  filter(specimen_type == "Bone Marrow NOS")

# git surv
fit <- survfit(Surv(days_to_last_known_disease_status, OS) ~ as.factor(L1_quantile_group), 
               data = clinical_data_addIDs_fil)

# Plot Kaplan-Meier curve
ggsurvplot(fit, data = clinical_data_addIDs_fil, conf.int = TRUE, pval = TRUE, 
           # fun = "pct", 
           test.for.trend=TRUE,
           surv.median.line = "hv")

fit <- survfit(Surv(days_to_last_known_disease_status, OS) ~ ZNF141group, data = clinical_data_addIDs_fil)
# Plot Kaplan-Meier curve
ggsurvplot(fit, data = clinical_data_addIDs_fil, conf.int = TRUE, pval = TRUE)

# filter for primary samples
clinical_data_addIDs_Primary <- clinical_data_addIDs_fil %>% filter(tumor_descriptor =="Primary")
# git surv
fit <- survfit(Surv(days_to_last_known_disease_status, OS) ~ L1_quantile_group, data = clinical_data_addIDs_Primary)
# Plot Kaplan-Meier curve
ggsurvplot(fit, data = clinical_data_addIDs_Primary, conf.int = TRUE, pval = TRUE, surv.median.line = "hv")

fit <- survfit(Surv(days_to_last_known_disease_status, OS) ~ ZNF141group, data = clinical_data_addIDs_Primary)
# Plot Kaplan-Meier curve
ggsurvplot(fit, data = clinical_data_addIDs_Primary, conf.int = TRUE, pval = TRUE, surv.median.line = "hv")

# sink()  # Turn off logging

# to get KZFP expression Quantile groups
dplyr::mutate(
  sampleID = factor(sampleID, levels = sampleID),
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



# get TCGA count mat
## project name, i.e. name for dataset
pro_i <- names(tumor_normal_project_expQuery)
query_i <- tumor_normal_project_expQuery[[pro_i]]
## sample info.
sampInfo_res <- getResults(query_i)

## raw counts and filter
raw_counts <- GDCprepare(query = query_i, 
                         directory = out_dir,
                         summarizedExperiment = FALSE)
raw_counts <- raw_counts[1:60660,] # there are 60660 genes in total

## Extract log2 transformed FPKM values
# cols <- grep("fpkm_unstranded_", colnames(raw_counts))
# dataGBMcomplete <- log2(as.data.frame(raw_counts)[,cols]+1)
## Extract raw counts, not FPKM
cols <- grep("^unstranded_MMRF", colnames(raw_counts))
rawCounts_df <- as.data.frame(raw_counts)[,cols]
## Add rownames and clean the colnames
row.names(rawCounts_df) <- paste(raw_counts$gene_name, raw_counts$gene_id, sep="_")
colnames(rawCounts_df) <- gsub("unstranded_", "", colnames(rawCounts_df))

## create DGElist
### get group info.
if (unique(clinical_data_addIDs$barcode == colnames(rawCounts_df))) {
  groups <- clinical_data_addIDs$L1_quantile_group
  groups[is.na(groups)] <- "nogroup"
}

dge_rawCounts <- DGEList(counts=rawCounts_df, 
                         genes=raw_counts[,1:3], #remove.zeros = TRUE,
                         group = groups)


## normalize counts with TMM method
dge_rawCounts <- calcNormFactors(dge_rawCounts, method = "TMM")
normalized_lcpm <- cpm(dge_rawCounts, normalized.lib.sizes = TRUE, log = TRUE)

boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF141"),] ~groups)

boxplot(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF382"),] ~groups)

plot(as.vector(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF649_"),]) ~clinical_data_addIDs$ZNF141)
plot(as.vector(normalized_lcpm[str_detect(rownames(normalized_lcpm), "ZNF84_"),]) ~clinical_data_addIDs$ZNF141)


