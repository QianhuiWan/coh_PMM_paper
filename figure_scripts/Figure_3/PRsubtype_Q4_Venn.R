
library(tidyverse)
library(readxl)
library(VennDiagram)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

output_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/mmrf_subtypes"

# load MMRF meta
MM_metadata <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/clinical_data_BMprimary_addID.tsv")

# load RNA sub-types
MM_RNAsubtype <- read_excel(path = paste0(output_dir, "/mmrf_RNA_subtypes.xlsx"))

MM_RNAsubtype_baseline <- MM_RNAsubtype %>% 
  dplyr::filter(Reason_For_Collection == "Baseline")
  
# Venn plot

table(MM_metadata$submitter_id %in% MM_RNAsubtype_baseline$Patient_ID) # 709 in common

MM_metadata_update <- MM_metadata %>% 
  dplyr::left_join(MM_RNAsubtype_baseline, by = c("submitter_id"="Patient_ID"))

write_tsv(MM_metadata_update, 
          file = paste0(output_dir, "/clinical_data_BMprimary_addID_addSubtypes.tsv"))

MM_metadata_update709 <- MM_metadata_update %>% 
  dplyr::filter(submitter_id %in% MM_RNAsubtype_baseline$Patient_ID)

MM_metadata_update_L1PA2top5 <- MM_metadata_update %>% 
  dplyr::filter(submitter_id %in% MM_RNAsubtype_baseline$Patient_ID) %>% 
  dplyr::arrange(desc(L1PA2)) %>% 
  # dplyr::slice(1:76)
  dplyr::filter(L1PA2 >= quantile(L1PA2, probs = 0.95))

table(MM_metadata_update709$RNA_Subtype_Name) %>% sort() # 51 PR
nrow(MM_metadata_update_L1PA2top5) # 36 in top 5%
table(MM_metadata_update_L1PA2top5$RNA_Subtype_Name) %>% sort() # 22 PR


Venn <- 
  draw.pairwise.venn(
    area1 = 51,
    area2 = 36,
    cross.area = 22,
    category = c("Patients\nin PR stubtype", 
                 "Top 5% patients ranked by\nthe number of active L1PA2s"),
    fill = brewer.pal(9, "Set1")[c(2, 1)],
    scaled = TRUE,
    alpha = 0.6,
    cex = 1.2, cat.cex = 0.7, 
    cat.fontface = 2,
    cat.fontfamily = rep("sans", 2),
    cat.just =list(c(0, -1), c(1, -1.5)),
    lty = "blank",
    cat.dist = c(0.05, 0.05)) %>% 
  as_ggplot()

ggsave(filename = paste0(output_dir, "/mmrf_RNA_subtype_PR_L1PA2high_Venn.pdf"),
       plot = Venn,
       device = "pdf", width = 9.5, 
       height = 9, units = "cm")




# Input values
PR <- 51
L1PA2 <- 36
Overlap <- 22
N_total <- 709

# Fill in the 2x2 contingency table
fisher_matrix <- matrix(c(
  Overlap,             # in PR and L1PA2
  PR - Overlap,        # in PR, not in L1PA2
  L1PA2 - Overlap,     # in L1PA2, not in PR
  N_total - PR - L1PA2 + Overlap  # in neither
), nrow = 2)

colnames(fisher_matrix) <- c("L1PA2", "Not_L1PA2")
rownames(fisher_matrix) <- c("PR", "Not_PR")

# Run Fisher's Exact Test
fisher.test(fisher_matrix)
pval <- fisher.test(fisher_matrix)$p.value

expected_overlap <- (PR * L1PA2) / N_total
enrichment_score <- Overlap / expected_overlap


# add p val to Venn

pdf(file =paste0(output_dir, "/mmrf_RNA_subtype_PR_L1PA2high_Venn2.pdf"), 
    width = 9.5/2.54, height = 9/2.54)
Venn
grid.text(
  label = paste0("Expected = ", round(expected_overlap, 1),
                 "\nFold Enr. = ", round(enrichment_score, 2),
                 "\nP = ", format.pval(pval, digits = 2, scientific = TRUE)),
  x = unit(0.3, "npc"), y = unit(0.8, "npc"), just = "left", gp = gpar(fontsize = 10)
)

dev.off()


