
# for encode samples/cell lines:
# we want to check KZFP genes expression first, especially ZNF141 and ZNF382

# load packages
library(tidyverse)
library(magrittr)

data_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_encode"
out_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_encode/output"

# H1
H1_tpm <- read_tsv(file = paste0(data_dir, "/H1/H1_RNAExpression.tsv"), 
                   skip = 1)
H1_tpm_KZFP <- H1_tpm %>% 
  dplyr::filter(`Gene symbol` %in% c("ZNF141", "ZNF382")) %>% 
  dplyr::select(ID:`Gene symbol`, `Assay title`, Assembly, 
                `Biosample sex`, `Biosample summary`) %>% 
  dplyr::mutate(cell_line = rep("H1", nrow(.)))

# H9
H9_tpm <- read_tsv(file = paste0(data_dir, "/H9/H9_RNAExpression.tsv"), 
                   skip = 1)
H9_tpm_KZFP <- H9_tpm %>% 
  dplyr::filter(`Gene symbol` %in% c("ZNF141", "ZNF382")) %>% 
  dplyr::select(ID:`Gene symbol`, `Assay title`, Assembly, 
                `Biosample sex`, `Biosample summary`) %>% 
  dplyr::mutate(cell_line = rep("H9", nrow(.)))

# K562
K562_tpm <- read_tsv(file = paste0(data_dir, "/K562/K562_RNAExpression.tsv"), 
                   skip = 1)
K562_tpm_KZFP <- K562_tpm %>% 
  dplyr::filter(`Gene symbol` %in% c("ZNF141", "ZNF382")) %>% 
  dplyr::select(ID:`Gene symbol`, `Assay title`, Assembly, 
                `Biosample sex`, `Biosample summary`) %>% 
  dplyr::mutate(cell_line = 
                  ifelse(str_detect(`Biosample summary`, "treated"), 
                         "IFN treated K562", "K562"))

# HepG2
HepG2_tpm <- read_tsv(file = paste0(data_dir, "/HepG2/HepG2_RNAExpression.tsv"), 
                   skip = 1)
HepG2_tpm_KZFP <- HepG2_tpm %>% 
  dplyr::filter(`Gene symbol` %in% c("ZNF141", "ZNF382")) %>% 
  dplyr::select(ID:`Gene symbol`, `Assay title`, Assembly, 
                `Biosample sex`, `Biosample summary`) %>% 
  dplyr::mutate(cell_line = rep("HepG2", nrow(.)))


# merge TPM for all cell line
all_tpm_KZFP <- rbind(H1_tpm_KZFP, H9_tpm_KZFP, K562_tpm_KZFP, HepG2_tpm_KZFP)

# plot TPM for ZNF141 and ZNF382 across all cell lines: barplot
all_tpm_KZFP %>% 
  ggplot(aes(x = cell_line, y = TPM, fill = `Gene symbol`))+
  geom_bar(position = "dodge") +
  labs(x="", y = "gene expression (TPM)", fill = "Genes")+
  theme_pubr(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p <- all_tpm_KZFP %>%
  group_by(cell_line, `Gene symbol`) %>%
  summarise(mean_TPM = mean(TPM),
            sd_TPM = sd(TPM),
            .groups = "drop") %>%
  ggplot(aes(x = cell_line, y = mean_TPM, fill = `Gene symbol`)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_TPM - sd_TPM, ymax = mean_TPM + sd_TPM),
                position = position_dodge(width = 0.9),
                width = 0.2) +
  labs(x = "", y = "Gene expression (TPM)", fill = "Genes") +
  theme_pubr(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


ggsave(paste0(out_dir, "/encode_cell_line_KZFP_expr.pdf"),
       plot = p,
       device = "pdf", width = 12, height = 10, units = "cm")


