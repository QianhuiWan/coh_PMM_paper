# 1. Create a new file `helpers.R` and enter

## Our App
library(tidyverse)
library(readr)
library(dplyr)
library(tibble)
library(magrittr)
# library(here)
# faa <- file.path("test_dat/faa.tsv") %>%
#   read_delim(delim = "\t")  %>%
#   column_to_rownames("FAAs") %>%
#   as.matrix()
# 
# sampleNames <- grep("(TN|HS)", colnames(faa), value = TRUE)
# treatments <- data.frame(
#   HeatShock = as.factor(grepl("HS", sampleNames)),
#   Zinc = as.factor(!grepl("Neg", sampleNames)),
#   row.names = sampleNames
# )
# groupColours <- list(
#   HeatShock = c(`TRUE` = rgb(1, 0, 0), `FALSE` = rgb(0.4, 0, 0)),
#   Zinc = c(`TRUE` = rgb(0.5, 0.5, 0.5), `FALSE` = rgb(0, 0, 0))
# )


mei_ins_only_DNAm_df <- file.path("test_dat/mei_ins_only_DNAm.tsv") %>% read_tsv()
# mei_ins_flank_DNAm_df <- read_tsv(here("ShinyApp_INS/test_dat/mei_ins_flank_DNAm.tsv"))
mei_ins_flank_DNAm_df <- file.path("test_dat/mei_ins_flank_DNAm.tsv") %>% read_tsv()

sample_order <- c("B1", "B2", "B3", 
                  "PMM14", "PMM18", "PMM17", "PMM13", "PMM3", 
                  "PMM12", "PMM16", "PMM15", "PMM11", "PMM4", "PMM6", 
                  "PMM9", "PMM1", "PMM2", "PMM7")

mei_ins_flank_DNAm_df$sample <- factor(mei_ins_flank_DNAm_df$sample, 
                                       levels = rev(sample_order))
mei_ins_only_DNAm_df$sample <- factor(mei_ins_only_DNAm_df$sample,
                                      levels = rev(sample_order))

INS_IDs <- unique(mei_ins_flank_DNAm_df$INS_ID) #n=291

# note: ins_262847 have no ins DNAm



## data for box plot
## get MEI group
mei_group_all <- unique(mei_ins_only_DNAm_df$mei_info)
mei_group_all <- str_remove_all(mei_group_all, "rest_") %>% unique()
mei_group_all <- str_remove(mei_group_all, "_MEI..+")

## get selected mei df
mei_ins_flank_DNAm_df_sel <- mei_ins_flank_DNAm_df %>%
  dplyr::mutate(location = flank_type) %>% 
  dplyr::select(sample:DNAm, location, INS_ID, group)

mei_ins_only_DNAm_df_sel <- mei_ins_only_DNAm_df %>%
  dplyr::mutate(location = coord) %>% 
  dplyr::select(sample:DNAm, location, INS_ID, group)

all_mei_box_df <- rbind(mei_ins_flank_DNAm_df_sel, mei_ins_only_DNAm_df_sel) %>% 
  dplyr::mutate(
    region = ifelse(!location %in% c("flank_up", "flank_down"), "INS", location),
    region = factor(region, levels = c("flank_up", "INS", "flank_down")))






