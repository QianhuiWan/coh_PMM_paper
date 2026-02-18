
# R2: update activation score: output from test2: add -O add --countReadPairs
# only keep "activation score" and "TEProf3" and "TElocal"

library(tidyverse)
library(plyranges)
library(rtracklayer)
library(ggplot2)
library(ggpubr)

# load rmsk repeat annotation
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

repeat_L1 <- repeat_GR %>% 
  plyranges::filter(repFamily =="L1") %>% 
  as.data.frame() %>% 
  dplyr::filter(end - (start-1) > 5900)

# 1. get truth TE tx ids
tx_te_info <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/teAnno_round3/useR_plyranges/Gencode_te_transcripts_allAnno_1based.tsv")
tx_te_info_addTxIDs <- tx_te_info %>% 
  dplyr::mutate(transcript_id_version = transcript_id) %>% 
  dplyr::mutate(transcript_id = str_replace(transcript_id_version, "\\.\\d+$", ""))

## read true TEtx for each sim sample and join to tetx annotations
tx_te_info_addTxIDs_addSamplesTruth <- tx_te_info_addTxIDs
for (i in 1:3) {
  sample_filePath <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/results_simulated/truth_counts/*/true_tetx_counts.tsv")
  s_tetx <- read_tsv(sample_filePath[i]) %>% 
    dplyr::mutate(transcript_id_beers2 = transcript_id) %>% 
    dplyr::mutate(transcript_id = str_extract(transcript_id, "ENST\\d+")) %>% 
    dplyr::group_by(transcript_id) %>% 
    dplyr::summarise(count = sum(count)) %>% 
    `colnames<-`(c("transcript_id", 
                   paste0("sample", i)))
  tx_te_info_addTxIDs_addSamplesTruth <- 
    tx_te_info_addTxIDs_addSamplesTruth %>% left_join(s_tetx, by = "transcript_id")
}

tx_te_info_addTxIDs_addSamplesTruth_noNA <- tx_te_info_addTxIDs_addSamplesTruth %>% 
  # keep rows where any of the three columns is not NA
  dplyr::filter(if_any(c(sample1, sample2, sample3), ~ !is.na(.)))



# 2. get salmonTE identified transcribed TE

tx_te_salmonTE <- list()
for (i in 1:3) {
  sample_filePath <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/results_simulated/results/tetx_detect/*/salmonTE/EXPR.csv")
  s_tetx <- read_csv(sample_filePath[i]) %>% 
    `colnames<-`(c("repName", "salmonTE"))
  tx_te_salmonTE[[paste0("salmonTE_sample", i)]] <- s_tetx
}

tx_te_salmonTE <- 
  bind_rows(tx_te_salmonTE, .id = "sample_name")


# 3. get TEProf3 identified transcribed TE
## read in gtf, join with tx_te_info_addTxIDs_addSamplesTruth_noNA

## left_join results with gtf

tx_te_TEProf3 <- list()
for (i in 1:3) {
  sample_gtfPath <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/results_simulated/results/tetx_detect/*/teprof3/teprof3_output/assembled/teprof3_output_TE_transcript_consensus.gtf")
  sample_quantPath <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/results_simulated/results/tetx_detect/*/teprof3/teprof3_output/assembled/teprof3_output_quantification.TE.tsv.gz")
  # load TE tx quant
  s_tetx <- read_tsv(sample_quantPath[i]) 
  # load TE tx gtf, left join with quant
  TEtx_GTFref <- import(sample_gtfPath[i]) %>% 
    as.data.frame(stringsAsFactors = FALSE) %>% 
    dplyr::filter(type == "transcript") %>% 
    dplyr::mutate(repName_geneName = gene_name) %>% 
    dplyr::select(-gene_name) %>% 
    dplyr::left_join(s_tetx, by = "transcript_id") %>% 
    dplyr::mutate(repName = str_remove(repName_geneName, 
                                       paste0("-", gene_name, "$")))

  tx_te_TEProf3[[paste0("TEProf3_sample", i)]] <- TEtx_GTFref
}

tx_te_TEProf3 <- 
  bind_rows(tx_te_TEProf3, .id = "sample_name")


# 4. get activation score identified transcribed TE
tx_te_activationScore <- list()
activeL1_path <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/results_simulated/results/tetx_detect/*/activation_score_addOtest/filtered_L1/*_L1score.tsv")
for (i in 1:length(activeL1_path)) {
  sample <-  str_remove_all(basename(activeL1_path[i]), "_.+")
  # tsv_name_i <- str_replace_all(basename(activeL1_path[i]), "^[^_]+_", "")
  tsv_name_i <- str_remove_all(basename(activeL1_path[i]), ".tsv")
  tsv_file <- read_tsv(file = activeL1_path[i])
  tsv_file$sampleID <- rep(sample, nrow(tsv_file))
  tx_te_activationScore[[tsv_name_i]] <- tsv_file
}

# Function to match classes to the first data frame
match_classes <- function(df, ref_df) {
  for (col in names(ref_df)) {
    if (col %in% names(df)) {
      df[[col]] <- as(df[[col]], class(ref_df[[col]]))
    }
  }
  return(df)
}

# Apply the function to the list
ref_df <- tx_te_activationScore[[1]] # Reference data frame (the first in the list)
tx_te_activationScore <- lapply(tx_te_activationScore, match_classes, ref_df = ref_df)

tx_te_activationScore <- tx_te_activationScore %>% 
  dplyr::bind_rows(.id = "group") %>% 
  dplyr::mutate(direction = str_remove_all(group, "_L1score")) %>% 
  dplyr::mutate(direction = str_remove_all(direction, "..+_"))


# TEtranscripts detected TEs (LINE1)
tx_te_TEtranscripts <- list()
for (i in 1:3) {
  sample_filePath <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/results_simulated/results/tetx_detect/*/TEtranscripts/TEtranscripts_out.cntTable")
  s_tetx <- read_tsv(sample_filePath[i]) %>% 
    dplyr::filter(!str_detect(`gene/TE`, "ENSG")) %>% 
    dplyr::filter(str_detect(`gene/TE`, "LINE")) %>% 
    dplyr::mutate(repName = str_split_fixed(`gene/TE`, ":", 3)[, 1]) %>% 
    `colnames<-`(c("gene/TE", "TEtranscripts", "repName")) %>% 
    dplyr::select(repName, TEtranscripts) %>% 
    distinct(repName, .keep_all = TRUE) # only 1 duplicate
  tx_te_TEtranscripts[[paste0("TEtranscripts_sample", i)]] <- s_tetx
}

tx_te_TEtranscripts <- 
  bind_rows(tx_te_TEtranscripts, .id = "sample_name")


# TElocal detected TEs (LINE1) #################################################
rmskAnn_TElocal <- read_tsv("/home/qwan/nonCondaPkg/TElocal/prebuilt_indices/GRCh38_GENCODE_rmsk_TE.gtf.locInd.locations") %>% 
  mutate(seqnames    = str_extract(`chromsome:start-stop:strand`, "^[^:]+"),
    start  = as.integer(str_match(`chromsome:start-stop:strand`, ":(\\d+)-")[,2]),
    end    = as.integer(str_match(`chromsome:start-stop:strand`, "-(\\d+):")[,2]),
    strand = str_extract(`chromsome:start-stop:strand`, "[+-]$"))

tx_te_TElocal <- list()
for (i in 1:3) {
  sample_filePath <- Sys.glob("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/results_simulated/results/tetx_detect/*/TElocal/TElocal_out.cntTable")
  s_tetx <- read_tsv(sample_filePath[i]) %>% 
    dplyr::filter(!str_detect(`gene/TE`, "ENSG")) %>% 
    dplyr::filter(str_detect(`gene/TE`, "LINE")) %>% 
    dplyr::mutate(rep_id = str_split_fixed(`gene/TE`, ":", 4)[, 1]) %>% 
    dplyr::mutate(repName = str_split_fixed(`gene/TE`, ":", 4)[, 2]) %>% 
    `colnames<-`(c("gene/TE", "TElocal", "rep_id", "repName")) %>% 
    dplyr::left_join(rmskAnn_TElocal, by = c("rep_id" = "TE"))
  tx_te_TElocal[[paste0("TElocal_sample", i)]] <- s_tetx
}

tx_te_TElocal <- 
  bind_rows(tx_te_TElocal, .id = "sample_name")


# 5. calculate F1
## repeat sub-family level
### merge with truth set, sample1
truth_l1 <- tx_te_info_addTxIDs_addSamplesTruth_noNA %>% 
  dplyr::filter(repFamily %in% "L1") %>% 
  dplyr::filter(te_id %in% repeat_L1$te_id)

# get method output for each sample function
getF1 <- function(sample){
  # salmonTE
tx_te_salmonTE_sample <- 
  tx_te_salmonTE[tx_te_salmonTE$sample_name==paste0("salmonTE_",sample),] %>% 
  dplyr::filter(repName %in% repeat_L1$repName)
# TEProf3
tx_te_TEProf3_sample <- 
  tx_te_TEProf3[tx_te_TEProf3$sample_name==paste0("TEProf3_", sample),
                c("sample_name","repName", "stringtie_tpm")] %>% 
  dplyr::group_by(repName) %>% 
  dplyr::summarise(stringtie_tpm = sum(stringtie_tpm))%>% 
  dplyr::filter(repName %in% repeat_L1$repName)
# activation score
tx_te_activationScore_sample <-
  tx_te_activationScore[tx_te_activationScore$sampleID==sample,
                c("sampleID","repName", "Ratio_Sense", "Ratio_Antisense")] %>% 
  dplyr::group_by(repName) %>% 
  dplyr::summarise(Ratio_Sense = sum(Ratio_Sense),
                   Ratio_Antisense = sum(Ratio_Antisense)) %>% 
  dplyr::mutate(Ratio_All = Ratio_Sense+Ratio_Antisense)%>% 
  dplyr::filter(repName %in% repeat_L1$repName)
# TEtranscripts
tx_te_TEtranscripts_sample <- 
  tx_te_TEtranscripts[tx_te_TEtranscripts$sample_name==paste0("TEtranscripts_",sample),] %>% 
  dplyr::filter(repName %in% repeat_L1$repName)

# TElocal
tx_te_TElocal_sample <- 
  tx_te_TElocal[tx_te_TElocal$sample_name == paste0("TElocal_",sample), 
                c("sample_name", "repName", "TElocal")] %>%
  dplyr::group_by(repName) %>% 
  dplyr::summarise(TElocal = sum(TElocal))%>% 
  dplyr::filter(repName %in% repeat_L1$repName)

l1_tp_sample <- tx_te_info_addTxIDs_addSamplesTruth_noNA %>% 
  dplyr::filter(repFamily %in% "L1") %>% 
  dplyr::filter(te_id %in% repeat_L1$te_id) %>% 
  dplyr::mutate(sample = !!sym(sample)) %>% 
  dplyr::group_by(te_id) %>% 
  dplyr::summarise(sample = sum(sample)) %>% 
  dplyr::left_join(repeat_L1, by = "te_id") %>% 
  dplyr::select(repName, sample) %>% 
  dplyr::mutate(sample_truth = ifelse(sample>=1, 1, 0)) %>% 
  dplyr::left_join(tx_te_salmonTE_sample, by = "repName") %>% 
  # dplyr::mutate(sample_salmonTE = ifelse(salmonTE>=1, 1, 0)) %>% 
  dplyr::mutate(sample_salmonTE = ifelse(is.na(salmonTE) | salmonTE < 1, 0, 1)) %>% 
  dplyr::left_join(tx_te_TEProf3_sample, by = "repName") %>% 
  dplyr::mutate(sample_TEProf3 = ifelse(is.na(stringtie_tpm) | stringtie_tpm < 1, 0, 1)) %>% 
  dplyr::left_join(tx_te_activationScore_sample, by = "repName") %>% 
  dplyr::mutate(sample_activationScore = ifelse(is.na(Ratio_All) | Ratio_All < 1, 0, 1)) %>% 
  dplyr::left_join(tx_te_TEtranscripts_sample, by = "repName") %>% 
  dplyr::mutate(sample_TEtranscripts = ifelse(is.na(TEtranscripts) | TEtranscripts < 1, 0, 1)) %>% 
  dplyr::left_join(tx_te_TElocal_sample, by = "repName") %>% 
  dplyr::mutate(sample_TElocal = ifelse(is.na(TElocal) | TElocal < 1, 0, 1)) 
  
f1_sample <- tibble(sample = sample,
                     TP_salmonTE = sum(l1_tp_sample$sample_salmonTE), 
                     FN_salmonTE = sum(l1_tp_sample$sample_salmonTE==0),
                     FP_salmonTE = sum(!(tx_te_salmonTE_sample$repName %in% l1_tp_sample$repName)),
                     precision_salmonTE = TP_salmonTE/(TP_salmonTE+FP_salmonTE),
                     recall_salmonTE =TP_salmonTE/(TP_salmonTE+FN_salmonTE),
                     F1_salmonTE = 2*(precision_salmonTE*recall_salmonTE)/(precision_salmonTE+recall_salmonTE),
                     
                     TP_TEProf3 = sum(l1_tp_sample$sample_TEProf3), 
                     FN_TEProf3 = sum(l1_tp_sample$sample_TEProf3==0),
                     FP_TEProf3 = sum(!(tx_te_TEProf3_sample$repName %in% l1_tp_sample$repName)),
                     precision_TEProf3 = TP_TEProf3/(TP_TEProf3+FP_TEProf3),
                     recall_TEProf3 =TP_TEProf3/(TP_TEProf3+FN_TEProf3),
                     F1_TEProf3 = 2*(precision_TEProf3*recall_TEProf3)/(precision_TEProf3+recall_TEProf3),
                     
                     TP_activationScore = sum(l1_tp_sample$sample_activationScore), 
                     FN_activationScore = sum(l1_tp_sample$sample_activationScore==0),
                     FP_activationScore = sum(!(tx_te_activationScore_sample$repName %in% l1_tp_sample$repName)),
                     precision_activationScore = TP_activationScore/(TP_activationScore+FP_activationScore),
                     recall_activationScore =TP_activationScore/(TP_activationScore+FN_activationScore),
                     F1_activationScore = 2*(precision_activationScore*recall_activationScore)/(precision_activationScore+recall_activationScore),
                    
                    TP_TEtranscripts = sum(l1_tp_sample$sample_TEtranscripts), 
                    FN_TEtranscripts = sum(l1_tp_sample$sample_TEtranscripts==0),
                    FP_TEtranscripts = sum(!(tx_te_TEtranscripts_sample$repName %in% l1_tp_sample$repName)),
                    precision_TEtranscripts = TP_TEtranscripts/(TP_TEtranscripts+FP_TEtranscripts),
                    recall_TEtranscripts = TP_TEtranscripts/(TP_TEtranscripts+FN_TEtranscripts),
                    F1_TEtranscripts = 2*(precision_TEtranscripts*recall_TEtranscripts)/(precision_TEtranscripts+recall_TEtranscripts),
                    
                    TP_TElocal = sum(l1_tp_sample$sample_TElocal), 
                    FN_TElocal = sum(l1_tp_sample$sample_TElocal==0),
                    FP_TElocal = sum(!(tx_te_TElocal_sample$repName %in% l1_tp_sample$repName)),
                    precision_TElocal = TP_TElocal/(TP_TElocal+FP_TElocal),
                    recall_TElocal = TP_TElocal/(TP_TElocal+FN_TElocal),
                    F1_TElocal = 2*(precision_TElocal*recall_TElocal)/(precision_TElocal+recall_TElocal),
                     )
return(f1_sample)
}

F1_all <- rbind(getF1(sample = "sample1"), 
                getF1(sample = "sample2"), 
                getF1(sample = "sample3")) 

F1_all_plot <- F1_all %>% 
  reshape2::melt() %>% 
  dplyr::mutate(score = str_replace_all(variable, "_..+", ""),
                method = str_replace_all(variable, ".+_", "")) %>% 
  dplyr::filter(method %in% c("activationScore", "TEProf3", 
                              "TEtranscripts", "TElocal")) %>%
  ggplot(aes(x=sample, y=value, fill = method))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # sets ymin=0
  facet_wrap(~score, scales = "free_y")+
  geom_bar(stat = "identity", position = "dodge")+
  theme_pubr(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/results_simulated/results/tetx_detect/benchmark_3methods_R4_addTElocal_repNameLevel.pdf",
       plot = F1_all_plot,
       device = "pdf", width = 15, height = 10, units = "cm")



## each repeat loci level only "TElocal", "TEProf3" and "activation score" can have this calculation
## repeat loci level comparison ++++++++++++++++++++++++++++++++++++++++++++++++
### function -- 
getF1_loci <- function(sample){
  # TEProf3
  keys <- c("seqnames", "start", "end", "strand")
  repeat_L1_GR <- repeat_L1 %>% as_granges()
  tx_te_TEProf3_sample <- tx_te_TEProf3 %>%
    filter(sample_name == paste0("TEProf3_", sample)) %>%
    transmute(
          seqnames = as.character(seqnames),
          start    = as.integer(start),
          end      = as.integer(end),
          strand   = as.character(strand),
          loci_tpm = as.numeric(stringtie_tpm)
          ) %>%
    group_by(across(all_of(keys))) %>%
    summarise(loci_tpm = sum(loci_tpm, na.rm = TRUE), .groups = "drop") %>%
    mutate(loci_tpm = replace_na(loci_tpm, 0)) %>% 
    as_granges() %>% 
    find_overlaps(repeat_L1_GR) %>% 
    as.data.frame() %>% 
    group_by(te_id) %>% 
    summarise(loci_tpm = sum(loci_tpm, na.rm = TRUE), .groups = "drop")
  
  # activation score
  tx_te_activationScore_sample <- tx_te_activationScore %>%
    filter(sampleID==sample) %>%
    transmute(
          seqnames = as.character(genoName), start    = as.integer(genoStart),
          end      = as.integer(genoEnd), strand   = as.character(strand),
          Ratio_Sense     = as.numeric(Ratio_Sense),
          Ratio_Antisense = as.numeric(Ratio_Antisense)
          ) %>%
    group_by(across(all_of(keys))) %>%
    summarise(Ratio_Sense = sum(Ratio_Sense, na.rm = TRUE), 
              Ratio_Antisense = sum(Ratio_Antisense, na.rm = TRUE), .groups = "drop") %>% 
    mutate( Ratio_Sense     = replace_na(Ratio_Sense, 0),
          Ratio_Antisense = replace_na(Ratio_Antisense, 0),
          Ratio_All       = Ratio_Sense + Ratio_Antisense) %>% 
    as_granges() %>% 
    find_overlaps(repeat_L1_GR) %>% 
    as.data.frame() %>% 
    group_by(te_id) %>% 
    summarise(Ratio_Sense = sum(Ratio_Sense, na.rm = TRUE), 
              Ratio_Antisense = sum(Ratio_Antisense, na.rm = TRUE), 
              .groups = "drop") %>% 
    mutate(Ratio_All  = Ratio_Sense + Ratio_Antisense)
  
  
  # TElocal
  tx_te_TElocal_sample <- tx_te_TElocal %>%
    filter(sample_name == paste0("TElocal_",sample)) %>%
    transmute(
          seqnames = as.character(seqnames),
          start    = as.integer(start),
          end      = as.integer(end),
          strand   = as.character(strand),
          TElocal = as.numeric(TElocal)
          ) %>%
    group_by(across(all_of(keys))) %>%
    summarise(TElocal = sum(TElocal, na.rm = TRUE), .groups = "drop") %>%
    mutate(TElocal = replace_na(TElocal, 0)) %>% 
    dplyr::filter(str_detect(seqnames, "chr")) %>% 
    as_granges() %>% 
    find_overlaps(repeat_L1_GR) %>% 
    as.data.frame() %>% 
    group_by(te_id) %>% 
    summarise(TElocal = sum(TElocal, na.rm = TRUE), .groups = "drop") 
  
  
  l1_tp_sample <- tx_te_info_addTxIDs_addSamplesTruth_noNA %>% 
    dplyr::filter(repFamily %in% "L1") %>% 
    dplyr::filter(te_id %in% repeat_L1$te_id) %>% 
    dplyr::mutate(sample = !!sym(sample)) %>% 
    dplyr::group_by(te_id) %>% 
    dplyr::summarise(sample = sum(sample)) %>% 
    dplyr::left_join(repeat_L1, by = "te_id") %>% 
    dplyr::select(te_id, sample) %>%
    dplyr::mutate(sample_truth = ifelse(sample>=1, 1, 0)) %>% 
    dplyr::left_join(tx_te_TEProf3_sample, by = "te_id") %>% 
    dplyr::mutate(sample_TEProf3 = ifelse(is.na(loci_tpm) | loci_tpm < 1, 0, 1)) %>% 
    dplyr::left_join(tx_te_activationScore_sample, by = "te_id") %>% 
    dplyr::mutate(sample_activationScore = ifelse(is.na(Ratio_All) | Ratio_All < 1, 0, 1)) %>% 
    dplyr::left_join(tx_te_TElocal_sample, by = "te_id") %>% 
    dplyr::mutate(sample_TElocal = ifelse(is.na(TElocal) | TElocal < 1, 0, 1)) 
  
  f1_sample <- tibble(sample = sample,
                      
                      TP_TEProf3 = sum(l1_tp_sample$sample_TEProf3), 
                      FN_TEProf3 = sum(l1_tp_sample$sample_TEProf3==0),
                      FP_TEProf3 = sum(!(tx_te_TEProf3_sample[tx_te_TEProf3_sample$loci_tpm!=0,]$te_id %in% l1_tp_sample$te_id)),
                      precision_TEProf3 = TP_TEProf3/(TP_TEProf3+FP_TEProf3),
                      recall_TEProf3 =TP_TEProf3/(TP_TEProf3+FN_TEProf3),
                      F1_TEProf3 = 2*(precision_TEProf3*recall_TEProf3)/(precision_TEProf3+recall_TEProf3),
                      
                      TP_activationScore = sum(l1_tp_sample$sample_activationScore), 
                      FN_activationScore = sum(l1_tp_sample$sample_activationScore==0),
                      FP_activationScore = sum(!(tx_te_activationScore_sample[tx_te_activationScore_sample$Ratio_All!=0,]$te_id %in% l1_tp_sample$te_id)),
                      precision_activationScore = TP_activationScore/(TP_activationScore+FP_activationScore),
                      recall_activationScore =TP_activationScore/(TP_activationScore+FN_activationScore),
                      F1_activationScore = 2*(precision_activationScore*recall_activationScore)/(precision_activationScore+recall_activationScore),
                      
                      TP_TElocal = sum(l1_tp_sample$sample_TElocal), 
                      FN_TElocal = sum(l1_tp_sample$sample_TElocal==0),
                      FP_TElocal = sum(!(tx_te_TElocal_sample[tx_te_TElocal_sample$TElocal!=0,]$te_id %in% l1_tp_sample$te_id)),
                      precision_TElocal = TP_TElocal/(TP_TElocal+FP_TElocal),
                      recall_TElocal = TP_TElocal/(TP_TElocal+FN_TElocal),
                      F1_TElocal = 2*(precision_TElocal*recall_TElocal)/(precision_TElocal+recall_TElocal),
  )
  return(f1_sample)
}

F1_loci_all <- rbind(getF1_loci(sample = "sample1"), 
                getF1_loci(sample = "sample2"), 
                getF1_loci(sample = "sample3")) 

F1_loci_all_plot <- F1_loci_all %>% 
  reshape2::melt() %>% 
  dplyr::mutate(score = str_replace_all(variable, "_..+", ""),
                method = str_replace_all(variable, ".+_", "")) %>% 
  # dplyr::filter(method %in% c("activationScore", "TEProf3", 
  #                             "TElocal")) %>%
  ggplot(aes(x=sample, y=value, fill = method))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # sets ymin=0
  facet_wrap(~score, scales = "free_y")+
  geom_bar(stat = "identity", position = "dodge")+
  theme_pubr(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA/RNA_hg38_p14_TEtx_benchmark_R2/results_simulated/results/tetx_detect/benchmark_3methods_R4_addTElocal_repLociLevel.pdf",
       plot = F1_loci_all_plot,
       device = "pdf", width = 15, height = 10, units = "cm")




