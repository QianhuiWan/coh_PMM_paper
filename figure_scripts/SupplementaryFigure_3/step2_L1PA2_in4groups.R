
# density plot and line/dot plot
# use s6_activationScore_R3 output
# devide samples in 4 groups based on active L1 subfamily/L1PA2 numbers

input_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_redo"
output_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/s6_activationScore_R3/filtered_L1_redo_res"

MM_metadata_update <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/clinical_data_BMprimary_addID_addSubtypes.tsv")


# .libPaths("/home/qwan/miniconda3/envs/coh/lib/R/library")
library(tidyverse)
library(data.table)
library(dplyr)
library(readr)
library(rlang)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(patchwork)
options(scipen=999)

# read in L1 tsv files for all MM samples
# read in tsv files
## get path of tsv files
activeL1_path <- list.files(input_dir, recursive = TRUE, 
                            pattern = paste0(".*\\.tsv$"), 
                            full.names = TRUE)

activeL1_list <- list()
for (i in 1:length(activeL1_path)) {
  sample <-  str_remove_all(basename(activeL1_path[i]), "_.+")
  # tsv_name_i <- str_replace_all(basename(activeL1_path[i]), "^[^_]+_", "")
  tsv_name_i <- str_remove_all(basename(activeL1_path[i]), ".tsv")
  tsv_file <- read_tsv(file = activeL1_path[i])
  tsv_file$sampleID <- rep(sample, nrow(tsv_file))
  activeL1_list[[tsv_name_i]] <- tsv_file
}

# Function to match classes for all dfs in the list to the first data frame
match_classes <- function(df, ref_df) {
  for (col in names(ref_df)) {
    if (col %in% names(df)) {
      df[[col]] <- as(df[[col]], class(ref_df[[col]]))
    }
  }
  return(df)
}

# Apply the function to the list
ref_df <- activeL1_list[[1]] # Reference data frame (the first in the list)
activeL1_list <- lapply(activeL1_list, match_classes, ref_df = ref_df)

activeL1_combined_df <- activeL1_list %>% 
  dplyr::bind_rows(.id = "group") 

# get default situation (named sense) L1 activation
activeL1_combined_sense <- activeL1_combined_df %>% 
  dplyr::filter(str_detect(group, "predict_sense"))

write_tsv(activeL1_combined_sense, 
          file = paste0(output_dir, "/pred_sense_L1s_MM.tsv"))

# get predicted antisense L1 activation
activeL1_combined_antisense <- activeL1_combined_df %>% 
  dplyr::filter(str_detect(group, "predict_antisense"))

write_tsv(activeL1_combined_antisense, 
          file = paste0(output_dir, "/pred_antisense_L1s_MM.tsv"))

# # get default situation (named sense) L1 activation
# activeL1_combined_sense <-
#   read_tsv(file = paste0(input_dir, "/pred_sense_L1s_step2R2_S6R3.tsv"))
# # get predicted antisense L1 activation
# activeL1_combined_antisense <-
#   read_tsv(file = paste0(input_dir, "/pred_antisense_L1s_step2R2_S6R3.tsv"))


# plot active pred. sense L1 distribution
pred_senseL1_num <- 
  activeL1_combined_sense %>% 
  group_by(sampleID) %>% 
  summarise(L1_num = dplyr::n(), .groups = 'drop') #%>% 
  # add this sample with 0 sense L1 in case 
  # bind_rows(data.frame(sampleID = "SRR1606171", L1_num = 0)) 

pred_senseL1subfamily_num <- 
  activeL1_combined_sense %>% 
  group_by(sampleID, repName) %>% 
  summarise(L1subfamily_num = dplyr::n(), .groups = 'drop') %>% 
  full_join(pred_senseL1_num, by = "sampleID") %>% 
  replace_na(list(
    repName = "L1_null",
    L1subfamily_num = 0,
    L1_num = 0 )) %>% 
  mutate(percentage = ifelse(L1_num > 0, L1subfamily_num / L1_num * 100, 0),
         # Reorder by sense L1_num
         sampleID = factor(sampleID, levels = unique(sampleID[order(L1_num)])),
         # Create a numeric index
         sample_index_senseL1 = as.numeric(factor(sampleID)))   

write_tsv(pred_senseL1subfamily_num, 
          file = paste0(output_dir, "/pred_senseL1_subFamily_num.tsv"))


# sense L1 subfamily plots
p1 <- 
  pred_senseL1subfamily_num %>% 
  ggplot(aes(x=L1subfamily_num))+
  geom_bar()+
  labs(x = "pred. active sense L1 subfamily number", y="patients number")+
  facet_wrap(~repName, scales = "free_y")+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(paste0(output_dir, "/pred_senseL1_subFamily_num_barPlot.pdf"), plot = p1,
       width = 12, height = 9, units = "in", limitsize = FALSE)


for (bar_position in c("fill","stack")) {
  if (bar_position == "fill") {
    ylab = "Percentage of L1 subfamily"
  }else{
    ylab = "Number of activated L1 subfamily"
  }
  p2 <- 
    pred_senseL1subfamily_num %>% 
    dplyr::mutate(
      group = ifelse(
        repName %in% c("L1HS", "L1PA2", "L1PA3", "L1PA4", 
                       "L1PA5", "L1PA6", "L1PA7", "L1PA8"), 
        repName, "others")) %>% 
    ggplot(aes(x = sampleID, y = L1subfamily_num, fill = group))+
    # geom_bar(stat = "identity", position = "stack")+
    geom_bar(stat = "identity", position = bar_position)+
    # scale_fill_brewer(palette = "Dark2") +
    scale_fill_viridis_d(option = "turbo") + #, direction = -1
    labs(   x = "Samples (reordered by increasing sense L1 number)",
            y = ylab,
            fill = "L1 subfamily")+
    theme_pubr()+
    theme(
      axis.text.x = element_blank(),    # Remove x-axis text labels
      axis.ticks.x = element_blank()   # Remove x-axis tick marks
    )
  
  ggsave(paste0(output_dir, 
                "/pred_senseL1_subFamily_num_barPlot_", bar_position,".pdf"), 
         plot = p2,
         width = 5.7, height = 3.8, units = "in", limitsize = FALSE)
}





## antisense (just based on prediction, not on stranded library)
pred_antisenseL1_num <- 
  activeL1_combined_antisense %>% 
  group_by(sampleID) %>% 
  summarise(L1_num = dplyr::n(), .groups = 'drop')


pred_antisenseL1subfamily_num <- 
  activeL1_combined_antisense %>% 
  group_by(sampleID, repName) %>% 
  summarise(L1subfamily_num =  dplyr::n(), .groups = 'drop') %>% 
  full_join(pred_antisenseL1_num, by = "sampleID") %>% 
  replace_na(list(
    repName = "L1_null",
    L1subfamily_num = 0,
    L1_num = 0 )) %>% 
  mutate(percentage = ifelse(L1_num > 0, L1subfamily_num / L1_num * 100, 0),
         # Reorder by antisense L1_num
         sampleID = 
           factor(sampleID, levels = unique(sampleID[order(L1_num)])),
         # Create a numeric index
         sample_index_antisenseL1 = as.numeric(factor(sampleID))
         )

# antisense L1 subfamily plots
p3 <- 
  pred_antisenseL1subfamily_num %>% 
  ggplot(aes(x=L1subfamily_num))+
  geom_bar()+
  labs(x = "pred. active antisense L1 subfamily number", y="patients number")+
  facet_wrap(~repName, scales = "free_y")+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(paste0(output_dir, "/pred_antisenseL1_subFamily_num_barPlot.pdf"), 
       plot = p3,
       width = 12, height = 9, units = "in", limitsize = FALSE)

## samples ordered by increasing antisense L1 number
for (bar_position in c("fill","stack")) {
  if (bar_position == "fill") {
    ylab = "Percentage of L1 subfamily"
  }else{
    ylab = "Number of activated L1 subfamily"
  }
  p4 <- 
    pred_antisenseL1subfamily_num %>% 
    dplyr::mutate(
      group = ifelse(
        repName %in% c("L1HS", "L1PA2", "L1PA3", "L1PA4", 
                       "L1PA5", "L1PA6", "L1PA7", "L1PA8"), 
        repName, "others")) %>% 
    ggplot(aes(x = sampleID, y = L1subfamily_num, fill = group))+
    # ggplot(aes(x = sampleID_sense, y = L1subfamily_num, fill = group))+
    geom_bar(stat = "identity", position = bar_position)+
    # scale_fill_brewer(palette = "Dark2") +
    scale_fill_viridis_d(option = "turbo") + #, direction = -1
    labs(   x = "Samples (reordered by increasing antisense L1 number)",
            y = ylab,
            fill = "L1 subfamily")+
    theme_pubr()+
    theme(
      axis.text.x = element_blank(),    # Remove x-axis text labels
      axis.ticks.x = element_blank()   # Remove x-axis tick marks
    )
  
  ggsave(paste0(output_dir, 
                "/pred_antisenseL1_subFamily_num_barPlot_",bar_position,".pdf"), 
         plot = p4,
         width = 5.7, height = 3.8, units = "in", limitsize = FALSE)
}


## samples ordered by increasing sense L1 number
for (bar_position in c("fill","stack")) {
  if (bar_position == "fill") {
    ylab = "Percentage of L1 subfamily"
  }else{
    ylab = "Number of activated L1 subfamily"
  }
  p4 <- 
    pred_antisenseL1subfamily_num %>% 
    dplyr::mutate(
      group = ifelse(
        repName %in% c("L1HS", "L1PA2", "L1PA3", "L1PA4", 
                       "L1PA5", "L1PA6", "L1PA7", "L1PA8"), 
        repName, "others")) %>% 
    # ggplot(aes(x = sampleID, y = L1subfamily_num, fill = group))+
    ggplot(aes(x = factor(sampleID, levels =levels(pred_senseL1subfamily_num$sampleID)), 
                          y = L1subfamily_num, fill = group))+
    geom_bar(stat = "identity", position = bar_position)+
    # scale_fill_brewer(palette = "Dark2") +
    scale_fill_viridis_d(option = "turbo") + #, direction = -1
    labs(   x = "Samples (reordered by increasing sense L1 number)",
            y = ylab,
            fill = "L1 subfamily")+
    theme_pubr()+
    theme(
      axis.text.x = element_blank(),    # Remove x-axis text labels
      axis.ticks.x = element_blank()   # Remove x-axis tick marks
    )
  
  ggsave(paste0(output_dir, 
                "/pred_antisenseL1_subFamily_num_barPlot_",bar_position,
                "_orderBySenseSampID.pdf"), 
         plot = p4,
         width = 5.7, height = 3.8, units = "in", limitsize = FALSE)
}




# sense and antisense together
all_activeL1s <- 
  bind_rows(pred_senseL1subfamily_num[1:5], pred_antisenseL1subfamily_num[1:5]) %>% 
  group_by(sampleID) %>% 
  summarise(L1_num = sum(L1subfamily_num), .groups = 'drop') %>% 
  left_join(unique(pred_senseL1subfamily_num[, c(1,6)]), by = "sampleID") %>% 
  left_join(unique(pred_antisenseL1subfamily_num[, c(1,6)]), by = "sampleID") 

all_activeL1subfamily <- 
  bind_rows(pred_senseL1subfamily_num, pred_antisenseL1subfamily_num) %>% 
  group_by(sampleID, repName) %>% 
  summarise(L1subfamily_num = sum(L1subfamily_num), .groups = 'drop') %>% 
  full_join(all_activeL1s, by = "sampleID") %>% 
  dplyr::filter(!repName == "L1_null") %>% 
  mutate(percentage = ifelse(L1_num > 0, L1subfamily_num / L1_num * 100, 0),
         # Reorder levels by all L1_num
         sampleID = 
           factor(sampleID, levels = unique(sampleID[order(L1_num)])),
         # Create a numeric index
         sample_index_allL1 = as.numeric(factor(sampleID)))


p5 <- 
  all_activeL1s %>% 
  ggplot(aes(x=L1_num))+
  # geom_bar(width = 1)+
  geom_density()+
  labs(x = "all active L1 number", y="density")+
  theme_pubr()

ggsave(paste0(output_dir, "/all_activeL1_num_densityPlot.pdf"), plot = p5,
       width = 4.5, height = 3, units = "in", limitsize = FALSE)


## samples ordered by all L1 number
plot_fun <- function(all_activeL1subfamily, order){
  if (order == "local") {
    all_activeL1subfamily$sampleID = all_activeL1subfamily$sampleID
  }else if (order == "sense") {
    all_activeL1subfamily$sampleID = 
      factor(all_activeL1subfamily$sampleID, 
             levels =levels(pred_senseL1subfamily_num$sampleID))
  }else if (order == "antisense") {
    all_activeL1subfamily$sampleID = 
      factor(all_activeL1subfamily$sampleID, 
             levels =levels(pred_antisenseL1subfamily_num$sampleID))
  }
  
  for (bar_position in c("fill","stack")) {
    if (bar_position == "fill") {
      ylab = "Percentage of L1 subfamily"
    }else{
      ylab = "Number of activated L1 subfamily"
    }
    p6 <- 
      all_activeL1subfamily %>% 
      dplyr::mutate(
        group = ifelse(
          repName %in% c("L1HS", "L1PA2", "L1PA3", "L1PA4", 
                         "L1PA5", "L1PA6", "L1PA7", "L1PA8"), 
          repName, "others")) %>% 
      ggplot(aes(x = sampleID, y = L1subfamily_num, fill = group))+
      geom_bar(stat = "identity", position = bar_position)+
      # scale_fill_brewer(palette = "Dark2") +
      scale_fill_viridis_d(option = "turbo") + #, direction = -1
      labs(   x = "Samples (reordered by increasing L1 number)",
              y = ylab,
              fill = "L1 subfamily")+
      theme_pubr()+
      theme(
        axis.text.x = element_blank(),    # Remove x-axis text labels
        axis.ticks.x = element_blank()   # Remove x-axis tick marks
      )
    
    ggsave(paste0(output_dir, 
                  "/all_L1_subFamily_num_barPlot_",bar_position,
                  "_orderBy", order, "SampID.pdf"), 
           plot = p6,
           width = 5.7, height = 3.8, units = "in", limitsize = FALSE)
  }
}

plot_fun(all_activeL1subfamily, order = "local")
plot_fun(all_activeL1subfamily, order = "sense")
plot_fun(all_activeL1subfamily, order = "antisense")





# Patients in 4 groups: quantile by sample index
# MM_subtypes <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_MMRF/mmrf_subtypes/clinical_data_BMprimary_addID_addSubtypes.tsv")


# MM_subtypes_fil <- MM_subtypes %>% 
#   dplyr::filter(Reason_For_Collection == "Baseline")

all_activeL1s_patientsIn4Groups <- 
  all_activeL1s %>% 
  dplyr::filter(sampleID %in% MM_metadata_update$SRR_id) %>% 
  arrange(L1_num) %>% 
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
    L1_num,
    breaks = quantile(L1_num, probs = seq(0, 1, 0.25), na.rm = TRUE),
    include.lowest = TRUE,
    labels = c("Q1", "Q2", "Q3", "Q4")
    ))

write_tsv(all_activeL1s_patientsIn4Groups, 
          file = paste0(output_dir, "/all_activeL1s_patientsIn4Groups_764.tsv"))


all_activeL1subfamily_patientsIn4Groups <- 
  all_activeL1subfamily %>%
  dplyr::select(-sample_index_senseL1, -sample_index_antisenseL1) %>% 
  dplyr::filter(sampleID %in% MM_metadata_update$SRR_id) %>% 
  dplyr::select(-L1_num) %>% 
  dplyr::left_join(all_activeL1s_patientsIn4Groups, by = "sampleID")

write_tsv(all_activeL1subfamily_patientsIn4Groups, 
          file = paste0(output_dir, "/all_activeL1_subfamily_patientsIn4Groups_764.tsv"))

all_activeL1PA2_patientsIn4Groups <- 
  all_activeL1subfamily %>% 
  dplyr::filter(sampleID %in% MM_metadata_update$SRR_id) %>% 
  pivot_wider(id_cols = c(sampleID, sample_index_antisenseL1), 
              names_from = repName, 
              values_from = L1subfamily_num, values_fill = 0) %>% 
  dplyr::select(sampleID, sample_index_antisenseL1, L1PA2) %>%  
  arrange(L1PA2) %>% 
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
write_tsv(all_activeL1PA2_patientsIn4Groups, 
          file = paste0(output_dir, "/all_activeL1PA2_patientsIn4Groups_764.tsv"))


# Patients in 4 groups: quantile by active L1 number
# Compute the L1 quantile positions (y-values)
L1_quantile_positions <- floor(quantile(all_activeL1s_patientsIn4Groups$L1_num, 
                                        probs = seq(0, 1, 0.25), na.rm = TRUE))

# Find the exact x-axis positions corresponding to these quantiles
x_positions <- sapply(L1_quantile_positions, function(y) {
  # Get the indices of rows where L1_num matches the quantile value
  indices <- which(all_activeL1s_patientsIn4Groups$L1_num == y)
  
  # If there are multiple matches, take the first one; otherwise return NA
  if (length(indices) > 0) {
    indices[length(indices)]  # Take the last occurrence
  } else {
    NA
  }
})




p7 <- 
  all_activeL1s_patientsIn4Groups %>% 
  dplyr::filter(sampleID %in% MM_metadata_update$SRR_id) %>% 
  ggplot(aes(x = sample_index, y = L1_num, color = L1_quantile_group)) +
  geom_point() +
  scale_x_continuous(
    breaks = x_positions,  # Quantile positions
    labels = x_positions   # Use quantile values as labels
  ) +
  labs(x = "Sample Index", y = "Active L1 Number", 
       color = "Quantile Group") +
  theme_pubr() +
  theme(
    axis.ticks.x = element_line()
  )

ggsave(paste0(output_dir, 
              "/all_activeL1s_patientsIn4Groups_useL1NumQuantile_dotPlot_764.pdf"), 
       plot = p7,
       width = 4.5, height = 3, units = "in", limitsize = FALSE)


p8 <- 
  all_activeL1subfamily_patientsIn4Groups %>% 
  ggplot(aes(x = sample_index, y = L1subfamily_num, color = L1_quantile_group)) +
  geom_point() +
  scale_x_continuous(
    breaks = x_positions,  # Quantile positions
    labels = x_positions   # Use quantile values as labels
  ) +
  facet_wrap(~repName)+
  labs(x = "Sample Index", y = "Active L1 Number", 
       color = "Quantile Group") +
  theme_pubr() +
  theme(
    axis.ticks.x = element_line()
  )

ggsave(paste0(output_dir, 
              "/all_activeL1subfamily_patientsIn4Groups_useL1NumQuantile_dotPlot_764.pdf"), 
       plot = p8,
       width = 12, height = 9, units = "in", limitsize = FALSE)

