# ============================================================
# Q1–Q4 vs metadata association test table (incl. MM subtype)
# - Continuous vars: Kruskal–Wallis (and optional pairwise Wilcoxon)
# - Categorical vars: Chi-squared (simulate p if sparse) / Fisher for 2x2
# - Outputs: CSV tables (ALL samples and MM-only), optional gtsummary Table-1
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(forcats)
    library(ggplot2)
    library(ggpubr)
})

# ----------------------------
# 0) LOAD DATA & USER INPUTS
# ----------------------------

# Path to the clinical data
# meta_file <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/clinical_data_BMprimary_addID_addSubtypes_uniqueMap.tsv"
# meta_file <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/clinical_data_BMprimary_addID_uniqueMap.tsv"
meta_file <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/MMRF_openAccess/clinical_data_BMprimary_addID_addL1_addSubtypes_uniqueMap.tsv"

if (!file.exists(meta_file)) {
    stop("Metadata file not found: ", meta_file)
}

# Load the data
meta <- readr::read_tsv(meta_file, show_col_types = FALSE)

# Definitions
sample_id_col <- "SRR_id" # or "barcode" depending on what you want to track
group_col <- "L1_quantile_group" # The pre-calculated Q1-Q4 column

# Variables to test
vars_to_test <- c(
    "age_at_index",
    "gender",
    "iss_stage",
    # "ecog_performance_status",
    "RNA_Subtype_Name", # or "subtype" depending on column name
    "Low purity",
    # "L1PA2"
    "L1_num"
)

# Output directory / file prefix
out_prefix <- "Q1Q4_metadata_assoc"

# Optional: lump rare subtype levels (recommended if many rare categories)
lump_subtype_levels <- TRUE
lump_subtype_col <- "RNA_Subtype_Name"
lump_n_keep <- 12

# Optional: pairwise Wilcoxon for continuous variables
do_pairwise_for_continuous <- TRUE
pairwise_p_adjust <- "BH"
pairwise_alpha <- 0.05

# ----------------------------
# 1) PREPARE DATAFRAME `df`
# ----------------------------

# Ensure the grouping column exists
if (!group_col %in% colnames(meta)) {
    stop("Grouping column '", group_col, "' not found in metadata.")
}

df <- meta %>%
    dplyr::filter(!is.na(.data[[group_col]])) %>%
    dplyr::mutate(
        Q = factor(.data[[group_col]], levels = c("Q1", "Q2", "Q3", "Q4")) # Ensure correct order
    )

# Optional subtype lumping
if (lump_subtype_col %in% colnames(df) && lump_subtype_levels) {
    df <- df %>%
        mutate(
            "{lump_subtype_col}" := forcats::fct_lump_n(.data[[lump_subtype_col]],
                n = lump_n_keep,
                other_level = "Other"
            )
        )
}

# Coerce specific columns if they are character but contain numbers (like ISS_stage sometimes)
# Note: ISS_stage is often I, II, III (roman), so as.numeric might fail if not converted.
# If ISS_stage is 1, 2, 3 it's fine. If "I", "II", "III", keep as factor.


# ----------------------------
# 2) TEST FUNCTIONS
# ----------------------------

# Categorical test helper:
# - If 2x2 and sparse -> Fisher
# - Else Chi-square; if sparse expected -> chi-square with simulated p
cat_test_p <- function(tab) {
    # Remove empty levels
    tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
    if (nrow(tab) < 2 || ncol(tab) < 2) {
        return(list(p = NA_real_, test = "NA (<2 levels)"))
    }

    # expected counts
    chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
    sparse <- any(chi$expected < 5)

    if (all(dim(tab) == c(2, 2)) && sparse) {
        p <- suppressWarnings(fisher.test(tab)$p.value)
        return(list(p = p, test = "Fisher's exact (2x2)"))
    }

    if (sparse) {
        p <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 20000)$p.value)
        return(list(p = p, test = "Chi-squared (simulated p)"))
    } else {
        p <- suppressWarnings(chisq.test(tab, correct = FALSE)$p.value)
        return(list(p = p, test = "Chi-squared"))
    }
}

# Main association function: one variable vs Q
assoc_one <- function(dat, group_col, var) {
    stopifnot(group_col %in% colnames(dat), var %in% colnames(dat))

    g <- dat[[group_col]]
    # Ensure g has all levels Q1-Q4
    if (!is.factor(g)) g <- factor(g, levels = c("Q1", "Q2", "Q3", "Q4"))

    x <- dat[[var]]

    keep <- !is.na(g) & !is.na(x)
    g_sub <- g[keep]
    x_sub <- x[keep]

    # Initialize empty result structure

    if (length(unique(g_sub)) < 2) {
        return(tibble(
            Variable = var, Level = NA_character_, Enriched_In = NA_character_, Test = NA_character_, P_value = NA_real_,
            Q1 = NA_character_, Q2 = NA_character_, Q3 = NA_character_, Q4 = NA_character_
        ))
    }

    # Continuous (numeric)
    if (is.numeric(x_sub)) {
        p_val <- suppressWarnings(kruskal.test(x_sub ~ g_sub)$p.value)

        # per-group summary
        summ_grouped <- tibble(g = g_sub, x = x_sub) %>%
            group_by(g) %>%
            dplyr::summarise(
                n = dplyr::n(),
                median = median(x, na.rm = TRUE),
                IQR = IQR(x, na.rm = TRUE),
                .groups = "drop"
            )

        # Identify enrichment (Group with max median)
        max_idx <- which.max(summ_grouped$median)
        enriched_group <- as.character(summ_grouped$g[max_idx])

        # Pivot for table
        summ <- summ_grouped %>%
            mutate(summary = paste0(round(median, 2), " (", round(IQR, 2), ")")) %>%
            select(g, summary) %>%
            # Ensure all levels present
            tidyr::complete(
                g = factor(c("Q1", "Q2", "Q3", "Q4"), levels = c("Q1", "Q2", "Q3", "Q4")),
                fill = list(summary = "NA")
            ) %>%
            tidyr::pivot_wider(names_from = g, values_from = summary)

        return(tibble(
            Variable = var,
            Level = "Median (IQR)",
            Enriched_In = enriched_group,
            Test = "Kruskal-Wallis",
            P_value = p_val
        ) %>% bind_cols(summ))
    }

    # Categorical
    # Table with all levels
    tab <- table(g_sub, x_sub)
    ct <- cat_test_p(tab)

    # Calculate n and %
    # tab is Group x Level.
    # prop.table(tab, 1) gives row proportions (sum(row) = 1)
    props <- prop.table(tab, 1) * 100

    # Enrichment: For a given Level (Column of tab), which Group (Row) has the highest %?
    # We need to iterate over Levels (Columns)
    enrich_vec <- character(ncol(props))
    for (j in seq_len(ncol(props))) {
        # props is (Groups x Levels)
        # Identify which Group (row) has max value for this Level (col)
        max_row_idx <- which.max(props[, j])
        max_group <- rownames(props)[max_row_idx]
        max_val <- round(props[max_row_idx, j], 1)
        enrich_vec[j] <- paste0(max_group, " (", max_val, "%)")
    }

    # Format: n (p%)
    mat_str <- matrix(paste0(tab, " (", round(props, 1), "%)"),
        nrow = nrow(tab), dimnames = dimnames(tab)
    )

    # Convert to df: Level | Q1 | Q2 | Q3 | Q4
    # Transpose so rows are Levels, cols are Groups
    df_cat <- as.data.frame.matrix(t(mat_str)) %>%
        tibble::rownames_to_column("Level")

    # Add metadata
    df_cat <- df_cat %>%
        mutate(
            Variable = var,
            Enriched_In = enrich_vec, # This aligns because enrich_vec was built iterating cols (Levels)
            Test = ct$test,
            P_value = ct$p
        ) %>%
        select(Variable, Level, Enriched_In, Test, P_value, everything())

    return(df_cat)
}

# Generate full association table
make_assoc_table <- function(dat, group_col, vars) {
    # Process each variable
    res_list <- lapply(vars, function(v) assoc_one(dat, group_col, v))

    # Bind rows
    res <- bind_rows(res_list) %>%
        mutate(
            FDR = p.adjust(P_value, method = "BH")
        ) %>%
        select(Variable, Level, Enriched_In, Test, P_value, FDR, everything())

    res
}

# Optional pairwise Wilcoxon for continuous vars (if global KW significant)
pairwise_continuous <- function(dat, group_col, var, p_adjust = "BH") {
    g <- dat[[group_col]]
    x <- dat[[var]]
    keep <- !is.na(g) & !is.na(x)
    g <- droplevels(g[keep])
    x <- x[keep]
    if (!is.numeric(x) || length(unique(g)) < 2) {
        return(NULL)
    }

    pw <- pairwise.wilcox.test(x, g, p.adjust.method = p_adjust)
    # convert to long tidy
    m <- pw$p.value
    if (is.null(m)) {
        return(NULL)
    }

    # Reshape matrix to long format safely
    # as.data.frame.table always returns Var1, Var2, Freq by default
    res <- as.data.frame(as.table(m))
    colnames(res) <- c("group1", "group2", "p")

    res %>%
        as_tibble() %>%
        filter(!is.na(p)) %>%
        mutate(variable = var, p_adj_method = p_adjust) %>%
        arrange(p)
}

# ----------------------------
# 3) RUN: ALL SAMPLES TABLE
# ----------------------------
vars_exist <- vars_to_test[vars_to_test %in% colnames(df)]
if (length(vars_exist) == 0) stop("None of vars_to_test exist in df. Check column names.")

res_all <- make_assoc_table(df, "Q", vars_exist)
print(res_all)

write.csv(res_all, paste0(out_prefix, "_ALL_L1.csv"), row.names = FALSE)

# Optional pairwise for continuous vars
if (do_pairwise_for_continuous) {
    cont_vars <- vars_exist[sapply(df[vars_exist], is.numeric)]
    # only do pairwise if global KW p < alpha
    sig_cont <- res_all %>%
        filter(Test == "Kruskal-Wallis", Variable %in% cont_vars, !is.na(P_value), P_value < pairwise_alpha) %>%
        pull(Variable)

    pw_list <- lapply(sig_cont, function(v) pairwise_continuous(df, "Q", v, pairwise_p_adjust))
    pw_all <- bind_rows(pw_list)
    if (nrow(pw_all) > 0) {
        write.csv(pw_all, paste0(out_prefix, "_ALL_pairwiseWilcox.csv"), row.names = FALSE)
    }
}

# # ----------------------------
# # 4) RUN: MM-ONLY TABLE (incl. subtype)
# # ----------------------------
# if (!is.null(disease_col) && disease_col %in% colnames(df)) {
#     df_mm <- df %>%
#         filter(.data[[disease_col]] == mm_value) %>%
#         droplevels()
#
#     if (nrow(df_mm) > 0) {
#         vars_mm <- vars_to_test[vars_to_test %in% colnames(df_mm)]
#         res_mm <- make_assoc_table(df_mm, "Q", vars_mm)
#         print(res_mm)
#
#         write.csv(res_mm, paste0(out_prefix, "_MM.csv"), row.names = FALSE)
#
#         if (do_pairwise_for_continuous) {
#             cont_vars_mm <- vars_mm[sapply(df_mm[vars_mm], is.numeric)]
#             sig_cont_mm <- res_mm %>%
#                 filter(Test == "Kruskal-Wallis", Variable %in% cont_vars_mm, !is.na(P_value), P_value < pairwise_alpha) %>%
#                 pull(Variable)
#
#             pw_list_mm <- lapply(sig_cont_mm, function(v) pairwise_continuous(df_mm, "Q", v, pairwise_p_adjust))
#             pw_mm <- bind_rows(pw_list_mm)
#             if (nrow(pw_mm) > 0) {
#                 write.csv(pw_mm, paste0(out_prefix, "_MM_pairwiseWilcox.csv"), row.names = FALSE)
#             }
#         }
#     } else {
#         message("No MM samples found with ", disease_col, " == ", mm_value, "; skipping MM-only table.")
#     }
# } else {
#     message("disease_col not provided or not in df; skipping MM-only table.")
# }

# ----------------------------
# 5) VISUALIZATION (Boxplots & Barplots)
# ----------------------------

pdf(paste0(out_prefix, "_plots.pdf"), width = 8, height = 6)

my_comparisons <- list(
    c("Q1", "Q2"), c("Q1", "Q3"), c("Q1", "Q4"),
    c("Q2", "Q3"), c("Q2", "Q4"),
    c("Q3", "Q4")
)

for (var in vars_exist) {
    # check if var is continuous or categorical
    val <- df[[var]]

    if (is.numeric(val)) {
        # Boxplot with Kruskal-Wallis p-value
        p <- ggplot(df, aes(x = Q, y = .data[[var]], fill = Q)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.2, alpha = 0.3) +
            # stat_compare_means(method = "kruskal.test", label.y.npc = "top") +
            ggpubr::stat_compare_means(
                comparisons = my_comparisons,
                method = "wilcox.test", # or "t.test"
                label = "p.signif"
            ) +
            theme_pubr() +
            labs(title = paste0(var, " vs Q1-Q4"), x = "Quantile", y = var) +
            theme(legend.position = "none")
        print(p)
    } else {
        # Stacked Barplot with Chi-square p-value
        # Filter NA
        df_sub <- df %>% filter(!is.na(.data[[var]]))

        p <- ggplot(df_sub, aes(x = Q, fill = .data[[var]])) +
            geom_bar(position = "fill") +
            # scale_fill_viridis_d(option = "D")+
            ggsci::scale_fill_igv() +
            scale_y_continuous(labels = scales::percent) +
            theme_pubr() +
            labs(title = paste0(var, " distribution by Q"), x = "Quantile", y = "Proportion")
        print(p)
    }
}

dev.off()

# ----------------------------
# 6) OPTIONAL: Pretty Table-1 using gtsummary
# ----------------------------
# Check if gtsummary is installed
# if (requireNamespace("gtsummary", quietly = TRUE)) {
#     library(gtsummary)
#
#     # Define vars to summarize
#     vars_tbl <- vars_exist
#
#     # Create Table 1
#     tbl1 <- df %>%
#         select(Q, all_of(vars_tbl)) %>%
#         tbl_summary(
#             by = Q,
#             statistic = list(
#                 all_continuous() ~ "{median} ({p25}, {p75})",
#                 all_categorical() ~ "{n} ({p}%)"
#             ),
#             missing = "no" # or "ifany"
#         ) %>%
#         add_p(
#             test = list(
#                 all_continuous() ~ "kruskal.test",
#                 all_categorical() ~ "chisq.test"
#             )
#         ) %>%
#         add_overall() %>%
#         bold_p() %>%
#         modify_header(label = "**Variable**")
#
#     # Print to console
#     print(tbl1)
#
#     # Save as HTML/Docx if desired (requires flextable for docx)
#     # as_gt(tbl1) %>% gt::gtsave(paste0(out_prefix, "_Table1.html"))
# } else {
#     message("Package 'gtsummary' not installed. Skipping Table 1 generation.")
# }


# ----------------------------
# 7) SUMMARY: N patients & Mean L1PA2
# ----------------------------
# if ("L1PA2" %in% colnames(df)) {
if ("L1_num" %in% colnames(df)) {
    l1_summary <- df %>%
        group_by(Q) %>%
        dplyr::summarise(
            n_patients = dplyr::n(),
            # mean_L1PA2 = mean(L1PA2, na.rm = TRUE),
            mean_L1 = mean(L1_num, na.rm = TRUE),
            # sd_L1PA2 = sd(L1PA2, na.rm = TRUE)
            sd_L1 = sd(L1_num, na.rm = TRUE)
        )

    cat("\n=== Summary: N patients and Mean L1PA2 per Group ===\n")
    print(l1_summary)

    write.csv(l1_summary, paste0(out_prefix, "_L1_summary.csv"), row.names = FALSE)
} else {
    message("L1PA2 column not found, skipping summary.")
}


# low purity vs L1 number cor
p_cor <- df %>%
    dplyr::filter(!is.na(`Low purity`)) %>%
    ggplot(aes(x = `Low purity`, y = L1_num)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", color = "blue", fill = "lightgray") +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    labs(
        x = "Tumor Purity", # Rename axis label for clarity
        y = "L1 Number",
        title = "Correlation: L1 Number vs Tumor Purity"
    ) +
    theme_pubr()

print(p_cor)

# If categorical:
df$Low_purity_group <- ifelse(df$`Low purity` > 0.5, "yes", "no")
p_lowpurity_group <- df %>%
    dplyr::filter(!is.na(`Low purity`)) %>%
    ggplot(aes(x = Low_purity_group, y = L1_num)) + # , fill =  L1_quantile_group
    geom_boxplot() +
    stat_compare_means(
        method = "wilcox.test", # or "t.test"
        label = "p.signif"
    ) +
    theme_pubr() +
    labs(x = "", y = "L1", fill = "Low purity")

df %>%
    dplyr::filter(!is.na(`Low purity`)) %>%
    ggplot(aes(x = L1_quantile_group, fill = Low_purity_group)) +
    geom_bar(position = "fill") +
    # geom_text(stat='count', aes(label=..count..), vjust=-0.5) + # Adds the count numbers on top
    theme_pubr() +
    labs(x = "Low purity group", y = "Count")

# 3. Save it
ggsave(paste0(out_prefix, "_L1_vs_Purity_group.pdf"), plot = p_lowpurity_group, width = 6, height = 5)


# ----------------------------
# 8) Subtype Specific Enrichment Analysis (One-vs-Rest)
# ----------------------------
target_col <- "RNA_Subtype_Name" # Or make this dynamic based on your needs

if (target_col %in% colnames(df)) {
    message(paste0("Calculating specific enrichment p-values for ", target_col, "..."))

    subtypes <- unique(df[[target_col]])
    subtypes <- subtypes[!is.na(subtypes)]

    subtype_res_list <- list()

    for (st in subtypes) {
        # Create binary vector for this subtype: TRUE if matches current subtype
        is_target_subtype <- df[[target_col]] == st

        # Test against each Q group (e.g. Q4 vs Not-Q4)
        for (q_grp in levels(df$Q)) {
            # Binary vector for this Q group: TRUE if matches current Q
            is_in_q <- df$Q == q_grp

            # Create 2x2 Contingency Table
            #                 In Group Q   Not in Group Q
            #    Is Subtype      a              b
            #    Not Subtype     c              d

            # table(Row, Col) -> table(is_target_subtype, is_in_q)
            # TRUE/TRUE = Is Subtype AND In Q
            tbl <- table(
                factor(is_target_subtype, levels = c(TRUE, FALSE)),
                factor(is_in_q, levels = c(TRUE, FALSE))
            )

            # Fisher's Exact Test
            # alternative="greater" tests for ENRICHMENT (is subtype more likely to be in Q?)
            if (all(dim(tbl) == c(2, 2))) {
                ft <- fisher.test(tbl, alternative = "greater")

                subtype_res_list[[paste(st, q_grp, sep = "_")]] <- tibble(
                    Variable = target_col,
                    Subtype = st,
                    Enriched_In_Group = q_grp,
                    P_value = ft$p.value,
                    Odds_Ratio = as.numeric(ft$estimate),
                    Count_In_Group = tbl[1, 1], # TRUE, TRUE
                    Total_In_Group = sum(tbl[, 1])
                )
            }
        }
    }

    if (length(subtype_res_list) > 0) {
        subtype_enrichment_df <- bind_rows(subtype_res_list) %>%
            mutate(FDR = p.adjust(P_value, method = "BH")) %>%
            arrange(P_value)

        cat("\n=== Subtype Enrichment P-values (Top 10) ===\n")
        print(head(subtype_enrichment_df, 10))

        write.csv(subtype_enrichment_df, paste0(out_prefix, "_", target_col, "_Enrichment_Pvals.csv"), row.names = FALSE)
    }
}
