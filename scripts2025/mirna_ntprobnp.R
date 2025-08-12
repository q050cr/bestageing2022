# ntprobnp - miRNA analysis
#
# This script performs analysis on miRNA data related to NT-proBNP levels
# Date: 2025-08-05

# DE analysis for manuscript was performed primarily in:
# scripts/001c_model_de_analysis_DET_MATRIX
#   + shortened script to recalc CI's January 2024 --> scripts/org/001c_model_de_analysis_DET_MATRIX_recalc_ci_and_matched.R

# ----- Setup -----
data_path_bestageing2022 <- "/mnt/nas185/reich/rockerprojects/bestageing2022"
data_path_BestAgeing <- "/mnt/nas185/reich/BestAgeing"


# Load required libraries
library(tidyverse)
library(limma)
library(edgeR)
library(ggplot2)
library(readxl)
library(writexl)
library(pROC)
library(readxl)
library(janitor)
library(glue)
library(ggthemes)
library(pROC)
# Load required libraries for ML predictions
library(workflows)
library(tune)
library(tidymodels)
library(colino)

source(
  file = glue("{data_path_bestageing2022}/scripts/helper/custom_ggplot_theme.R")
)

# functions ---------------------------------------------------------------
convert_mir_name <- function(name) {
  # Replace 'mir' with 'miR'
  name <- gsub("mir", "miR", name)

  # Replace underscores with hyphens
  name <- gsub("_", "-", name)

  return(name)
}
convert_mir_name_V <- Vectorize(convert_mir_name)

mode_function <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

remove_digits <- function(vec) {
  # Use regex to replace patterns like '-DIGITS' with an empty string
  return(sub("-\\d+$", "", vec))
}

# Define a function for the permutation test
permute_test <- function(data, column_name, nperm = 1000) {
  # Split data by disease status
  group1 <- data[data$disease == 'control', column_name]
  group2 <- data[data$disease != 'control', column_name] # Assuming all other are disease

  # Compute observed t-statistic
  observed_t <- abs(t.test(group1, group2)$statistic)

  # INIT vector with length "nperm" to store permuted t-values
  permuted_t <- numeric(nperm)

  # Loop for permutations
  for (i in 1:nperm) {
    shuffled_disease <- sample(data$disease) # Shuffle disease status!!
    perm_group1 <- data[shuffled_disease == 'control', column_name]
    perm_group2 <- data[shuffled_disease != 'control', column_name]
    permuted_t[i] <- abs(t.test(perm_group1, perm_group2)$statistic)
  }

  # Compute p-value as proportion of times permuted t-values exceed observed t-value
  p_value <- mean(permuted_t >= observed_t)

  return(list(observed_t = observed_t, p_value = p_value))
}


# load data ---------------------------------------------------------------

path2dataprocessed_dcm <- glue(
  glue(
    "{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_dcm_data01.rds"
  )
)
path2dataprocessed_hfref <- glue(
  glue(
    "{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_hfref_data01.rds"
  )
)
data_processed_dcm <- readRDS(path2dataprocessed_dcm)
data_processed_hfref <- readRDS(path2dataprocessed_hfref)

# Load metadata with correct columns for each group
meta_dcm <- read_excel("/mnt/nas185/reich/BestAgeing/data/pheno_dcm.xlsx") %>%
  mutate(
    BestAgeingCode,
    phenodat = "dcm",
    ntprobnp = as.numeric(`DCM-Lab-ProBNP`),
    crea = as.numeric(`DCM-Lab-Cr`),
    lvdil_henry = as.numeric(`DCM-ECH-LVDil-Henry`),
    lvdil_lvedd = as.numeric(`DCM-ECH-LVDil-LVEDD`),
    lvsyst_ef = as.numeric(`DCM-ECH-LVSyst-EF`),
    symptoms = `DCM-PheDef-Symptoms`,
    wbc = as.numeric(`DCM-Lab-WBC`),
    .keep = "none"
  )

meta_hfref <- read_excel(
  "/mnt/nas185/reich/BestAgeing/data/pheno_hfref.xlsx"
) %>%
  mutate(
    BestAgeingCode,
    phenodat = "hfref",
    ntprobnp = as.numeric(`SD-Lab-ProBNP`),
    crea = as.numeric(`SD-Lab-Cr`),
    lvdil_henry = as.numeric(`SD-ECH-LVDil-Henry`),
    lvdil_lvedd = as.numeric(`SD-ECH-LVDil-LVEDD`),
    lvsyst_ef = as.numeric(`SD-ECH-LVSyst-EF`),
    dyspnea = `SD-CliExam-Dyspnea`,
    wbc = as.numeric(`SD-Lab-WBC`),
    .keep = "none"
  )

meta_control <- read_excel(
  "/mnt/nas185/reich/BestAgeing/data/pheno_controls.xlsx"
) %>%
  mutate(
    BestAgeingCode,
    phenodat = "control",
    ntprobnp = as.numeric(`Controls-Lab-ProBNP`),
    crea = as.numeric(`Controls-Lab-Cr`),
    .keep = "none"
  )

meta_all <- bind_rows(meta_dcm, meta_hfref, meta_control)

# Merge with processed data and add preprocessing column
dcm_merged <- data_processed_dcm %>%
  left_join(meta_all, by = c("pat_id" = "BestAgeingCode")) %>%
  mutate(preprocessing = "dcm_pipeline")

hfref_merged <- data_processed_hfref %>%
  left_join(meta_all, by = c("pat_id" = "BestAgeingCode")) %>%
  mutate(preprocessing = "hfref_pipeline")

# Data CLEANING for NT-proBNP analysis --------------------------------
clean_ntprobnp_data <- function(data, disease_name) {
  cat("Initial data dimensions:", dim(data), "\n")

  # 1) Remove controls with NT-proBNP > 450
  controls_high_bnp <- data$disease == "control" &
    !is.na(data$ntprobnp) &
    data$ntprobnp > 450
  n_removed_controls <- sum(controls_high_bnp, na.rm = TRUE)

  if (n_removed_controls > 0) {
    data <- data[!controls_high_bnp, ]
    cat("Removed", n_removed_controls, "controls with NT-proBNP > 450\n")
  } else {
    cat("No controls with NT-proBNP > 450 found\n")
  }

  # 2) For GLM analysis: use mean imputation for controls with missing NT-proBNP
  control_ntprobnp_mean <- mean(
    data$ntprobnp[data$disease == "control"],
    na.rm = TRUE
  )
  cat(
    "Control NT-proBNP mean for imputation:",
    round(control_ntprobnp_mean, 2),
    "\n"
  )

  # Impute missing NT-proBNP values in controls only
  missing_control_bnp <- data$disease == "control" & is.na(data$ntprobnp)
  n_imputed <- sum(missing_control_bnp, na.rm = TRUE)
  data$ntprobnp[missing_control_bnp] <- control_ntprobnp_mean

  if (n_imputed > 0) {
    cat(
      "Imputed NT-proBNP for",
      n_imputed,
      "controls using mean =",
      round(control_ntprobnp_mean, 2),
      "\n"
    )
  }

  # Store imputation value for later filtering in enhanced analyses
  # This allows us to remove imputed values when we want actual measurements only
  # We store it as an attribute to make the filtering dynamic and transparent
  attr(data, "control_imputation_value") <- control_ntprobnp_mean

  cat("Imputation value stored for dynamic filtering in enhanced analyses\n")

  # 3) Remove all patients (cases and controls) with missing NT-proBNP for GLM analysis
  data_glm <- data[!is.na(data$ntprobnp), ]
  n_removed_missing <- nrow(data) - nrow(data_glm)
  if (n_removed_missing > 0) {
    cat(
      "Removed",
      n_removed_missing,
      "patients with missing NT-proBNP for GLM analysis\n"
    )
  }

  cat("Final data dimensions for GLM:", dim(data_glm), "\n")
  cat(
    "NT-proBNP availability - Cases:",
    sum(!is.na(data_glm$ntprobnp[data_glm$disease == disease_name])),
    "Controls:",
    sum(!is.na(data_glm$ntprobnp[data_glm$disease == "control"])),
    "\n"
  )

  return(list(full_data = data, glm_data = data_glm))
}

# Clean data for both datasets
dcm_cleaned <- clean_ntprobnp_data(dcm_merged, "dcm")
hfref_cleaned <- clean_ntprobnp_data(hfref_merged, "hfref")

# Ensure imputation values are stored (robust handling)
# This handles cases where the function was updated after data processing
if (is.null(attr(dcm_cleaned$full_data, "control_imputation_value"))) {
  # Calculate the actual imputation value used
  dcm_actual_imputation <- mean(
    dcm_cleaned$full_data$ntprobnp[dcm_cleaned$full_data$disease == "control"],
    na.rm = TRUE
  )
  attr(
    dcm_cleaned$full_data,
    "control_imputation_value"
  ) <- dcm_actual_imputation
  cat(
    "Added actual imputation value for DCM:",
    round(dcm_actual_imputation, 2),
    "\n"
  )
}
if (is.null(attr(hfref_cleaned$full_data, "control_imputation_value"))) {
  # Calculate the actual imputation value used
  hfref_actual_imputation <- mean(
    hfref_cleaned$full_data$ntprobnp[
      hfref_cleaned$full_data$disease == "control"
    ],
    na.rm = TRUE
  )
  attr(
    hfref_cleaned$full_data,
    "control_imputation_value"
  ) <- hfref_actual_imputation
  cat(
    "Added actual imputation value for HFrEF:",
    round(hfref_actual_imputation, 2),
    "\n"
  )
}

# ------------------#
# DE Analysis Function (adapted from 001c_model_de_analysis_DET_MATRIX_recalc_ci_and_matched.R) ---------------------------
# ------------------#

run_de_analysis <- function(data_full, data_glm, disease_name) {
  # Use full data for univariate analysis, GLM data for multivariate
  data <- data_full

  # Remove columns not needed for miRNA analysis
  all_filtered_mirnas <- data %>%
    select(
      -any_of(c(
        "disease",
        "age",
        "sex",
        "phenodat",
        "ntprobnp",
        "crea",
        "preprocessing"
      ))
    ) %>%
    select(-contains("lvdil")) %>%
    select(-contains("lvsyst")) %>%
    select(-contains("dyspnea")) %>%
    select(-contains("symptoms")) %>%
    select(-contains("wbc"))

  cat("Filtered miRNA data dimensions:", dim(all_filtered_mirnas), "\n")
  cat("GLM data dimensions:", dim(data_glm), "\n")

  # Initialize result vectors
  n_mirnas <- ncol(all_filtered_mirnas) - 1 # Subtract pat_id column
  aucs <- aucs_lowerci <- aucs_upperci <- auc_glm <- auc_glm_ntprobnp <- auc_glm_ntprobnp_crea <- rep(
    NA,
    n_mirnas
  )
  pval_glm <- pval_glm_ntprobnp <- pval_glm_ntprobnp_crea <- log2fc <- rep(
    NA,
    n_mirnas
  )
  name_mir <- rep(NA, n_mirnas)

  # NT-proBNP only model for comparison
  data_glm_complete <- data_glm[
    !is.na(data_glm$ntprobnp) & !is.na(data_glm$age) & !is.na(data_glm$sex),
  ]

  f_ntprobnp_only <- as.formula("disease ~ ntprobnp + age + sex")
  glm_ntprobnp_only <- tryCatch(
    {
      glm(f_ntprobnp_only, data = data_glm_complete, family = binomial())
    },
    error = function(e) NULL
  )

  auc_ntprobnp_only <- NA
  if (!is.null(glm_ntprobnp_only)) {
    pred_ntprobnp_only <- predict(glm_ntprobnp_only, type = "response")
    auc_ntprobnp_only <- as.numeric(
      roc(data_glm_complete$disease, pred_ntprobnp_only, quiet = TRUE)$auc
    )
    cat(
      "NT-proBNP + age + sex baseline AUC:",
      round(auc_ntprobnp_only, 3),
      "\n"
    )
  }

  # Indexing for univariate analysis
  cont_index <- data$disease == "control"
  case_index <- data$disease == disease_name

  pb <- txtProgressBar(min = 0, max = n_mirnas, style = 3)

  for (i in 1:n_mirnas) {
    mirna_col <- i + 4 # Skip pat_id column
    cont <- as.numeric(data[cont_index, mirna_col, drop = TRUE])
    case <- as.numeric(data[case_index, mirna_col, drop = TRUE])

    # Remove NAs
    cont <- cont[!is.na(cont)]
    case <- case[!is.na(case)]

    # Skip if insufficient data
    if (length(cont) < 2 || length(case) < 2) {
      name_mir[i] <- names(data)[mirna_col]
      next
    }

    # Basic stats
    log2fc[i] <- mean(case) - mean(cont)
    name_mir[i] <- names(data)[mirna_col]

    # Univariate AUC with CI
    roc_obj <- roc(controls = cont, cases = case, quiet = TRUE)
    aucs[i] <- as.numeric(roc_obj$auc)
    auc_ci <- ci(roc_obj, quiet = TRUE)
    aucs_lowerci[i] <- auc_ci[1]
    aucs_upperci[i] <- auc_ci[3]

    # GLM models using cleaned data
    mirna_name <- names(data)[mirna_col]

    # Check if miRNA exists in GLM data
    if (!mirna_name %in% names(data_glm)) {
      next
    }

    # Model 1: miRNA + age + sex
    data_complete <- data_glm[
      !is.na(data_glm[[mirna_name]]) &
        !is.na(data_glm$age) &
        !is.na(data_glm$sex),
    ]

    if (nrow(data_complete) > 10) {
      f1 <- as.formula(paste("disease ~", mirna_name, "+ age + sex"))
      glm1 <- tryCatch(
        {
          glm(f1, data = data_complete, family = binomial())
        },
        error = function(e) NULL
      )

      if (!is.null(glm1)) {
        pred1 <- predict(glm1, type = "response")
        auc_glm[i] <- as.numeric(
          roc(data_complete$disease, pred1, quiet = TRUE)$auc
        )
        pval_glm[i] <- coef(summary(glm1))[2, 4]
      }

      # Model 2: miRNA + age + sex + ntprobnp (incremental benefit analysis)
      data_complete_nt <- data_complete[!is.na(data_complete$ntprobnp), ]

      if (nrow(data_complete_nt) > 10) {
        f2 <- as.formula(paste(
          "disease ~",
          mirna_name,
          "+ age + sex + ntprobnp"
        ))
        glm2 <- tryCatch(
          {
            glm(f2, data = data_complete_nt, family = binomial())
          },
          error = function(e) NULL
        )

        if (!is.null(glm2)) {
          pred2 <- predict(glm2, type = "response")
          auc_glm_ntprobnp[i] <- as.numeric(
            roc(data_complete_nt$disease, pred2, quiet = TRUE)$auc
          )
          pval_glm_ntprobnp[i] <- coef(summary(glm2))[2, 4]
        }

        # Model 3: miRNA + age + sex + ntprobnp + crea (sensitivity analysis)
        if (!all(is.na(data_complete_nt$crea))) {
          data_complete_crea <- data_complete_nt[
            !is.na(data_complete_nt$crea),
          ]

          if (nrow(data_complete_crea) > 10) {
            f3 <- as.formula(paste(
              "disease ~",
              mirna_name,
              "+ age + sex + ntprobnp + crea"
            ))
            glm3 <- tryCatch(
              {
                glm(f3, data = data_complete_crea, family = binomial())
              },
              error = function(e) NULL
            )

            if (!is.null(glm3)) {
              pred3 <- predict(glm3, type = "response")
              auc_glm_ntprobnp_crea[i] <- as.numeric(
                roc(data_complete_crea$disease, pred3, quiet = TRUE)$auc
              )
              pval_glm_ntprobnp_crea[i] <- coef(summary(glm3))[2, 4]
            }
          }
        }
      }
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  # Compile results
  results <- tibble(
    miRNA = name_mir,
    log2FoldChange = log2fc,
    auc_univariate = aucs,
    auc_lowerci = aucs_lowerci,
    auc_upperci = aucs_upperci,
    auc_glm = auc_glm,
    auc_glm_ntprobnp = auc_glm_ntprobnp,
    auc_glm_ntprobnp_crea = auc_glm_ntprobnp_crea,
    pval_glm = pval_glm,
    pval_glm_ntprobnp = pval_glm_ntprobnp,
    pval_glm_ntprobnp_crea = pval_glm_ntprobnp_crea,
    auc_ntprobnp_baseline = auc_ntprobnp_only # Add baseline for comparison
  ) %>%
    mutate(
      auc_univariate = ifelse(
        auc_univariate < 0.5,
        1 - auc_univariate,
        auc_univariate
      ),
      padj_glm = p.adjust(pval_glm, method = "holm"),
      padj_glm_ntprobnp = p.adjust(pval_glm_ntprobnp, method = "holm"),
      padj_glm_ntprobnp_crea = p.adjust(
        pval_glm_ntprobnp_crea,
        method = "holm"
      ),
      incremental_benefit_vs_baseline = auc_glm - auc_ntprobnp_baseline, # miRNA vs NT-proBNP baseline
      incremental_benefit_combined = auc_glm_ntprobnp - auc_ntprobnp_baseline, # combined vs NT-proBNP baseline
      incremental_benefit_miRNA_added = auc_glm_ntprobnp - auc_glm, # adding miRNA to baseline
      incremental_benefit_crea = auc_glm_ntprobnp_crea - auc_glm_ntprobnp
    ) %>%
    arrange(desc(auc_univariate))

  return(results)
}

# Run DE analysis for both datasets
cat("Running DE analysis for DCM...\n")
dcm_results <- run_de_analysis(
  dcm_cleaned$full_data,
  dcm_cleaned$glm_data,
  "dcm"
)

cat("Running DE analysis for HFrEF...\n")
hfref_results <- run_de_analysis(
  hfref_cleaned$full_data,
  hfref_cleaned$glm_data,
  "hfref"
)

# NT-proBNP baseline performance (using original univariate approach)
ntprobnp_performance <- function(data, disease_name) {
  if (all(is.na(data$ntprobnp))) {
    return(NA)
  }

  data_clean <- data %>% filter(!is.na(ntprobnp))
  cont <- data_clean[data_clean$disease == "control", "ntprobnp", drop = TRUE]
  case <- data_clean[
    data_clean$disease == disease_name,
    "ntprobnp",
    drop = TRUE
  ]

  roc_obj <- roc(controls = cont, cases = case, quiet = TRUE)
  return(as.numeric(roc_obj$auc))
}

ntprobnp_auc_dcm <- ntprobnp_performance(dcm_cleaned$full_data, "dcm")
ntprobnp_auc_hfref <- ntprobnp_performance(hfref_cleaned$full_data, "hfref")

cat(paste("NT-proBNP univariate AUC DCM:", round(ntprobnp_auc_dcm, 3), "\n"))
cat(paste(
  "NT-proBNP univariate AUC HFrEF:",
  round(ntprobnp_auc_hfref, 3),
  "\n"
))

# Display top results with incremental benefit
cat("\nTop 10 miRNAs DCM (with incremental benefit analysis):\n")
dcm_top <- dcm_results %>%
  filter(!is.na(auc_univariate) & miRNA != "disease") %>%
  select(
    miRNA,
    auc_univariate,
    auc_glm,
    auc_glm_ntprobnp,
    auc_ntprobnp_baseline,
    incremental_benefit_vs_baseline,
    incremental_benefit_combined,
    padj_glm
  ) %>%
  head(10)
print(dcm_top)

cat("\nTop 10 miRNAs HFrEF (with incremental benefit analysis):\n")
hfref_top <- hfref_results %>%
  filter(!is.na(auc_univariate) & miRNA != "disease") %>%
  select(
    miRNA,
    auc_univariate,
    auc_glm,
    auc_glm_ntprobnp,
    auc_ntprobnp_baseline,
    incremental_benefit_vs_baseline,
    incremental_benefit_combined,
    padj_glm
  ) %>%
  head(10)
print(hfref_top)

# Save results
if (exists("dcm_results")) {
  write_csv(dcm_results, "output/dcm_mirna_ntprobnp_results.csv")
  cat("DCM results saved to output/dcm_mirna_ntprobnp_results.csv\n")
}

if (exists("hfref_results")) {
  write_csv(hfref_results, "output/hfref_mirna_ntprobnp_results.csv")
  cat("HFrEF results saved to output/hfref_mirna_ntprobnp_results.csv\n")
}

# Summary statistics
cat("\n=== SUMMARY ===\n")
cat("NT-proBNP baseline AUC (univariate):\n")
cat("  DCM:", round(ntprobnp_auc_dcm, 3), "\n")
cat("  HFrEF:", round(ntprobnp_auc_hfref, 3), "\n")

# Get NT-proBNP multivariate baseline from results
ntprobnp_baseline_dcm <- unique(dcm_results$auc_ntprobnp_baseline)[1]
ntprobnp_baseline_hfref <- unique(hfref_results$auc_ntprobnp_baseline)[1]

if (!is.na(ntprobnp_baseline_dcm)) {
  cat("NT-proBNP + age + sex baseline AUC:\n")
  cat("  DCM:", round(ntprobnp_baseline_dcm, 3), "\n")
  cat("  HFrEF:", round(ntprobnp_baseline_hfref, 3), "\n")
}

cat("\nBest miRNA univariate AUC:\n")
if (exists("dcm_results")) {
  best_dcm_auc <- max(
    dcm_results$auc_univariate[dcm_results$miRNA != "disease"],
    na.rm = TRUE
  )
  cat("  DCM:", round(best_dcm_auc, 3), "\n")
}
if (exists("hfref_results")) {
  best_hfref_auc <- max(
    hfref_results$auc_univariate[hfref_results$miRNA != "disease"],
    na.rm = TRUE
  )
  cat("  HFrEF:", round(best_hfref_auc, 3), "\n")
}

# Additional Analysis: miRNA-NT-proBNP Correlation Analysis ----------------------------------------
correlation_analysis <- function(data, disease_name) {
  # Use only patients with available NT-proBNP (no imputation)
  data_corr <- data %>% filter(!is.na(ntprobnp), ntprobnp > 0)

  # Log-transform NT-proBNP to handle skewness
  data_corr$log_ntprobnp <- log(data_corr$ntprobnp)

  cat("\nCorrelation analysis for", disease_name, "\n")
  cat("Samples with available NT-proBNP:", nrow(data_corr), "\n")
  cat("NT-proBNP range:", round(range(data_corr$ntprobnp), 1), "\n")

  # Get miRNA columns only
  mirna_cols <- grep("hsa_", names(data_corr), value = TRUE)

  # Calculate correlations
  correlations <- map_dfr(mirna_cols, function(mir) {
    if (sum(!is.na(data_corr[[mir]])) < 10) {
      return(NULL)
    }

    cor_test <- cor.test(
      data_corr[[mir]],
      data_corr$log_ntprobnp,
      method = "spearman"
    )
    tibble(
      miRNA = mir,
      correlation = cor_test$estimate,
      p_value = cor_test$p.value,
      n_samples = sum(!is.na(data_corr[[mir]]))
    )
  }) %>%
    filter(!is.na(correlation)) %>%
    mutate(
      abs_correlation = abs(correlation),
      padj = p.adjust(p_value, method = "holm")
    ) %>%
    arrange(desc(abs_correlation))

  return(list(correlations = correlations, data = data_corr))
}

# Run correlation analysis
dcm_corr <- correlation_analysis(dcm_cleaned$full_data, "DCM")
hfref_corr <- correlation_analysis(hfref_cleaned$full_data, "HFrEF")

# Head-to-head performance comparison
performance_comparison <- function(results, ntprobnp_auc, disease_name) {
  top_mirnas <- results %>%
    filter(miRNA != "disease", !is.na(auc_univariate)) %>%
    slice_max(auc_univariate, n = 5)

  cat("\n=== HEAD-TO-HEAD COMPARISON:", disease_name, "===\n")
  cat("NT-proBNP AUC:", round(ntprobnp_auc, 3), "\n")
  cat("Top 5 miRNAs vs NT-proBNP:\n")

  for (i in 1:min(5, nrow(top_mirnas))) {
    mir_auc <- top_mirnas$auc_univariate[i]
    mir_name <- top_mirnas$miRNA[i]
    difference <- mir_auc - ntprobnp_auc
    cat(sprintf("  %s: AUC=%.3f (Δ=%.3f)\n", mir_name, mir_auc, difference))
  }

  return(top_mirnas)
}

dcm_comparison <- performance_comparison(dcm_results, ntprobnp_auc_dcm, "DCM")
hfref_comparison <- performance_comparison(hfref_results, ntprobnp_auc_hfref, "HFrEF")

# Display correlation results -----------------------------------------------------------------------------
cat("\n=== CORRELATION WITH LOG(NT-proBNP) ===\n")
cat("Top 10 correlated miRNAs - DCM:\n")
print(dcm_corr$correlations %>% head(10) %>% select(miRNA, correlation, padj))

cat("\nTop 10 correlated miRNAs - HFrEF:\n")
print(hfref_corr$correlations %>% head(10) %>% select(miRNA, correlation, padj))

# Save correlation results
write_csv(dcm_corr$correlations, "output/dcm_ntprobnp_correlations.csv")
write_csv(hfref_corr$correlations, "output/hfref_ntprobnp_correlations.csv")

# Enhanced Head-to-Head Performance Comparison with Confidence Intervals -----------------------------------
#
# IMPORTANT: This analysis uses only ACTUAL NT-proBNP measurements (not imputed values)
# to provide the most accurate comparison between miRNAs and NT-proBNP.
#
# Rationale: In clinical practice, physicians would have actual NT-proBNP values,
# so we want to compare miRNAs against real biomarker measurements, not imputed ones.
#
enhanced_performance_comparison <- function(
  data,
  results,
  ntprobnp_auc,
  disease_name
) {
  # Use data with actual NT-proBNP measurements only
  data_original <- data %>%
    filter(!is.na(ntprobnp))

  # Remove imputed control values using the stored imputation value
  imputation_value <- attr(data, "control_imputation_value")
  if (!is.null(imputation_value)) {
    data_original <- data_original %>%
      filter(!(disease == "control" & abs(ntprobnp - imputation_value) < 0.01))
    cat(
      "Removed imputed control values (",
      round(imputation_value, 2),
      " pg/mL)\n"
    )
  }

  cont <- data_original[
    data_original$disease == "control",
    "ntprobnp",
    drop = TRUE
  ]
  case <- data_original[
    data_original$disease == disease_name,
    "ntprobnp",
    drop = TRUE
  ]

  # Check if we have enough data
  if (length(cont) < 5 || length(case) < 5) {
    cat(
      "Insufficient non-imputed data for enhanced comparison in",
      disease_name,
      "\n"
    )
    cat("Controls:", length(cont), "Cases:", length(case), "\n")
    return(tibble())
  }

  ntprobnp_roc <- roc(controls = cont, cases = case, quiet = TRUE)
  ntprobnp_ci <- ci(ntprobnp_roc, quiet = TRUE)

  top_mirnas <- results %>%
    filter(miRNA != "disease", !is.na(auc_univariate)) %>%
    slice_max(auc_univariate, n = 5)

  cat("\n=== ENHANCED HEAD-TO-HEAD COMPARISON:", disease_name, "===\n")
  cat(sprintf(
    "NT-proBNP: AUC=%.3f (95%% CI: %.3f-%.3f)\n",
    ntprobnp_auc,
    ntprobnp_ci[1],
    ntprobnp_ci[3]
  ))
  cat("\nTop 5 miRNAs with confidence intervals:\n")

  comparison_results <- tibble()

  for (i in 1:min(5, nrow(top_mirnas))) {
    mir_auc <- top_mirnas$auc_univariate[i]
    mir_ci_lower <- top_mirnas$auc_lowerci[i]
    mir_ci_upper <- top_mirnas$auc_upperci[i]
    mir_name <- top_mirnas$miRNA[i]

    # Check CI overlap
    ci_overlap <- (mir_ci_lower <= ntprobnp_ci[3]) &
      (mir_ci_upper >= ntprobnp_ci[1])
    overlap_status <- ifelse(
      ci_overlap,
      "Overlapping CIs",
      "Non-overlapping CIs"
    )

    cat(sprintf(
      "  %s: AUC=%.3f (95%% CI: %.3f-%.3f) - %s\n",
      mir_name,
      mir_auc,
      mir_ci_lower,
      mir_ci_upper,
      overlap_status
    ))

    comparison_results <- bind_rows(
      comparison_results,
      tibble(
        miRNA = mir_name,
        mirna_auc = mir_auc,
        mirna_ci_lower = mir_ci_lower,
        mirna_ci_upper = mir_ci_upper,
        ntprobnp_auc = ntprobnp_auc,
        ntprobnp_ci_lower = ntprobnp_ci[1],
        ntprobnp_ci_upper = ntprobnp_ci[3],
        ci_overlap = ci_overlap,
        difference = mir_auc - ntprobnp_auc
      )
    )
  }

  return(comparison_results)
}

# Run enhanced analyses
cat("\n", rep("=", 60), "\n")
cat("ENHANCED PERFORMANCE ANALYSES\n")
cat(rep("=", 60), "\n")

# Enhanced head-to-head comparison
dcm_enhanced <- enhanced_performance_comparison(
  data = dcm_cleaned$full_data,
  results = dcm_results,
  ntprobnp_auc = ntprobnp_auc_dcm,
  disease_name = "dcm"
)
hfref_enhanced <- enhanced_performance_comparison(
  hfref_cleaned$full_data,
  hfref_results,
  ntprobnp_auc_hfref,
  "hfref"
)


# NT-proBNP Gray Zone Analysis -----------------------------------------------------------------------
gray_zone_analysis <- function(data, results, disease_name) {
  # Use data with actual NT-proBNP measurements only
  data_original <- data %>%
    filter(!is.na(ntprobnp))

  # Remove imputed control values using the stored imputation value
  imputation_value <- attr(data, "control_imputation_value")
  if (!is.null(imputation_value)) {
    data_original <- data_original %>%
      filter(!(disease == "control" & abs(ntprobnp - imputation_value) < 0.01))
    cat(
      "Removed imputed control values (",
      round(imputation_value, 2),
      " pg/mL)\n"
    )
  }

  if (nrow(data_original) < 20) {
    cat(
      "Insufficient non-imputed data for gray zone analysis in",
      disease_name,
      "\n"
    )
    return(tibble())
  }

  # Calculate quartiles based on non-imputed data
  data_dcm_only <- data_original %>% filter(disease == disease_name)
  q25 <- quantile(data_dcm_only$ntprobnp, 0.25, na.rm = TRUE)
  median <- median(data_dcm_only$ntprobnp, na.rm = TRUE)
  q75 <- quantile(data_dcm_only$ntprobnp, 0.75, na.rm = TRUE)

  # Gray zone (lower than MEDIAN) ----------------------------------------
  gray_zone_data <- data_original %>%
    filter(ntprobnp < median)

  cat("\n=== GRAY ZONE ANALYSIS:", disease_name, "===\n")
  cat(sprintf(
    "NT-proBNP %s: %.1f [%.1f - %.1f] pg/mL (non-imputed data)\n",
    disease_name,
    median,
    q25,
    q75
  ))
  cat(sprintf(
    "Samples in analysis: %d (%.1f%% of non-imputed data)\n",
    nrow(gray_zone_data),
    100 * nrow(gray_zone_data) / nrow(data_original)
  ))

  # Test top DE miRNAs in gray zone ------------------------------
  top_mirnas <- results %>%
    filter(miRNA != "disease", !is.na(auc_univariate)) %>%
    slice_max(auc_univariate, n = 10)

  gray_zone_results <- tibble()

  for (i in 1:nrow(top_mirnas)) {
    mir_name <- top_mirnas$miRNA[i]

    if (!mir_name %in% names(gray_zone_data)) {
      next
    }

    # Prepare data for GLM analysis in gray zone
    gray_glm_data <- gray_zone_data %>%
      filter(!is.na(.data[[mir_name]]), !is.na(age), !is.na(sex)) %>%
      mutate(disease_binary = ifelse(disease == disease_name, 1, 0))

    if (nrow(gray_glm_data) < 20) {
      next
    }

    # GLM Models in Gray Zone with incremental value analysis
    # Model 1: miRNA + age + sex
    f1 <- as.formula(paste("disease_binary ~", mir_name, "+ age + sex"))
    glm1 <- tryCatch(
      glm(f1, data = gray_glm_data, family = binomial()),
      error = function(e) NULL
    )

    # Model 2: NT-proBNP + age + sex
    f2 <- as.formula("disease_binary ~ ntprobnp + age + sex")
    glm2 <- tryCatch(
      glm(f2, data = gray_glm_data, family = binomial()),
      error = function(e) NULL
    )

    # Model 3: Combined miRNA + NT-proBNP + age + sex
    f3 <- as.formula(paste(
      "disease_binary ~",
      mir_name,
      "+ ntprobnp + age + sex"
    ))
    glm3 <- tryCatch(
      glm(f3, data = gray_glm_data, family = binomial()),
      error = function(e) NULL
    )

    if (!is.null(glm1) && !is.null(glm2) && !is.null(glm3)) {
      # Calculate AUCs
      auc1 <- as.numeric(
        roc(
          gray_glm_data$disease_binary,
          predict(glm1, type = "response"),
          quiet = TRUE
        )$auc
      )
      auc2 <- as.numeric(
        roc(
          gray_glm_data$disease_binary,
          predict(glm2, type = "response"),
          quiet = TRUE
        )$auc
      )
      auc3 <- as.numeric(
        roc(
          gray_glm_data$disease_binary,
          predict(glm3, type = "response"),
          quiet = TRUE
        )$auc
      )

      # Statistical summary
      cat(sprintf("  %s Gray Zone Incremental Analysis:\n", mir_name))
      cat(sprintf("    miRNA + age + sex:           AUC=%.3f\n", auc1))
      cat(sprintf("    NT-proBNP + age + sex:       AUC=%.3f\n", auc2))
      cat(sprintf("    Combined + age + sex:        AUC=%.3f\n", auc3))
      cat(sprintf(
        "    Incremental benefit:         Δ=%.3f (p=%.3f)\n",
        auc3 - auc2,
        anova(glm2, glm3, test = "Chisq")$`Pr(>Chi)`[2]
      ))

      gray_zone_results <- bind_rows(
        gray_zone_results,
        tibble(
          miRNA = mir_name,
          auc_mirna_model = auc1,
          auc_ntprobnp_model = auc2,
          auc_combined_model = auc3,
          incremental_benefit = auc3 - auc2,
          statistical_test_p = anova(glm2, glm3, test = "Chisq")$`Pr(>Chi)`[2],
          n_gray_zone = nrow(gray_glm_data),
          n_cases_gray = sum(gray_glm_data$disease_binary == 1),
          n_controls_gray = sum(gray_glm_data$disease_binary == 0)
        )
      )
    }
  }

  return(gray_zone_results)
}


# Gray zone analysis
dcm_gray <- gray_zone_analysis(dcm_cleaned$full_data, dcm_results, "dcm")
hfref_gray <- gray_zone_analysis(hfref_cleaned$full_data, hfref_results, "hfref")

# Save enhanced analysis results
if (nrow(dcm_enhanced) > 0) {
  write_csv(dcm_enhanced, "output/dcm_enhanced_comparison.csv")
}
if (nrow(hfref_enhanced) > 0) {
  write_csv(hfref_enhanced, "output/hfref_enhanced_comparison.csv")
}
if (nrow(dcm_gray) > 0) {
  write_csv(dcm_gray, "output/dcm_gray_zone_analysis.csv")
}
if (nrow(hfref_gray) > 0) {
  write_csv(hfref_gray, "output/hfref_gray_zone_analysis.csv")
}


# Clinical Reclassification Analysis using CLINICAL NT-proBNP thresholds 125-450ng/ml ----------------------
reclassification_analysis <- function(data, results, disease_name) {
  # Use data with actual NT-proBNP measurements only
  data_original <- data %>%
    filter(!is.na(ntprobnp))

  # Remove imputed control values using the stored imputation value
  imputation_value <- attr(data, "control_imputation_value")
  if (!is.null(imputation_value)) {
    data_original <- data_original %>%
      filter(!(disease == "control" & abs(ntprobnp - imputation_value) < 0.01))
    cat(
      "Removed imputed control values (",
      round(imputation_value, 2),
      " pg/mL)\n"
    )
  }

  if (nrow(data_original) < 20) {
    cat(
      "Insufficient non-imputed data for reclassification analysis in",
      disease_name,
      "\n"
    )
    return(list())
  }

  # Clinical NT-proBNP thresholds
  rule_out_threshold <- 125 # Rule out HF if < 125 pg/mL
  rule_in_threshold <- 450 # Rule in HF if > 450 pg/mL
  # Gray zone: 125-450 pg/mL

  cat("\n=== CLINICAL RECLASSIFICATION ANALYSIS:", disease_name, "===\n")
  cat("NT-proBNP Clinical Thresholds:\n")
  cat("  Rule-out (<125 pg/mL): Low probability of HF\n")
  cat("  Gray zone (125-450 pg/mL): Intermediate probability\n")
  cat("  Rule-in (>450 pg/mL): High probability of HF\n")

  # Classify patients by NT-proBNP thresholds
  data_classified <- data_original %>%
    mutate(
      ntprobnp_class = case_when(
        ntprobnp < rule_out_threshold ~ "Rule-out (<125)",
        ntprobnp >= rule_out_threshold & ntprobnp <= rule_in_threshold ~ "Gray zone (125-450)",
        ntprobnp > rule_in_threshold ~ "Rule-in (>450)",
        TRUE ~ NA_character_
      ),
      is_case = ifelse(disease == disease_name, 1, 0)
    )

  # Summary by zones
  zone_summary <- data_classified %>%
    group_by(ntprobnp_class) %>%
    summarise(
      n_total = n(),
      n_cases = sum(is_case),
      n_controls = sum(1 - is_case),
      case_percentage = round(100 * n_cases / n_total, 1),
      .groups = "drop"
    )

  cat("\nPatient distribution by NT-proBNP zones:\n")
  print(zone_summary)

  # Get top miRNAs for reclassification analysis
  top_mirnas <- results %>%
    filter(miRNA != "disease", !is.na(auc_univariate)) %>%
    slice_max(auc_univariate, n = 3)

  reclassification_results <- list()

  for (i in 1:nrow(top_mirnas)) {
    mir_name <- top_mirnas$miRNA[i]

    if (!mir_name %in% names(data_classified)) {
      next
    }

    # Focus on gray zone patients (most clinically relevant)
    gray_zone_patients <- data_classified %>%
      filter(
        ntprobnp_class == "Gray zone (125-450)",
        !is.na(.data[[mir_name]]),
        !is.na(age),
        !is.na(sex)
      )

    if (nrow(gray_zone_patients) < 10) {
      next
    }

    # Build GLM model for reclassification
    f_model <- as.formula(paste("is_case ~", mir_name, "+ age + sex"))
    glm_model <- tryCatch(
      glm(f_model, data = gray_zone_patients, family = binomial()),
      error = function(e) NULL
    )

    if (!is.null(glm_model)) {
      # Get predicted probabilities
      gray_zone_patients$pred_prob <- predict(glm_model, type = "response")

      # Define reclassification thresholds based on predicted probability
      # High confidence: >0.7 probability
      # Low confidence: <0.3 probability
      # Intermediate: 0.3-0.7 probability

      gray_zone_reclassified <- gray_zone_patients %>%
        mutate(
          mirna_reclass = case_when(
            pred_prob > 0.7 ~ "High risk (miRNA)",
            pred_prob < 0.3 ~ "Low risk (miRNA)",
            TRUE ~ "Intermediate (miRNA)"
          ),
          # Create reclassification categories
          reclass_category = case_when(
            mirna_reclass == "High risk (miRNA)" & is_case == 1 ~ "Correct rule-in",
            mirna_reclass == "High risk (miRNA)" & is_case == 0 ~ "Incorrect rule-in",
            mirna_reclass == "Low risk (miRNA)" & is_case == 0 ~ "Correct rule-out",
            mirna_reclass == "Low risk (miRNA)" & is_case == 1 ~ "Incorrect rule-out",
            TRUE ~ "No reclassification"
          )
        )

      # Calculate reclassification metrics
      reclass_summary <- gray_zone_reclassified %>%
        count(reclass_category) %>%
        mutate(percentage = round(100 * n / nrow(gray_zone_reclassified), 1))

      # Net Reclassification Improvement (NRI)
      correct_reclass <- sum(reclass_summary$n[
        reclass_summary$reclass_category %in%
          c("Correct rule-in", "Correct rule-out")
      ])
      incorrect_reclass <- sum(reclass_summary$n[
        reclass_summary$reclass_category %in%
          c("Incorrect rule-in", "Incorrect rule-out")
      ])

      nri <- (correct_reclass - incorrect_reclass) /
        nrow(gray_zone_reclassified)

      cat(sprintf(
        "\n%s Reclassification in Gray Zone (n=%d):\n",
        mir_name,
        nrow(gray_zone_reclassified)
      ))
      print(reclass_summary)
      cat(sprintf("Net Reclassification Improvement: %.3f\n", nri))

      reclassification_results[[mir_name]] <- list(
        data = gray_zone_reclassified,
        summary = reclass_summary,
        nri = nri,
        zone_summary = zone_summary
      )
    }
  }

  return(reclassification_results)
}

# Reclassification Visualization Function
create_reclassification_plots <- function(reclass_results, disease_name) {
  if (length(reclass_results) == 0) {
    cat("No reclassification results to plot for", disease_name, "\n")
    return(NULL)
  }

  plots <- list()

  for (mir_name in names(reclass_results)) {
    result <- reclass_results[[mir_name]]

    # Plot 1: Reclassification summary bar plot
    p1 <- result$summary %>%
      mutate(
        reclass_category = factor(
          reclass_category,
          levels = c(
            "Correct rule-in",
            "Correct rule-out",
            "Incorrect rule-in",
            "Incorrect rule-out",
            "No reclassification"
          )
        ),
        fill_color = case_when(
          reclass_category %in% c("Correct rule-in", "Correct rule-out") ~ "Improved",
          reclass_category %in% c("Incorrect rule-in", "Incorrect rule-out") ~ "Worsened",
          TRUE ~ "Unchanged"
        )
      ) %>%
      ggplot(aes(x = reclass_category, y = n, fill = fill_color)) +
      geom_col() +
      geom_text(
        aes(label = paste0(n, " (", percentage, "%)")),
        vjust = -0.3,
        size = 3
      ) +
      scale_fill_manual(
        values = c(
          "Improved" = "#2E8B57",
          "Worsened" = "#DC143C",
          "Unchanged" = "#696969"
        ),
        name = "Reclassification"
      ) +
      labs(
        title = paste("Clinical Reclassification by", mir_name),
        subtitle = paste(
          disease_name,
          "- Gray Zone Patients (NT-proBNP 125-450 pg/mL)"
        ),
        x = "Reclassification Category",
        y = "Number of Patients",
        caption = paste("NRI =", round(result$nri, 3))
      ) +
      theme_minimal(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )

    # Plot 2: Probability distribution by actual disease status
    p2 <- result$data %>%
      mutate(
        disease_status = ifelse(
          is_case == 1,
          paste(disease_name, "cases"),
          "Controls"
        )
      ) %>%
      ggplot(aes(x = pred_prob, fill = disease_status)) +
      geom_histogram(alpha = 0.7, bins = 20, position = "identity") +
      geom_vline(xintercept = c(0.3, 0.7), linetype = "dashed", color = "red") +
      scale_fill_manual(
        values = c("#1f77b4", "#ff7f0e"),
        name = "Patient Type"
      ) +
      labs(
        title = paste("Predicted Probability Distribution -", mir_name),
        subtitle = paste(disease_name, "- Gray Zone Patients"),
        x = "Predicted Probability of Disease",
        y = "Count",
        caption = "Dashed lines: Reclassification thresholds (0.3, 0.7)"
      ) +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )

    # Plot 3: NT-proBNP vs miRNA scatter with reclassification
    p3 <- result$data %>%
      mutate(
        disease_status = ifelse(
          is_case == 1,
          paste(disease_name, "cases"),
          "Controls"
        ),
        reclass_simple = case_when(
          reclass_category %in% c("Correct rule-in", "Correct rule-out") ~ "Improved classification",
          reclass_category %in% c("Incorrect rule-in", "Incorrect rule-out") ~ "Worsened classification",
          TRUE ~ "No change"
        )
      ) %>%
      ggplot(aes(
        x = ntprobnp,
        y = .data[[mir_name]],
        color = disease_status,
        shape = reclass_simple
      )) +
      geom_point(size = 2.5, alpha = 0.8) +
      geom_hline(
        yintercept = median(result$data[[mir_name]], na.rm = TRUE),
        linetype = "dashed",
        alpha = 0.5
      ) +
      geom_vline(
        xintercept = c(125, 450),
        linetype = "solid",
        alpha = 0.3,
        size = 1
      ) +
      scale_color_manual(
        values = c("#1f77b4", "#ff7f0e"),
        name = "Patient Type"
      ) +
      scale_shape_manual(
        values = c(16, 17, 15),
        name = "Reclassification"
      ) +
      labs(
        title = paste("NT-proBNP vs", mir_name, "- Reclassification Results"),
        subtitle = paste(disease_name, "- Gray Zone Analysis"),
        x = "NT-proBNP (pg/mL)",
        y = paste(mir_name, "Expression"),
        caption = "Vertical lines: NT-proBNP thresholds (125, 450 pg/mL)"
      ) +
      theme_economist_white() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.position = "bottom"
      )

    plots[[mir_name]] <- list(summary = p1, distribution = p2, scatter = p3)
  }

  return(plots)
}


# Run Reclassification Analysis
cat("\n", rep("=", 60), "\n")
cat("CLINICAL RECLASSIFICATION ANALYSIS\n")
cat(rep("=", 60), "\n")

# Run reclassification analysis
dcm_reclass <- reclassification_analysis(dcm_cleaned$full_data, dcm_results, "dcm")
hfref_reclass <- reclassification_analysis(hfref_cleaned$full_data, hfref_results, "hfref")

# Create and save plots
if (length(dcm_reclass) > 0) {
  dcm_plots <- create_reclassification_plots(dcm_reclass, "DCM")

  # Save plots for each miRNA
  for (mir_name in names(dcm_plots)) {
    # Clean miRNA name for filename
    clean_name <- gsub("[^A-Za-z0-9]", "_", mir_name)

    ggsave(
      filename = paste0(
        "figures/dcm_",
        clean_name,
        "_reclassification_summary.png"
      ),
      plot = dcm_plots[[mir_name]]$summary,
      width = 10,
      height = 6,
      dpi = 300
    )

    ggsave(
      filename = paste0(
        "figures/dcm_",
        clean_name,
        "_probability_distribution.png"
      ),
      plot = dcm_plots[[mir_name]]$distribution,
      width = 10,
      height = 6,
      dpi = 300
    )

    ggsave(
      filename = paste0(
        "figures/dcm_",
        clean_name,
        "_scatter_reclassification.png"
      ),
      plot = dcm_plots[[mir_name]]$scatter,
      width = 12,
      height = 8,
      dpi = 300
    )
  }

  # Save reclassification data
  for (mir_name in names(dcm_reclass)) {
    clean_name <- gsub("[^A-Za-z0-9]", "_", mir_name)
    write_csv(
      dcm_reclass[[mir_name]]$data,
      paste0("output/dcm_", clean_name, "_reclassification_data.csv")
    )
    write_csv(
      dcm_reclass[[mir_name]]$summary,
      paste0("output/dcm_", clean_name, "_reclassification_summary.csv")
    )
  }
}

if (length(hfref_reclass) > 0) {
  hfref_plots <- create_reclassification_plots(hfref_reclass, "HFrEF")

  # Save plots for each miRNA
  for (mir_name in names(hfref_plots)) {
    # Clean miRNA name for filename
    clean_name <- gsub("[^A-Za-z0-9]", "_", mir_name)

    ggsave(
      filename = paste0(
        "figures/hfref_",
        clean_name,
        "_reclassification_summary.png"
      ),
      plot = hfref_plots[[mir_name]]$summary,
      width = 10,
      height = 6,
      dpi = 300
    )

    ggsave(
      filename = paste0(
        "figures/hfref_",
        clean_name,
        "_probability_distribution.png"
      ),
      plot = hfref_plots[[mir_name]]$distribution,
      width = 10,
      height = 6,
      dpi = 300
    )

    ggsave(
      filename = paste0(
        "figures/hfref_",
        clean_name,
        "_scatter_reclassification.png"
      ),
      plot = hfref_plots[[mir_name]]$scatter,
      width = 12,
      height = 8,
      dpi = 300
    )
  }

  # Save reclassification data
  for (mir_name in names(hfref_reclass)) {
    clean_name <- gsub("[^A-Za-z0-9]", "_", mir_name)
    write_csv(
      hfref_reclass[[mir_name]]$data,
      paste0("output/hfref_", clean_name, "_reclassification_data.csv")
    )
    write_csv(
      hfref_reclass[[mir_name]]$summary,
      paste0("output/hfref_", clean_name, "_reclassification_summary.csv")
    )
  }
}


saveRDS(dcm_cleaned, "revision2025/data/dcm_cleaned_data.rds")
saveRDS(hfref_cleaned, "revision2025/data/hfref_cleaned_data.rds")
