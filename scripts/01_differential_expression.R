# Differential miRNA Expression Analysis
#
# This script performs differential expression analysis on miRNA data
# Author: Christoph Reich
# Date: 2025-05-20 (updated)

# ----- Setup -----
# Load configuration and utilities
source("scripts/config/config.R")
source("scripts/helper/utils.R")

# Get paths
paths <- get_project_paths()

# Load required libraries
load_project_libraries(c(
  "tidyverse",
  "limma",
  "edgeR",
  "ggplot2",
  "readxl",
  "writexl",
  "pROC"
))

# ----- Parameters -----
# Analysis parameters - edit these as needed
DISEASE <- "CAD" # Choose from available diseases
MATCHED <- TRUE # Use matched samples?
NORMALIZE <- TRUE # Normalize expression data?
FILTER_LOW_EXPRESSED <- TRUE # Filter lowly expressed miRNAs?

# ----- Data Import -----
# Import data
clinical_data <- import_project_data("clinical")
mirna_data <- import_project_data("mirna")
metadata <- import_project_data("metadata")

# ----- Data Preparation -----
# Prepare data for analysis
prepare_data <- function(clinical_data, mirna_data, metadata, disease, matched = FALSE) {
  # Join clinical data with miRNA data
  merged_data <- inner_join(clinical_data, mirna_data, by = "patient_id")

  # Filter by disease
  disease_data <- merged_data %>%
    filter(Disease == disease)

  # Match samples if required
  if (matched) {
    # Create matched groups based on age, sex, etc.
    matched_data <- matchIt::matchIt(
      Disease ~ Age + Sex + BMI,
      data = disease_data,
      method = "nearest",
      ratio = 1
    )

    # Extract matched samples
    disease_data <- matchIt::match.data(matched_data)
  }

  return(disease_data)
}

# Prepare data
analysis_data <- prepare_data(
  clinical_data = clinical_data,
  mirna_data = mirna_data,
  metadata = metadata,
  disease = DISEASE,
  matched = MATCHED
)

# ----- Differential Expression Analysis -----
# Perform differential expression analysis
run_de_analysis <- function(analysis_data, normalize = TRUE, filter_low = TRUE) {
  # Extract expression matrix
  expr_matrix <- analysis_data %>%
    select(starts_with("hsa-")) %>%
    as.matrix()

  # Create design matrix
  design <- model.matrix(~Disease, data = analysis_data)

  # Create DGEList object
  dge <- DGEList(counts = expr_matrix)

  # Filter low expressed miRNAs if required
  if (filter_low) {
    keep <- filterByExpr(dge, design)
    dge <- dge[keep, ]
  }

  # Normalize if required
  if (normalize) {
    dge <- calcNormFactors(dge)
  }

  # Fit model
  fit <- lmFit(dge, design)
  fit <- eBayes(fit)

  # Get results
  results <- topTable(fit, coef = 2, number = Inf)

  return(list(
    dge = dge,
    fit = fit,
    results = results
  ))
}

# Run DE analysis
de_results <- run_de_analysis(
  analysis_data = analysis_data,
  normalize = NORMALIZE,
  filter_low = FILTER_LOW_EXPRESSED
)

# ----- Visualization -----
# Create volcano plot
create_volcano_plot <- function(de_results, title = "Volcano Plot", fdr_threshold = 0.05, fc_threshold = 1) {
  results_df <- de_results$results %>%
    rownames_to_column("miRNA") %>%
    mutate(
      significant = adj.P.Val < fdr_threshold & abs(logFC) > fc_threshold,
      direction = case_when(
        logFC > fc_threshold & adj.P.Val < fdr_threshold ~ "Up",
        logFC < -fc_threshold & adj.P.Val < fdr_threshold ~ "Down",
        TRUE ~ "Not Significant"
      )
    )

  # Create plot
  volcano_plot <- ggplot(results_df, aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("Down" = "blue", "Up" = "red", "Not Significant" = "gray")) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed") +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed") +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-log10(FDR)",
      color = "Regulation"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  return(volcano_plot)
}

# Create volcano plot
volcano_plot <- create_volcano_plot(
  de_results = de_results,
  title = paste("Differential Expression:", DISEASE)
)

# Save plot
if (SAVE_FILES) {
  output_dir <- file.path(paths$output_path, "de_results")
  ensure_dir_exists(output_dir)

  # Save plot
  ggsave(
    file.path(output_dir, paste0(format(Sys.Date(), "%Y%m%d"), "_volcano_plot_", DISEASE, ".png")),
    volcano_plot,
    width = 10,
    height = 8,
    dpi = 300
  )

  # Save DE results
  save_project_object(
    de_results$results %>% rownames_to_column("miRNA"),
    name = paste0("de_results_", DISEASE),
    type = ifelse(MATCHED, "matched", "all"),
    extension = "csv",
    dir = output_dir
  )
}

# Print top DE miRNAs
top_de <- de_results$results %>%
  rownames_to_column("miRNA") %>%
  arrange(adj.P.Val) %>%
  head(20)

print(top_de)

# ----- Session Info -----
# Print session information for reproducibility
print(sessionInfo())
