### INFO ----------------------------------------------------------------------
# Differential miRNA Expression Analysis with Detection Matrix Filtering
# Author: Christoph Reich
# Date: 2025-05-20
#
# This script performs differential expression analysis with detection matrix filtering
# It creates plots: "fig01vogel2013", "fig02vogel2013"

# Load configuration
source("scripts/config/config.R")

# Get project paths
paths <- get_project_paths()
data_path_bestageing2022 <- paths$project_path
data_path_BestAgeing <- paths$data_path_BestAgeing
lib_path <- paths$lib_path

# Global settings
SAVE.files <- TRUE
runTests <- TRUE
runqqnorm <- FALSE # already done before
filterDetMatrix <- TRUE # Use detection matrix filtering

# Helper functions
convert_mir_name <- function(name) {
  # Replace 'mir' with 'miR'
  name <- gsub("mir", "miR", name)

  # Replace underscores with hyphens
  name <- gsub("_", "-", name)

  return(name)
}
convert_mir_name_V <- Vectorize(convert_mir_name)

# Load required libraries
required_packages <- c(
  # Data manipulation
  "readxl", "janitor", "glue", "dplyr", "tidyr", "tibble", "stringr", "purrr",
  # Visualization
  "ggplot2", "ggrepel", "ggforce", "pheatmap", "corrplot",
  # Analysis
  "limma", "edgeR", "sva", "WGCNA", "dendextend",
  # Tables
  "flextable", "officer"
)

# Load packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, lib = lib_path)
    library(pkg, character.only = TRUE)
  }
}

# Load helper functions
source(file.path(data_path_bestageing2022, "scripts/helper/custom_ggplot_theme.R"))

#-----------------------------------------------------------------------------
# Load data
#-----------------------------------------------------------------------------

# Load miRNA data
load(file = glue("{data_path_BestAgeing}/data/mirnas.rda")) # "UKL-HD" n=765
load(file = glue("{data_path_BestAgeing}/data/data.rda")) # "UKL-HD" n=731

# Clean miRNA data
all_mirnas <- clean_names(mirnas)
names(all_mirnas) <- str_replace(string = names(all_mirnas), pattern = "mi_r", replacement = "mir")
mirnas_disease <- clean_names(data)
names(mirnas_disease) <- str_replace(string = names(mirnas_disease), pattern = "mi_r", replacement = "mir")
rm(data, mirnas)

# Load sequencing metadata
hbdx_metadat <- clean_names(readxl::read_excel(path = glue("{data_path_bestageing2022}/data/kahraman2023/230907_annotation_chrstoph_reich.xlsx")))
rin_mean <- mean(as.numeric(hbdx_metadat$rin), na.rm = TRUE)
rin_sd <- sd(as.numeric(hbdx_metadat$rin), na.rm = TRUE)

# Load detection matrix
det_mat_all_mirnas <- read.table(glue("{data_path_bestageing2022}/data/kahraman2023/det_mat_all_mirnas.txt")) %>%
  t() %>%
  as.data.frame() %>%
  clean_names()
colnames(det_mat_all_mirnas) <- str_replace(string = colnames(det_mat_all_mirnas), pattern = "mi_r", replacement = "mir")
rownames(det_mat_all_mirnas) <- gsub("\\.", "-", rownames(det_mat_all_mirnas))
det_mat_all_mirnas <- det_mat_all_mirnas %>%
  rownames_to_column(var = "pat_id") %>%
  as_tibble()

# Process metadata and handle duplicates
unique_slide_ids <- length(unique(hbdx_metadat$slide_id))
arrays_per_slide <- nrow(hbdx_metadat) / unique_slide_ids
duplicate_ids <- hbdx_metadat[duplicated(hbdx_metadat$customer_id), ] %>% pull(customer_id)
hbdx_metadat <- hbdx_metadat %>% distinct(customer_id, .keep_all = TRUE)

# Load project metadata
load(glue("{data_path_BestAgeing}/data/clean_all_meta.rda"))
clean_all_meta <- clean_all_meta %>% mutate(age = ifelse(age < 18, NA, age))

# Load control and disease IDs
control_ids <- read_excel(glue("{data_path_BestAgeing}/data/pheno_controls.xlsx")) %>%
  dplyr::pull(BestAgeingCode)
control_ids <- control_ids[control_ids != "UKL-HD-00318"] # In CAD dataset

# Define diseases and analysis types
diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>%
  filter(analysis == "selected")

#-----------------------------------------------------------------------------
# Main analysis loop
#-----------------------------------------------------------------------------

for (i in 1:nrow(all_combis)) {
  # Display info header
  cat("\n\n")
  print(glue("|||------------------------------------------------------------------------------------------------------------|||"))
  print(glue(
    "|||-----------------------Start DE analysis for disease: ", stringr::str_to_upper(all_combis$diseases[i]),
    ", and selection of miRNAs: ", stringr::str_to_upper(all_combis$analysis[i]), " -----------------------|||"
  ))
  print(glue("|||------------------------------------------------------------------------------------------------------------|||"))
  cat("\n")

  # Load disease-specific IDs
  filename <- glue("{data_path_BestAgeing}/data/pheno_{all_combis$diseases[i]}.xlsx")
  disease_vector <- paste0(all_combis$diseases[i], "_ids")
  assign(
    x = disease_vector,
    value = read_excel(filename) %>%
      dplyr::pull(BestAgeingCode)
  )

  # Create disease identification dataframe
  disease_ident_df <- tibble(
    pat_id = c(eval(as.symbol(disease_vector)), control_ids),
    disease = c(
      rep(all_combis$diseases[i], length(eval(as.symbol(disease_vector)))),
      rep("control", length(control_ids))
    )
  )

  # Prepare data for analysis
  data01 <- all_mirnas %>%
    # Filter patients
    filter(pat_id %in% !!sym(disease_vector) | pat_id %in% control_ids) %>%
    # Join metadata
    left_join(clean_all_meta %>% select(patID, age, sex), by = c("pat_id" = "patID")) %>%
    left_join(disease_ident_df, by = "pat_id") %>%
    relocate(c(disease, age, sex), .after = pat_id) %>%
    mutate(
      disease = factor(disease),
      sex = factor(sex)
    ) %>%
    as_tibble()

  # Join detection matrix
  data01 <- data01 %>%
    left_join(det_mat_all_mirnas, by = "pat_id")

  # Separate detection matrix columns
  data01_miRNAs <- data01 %>%
    select(starts_with("hsa_mir"))

  # Count miRNAs before filtering
  n_mirnas_before <- ncol(data01_miRNAs) / 2

  # Get column names for miRNAs and detection values
  miRNA_cols <- colnames(data01_miRNAs)[1:n_mirnas_before]
  detection_cols <- colnames(data01_miRNAs)[(n_mirnas_before + 1):(2 * n_mirnas_before)]

  # Apply detection matrix filtering if enabled
  if (filterDetMatrix) {
    # Filter miRNAs based on detection threshold
    detection_matrix <- data01 %>% select(all_of(detection_cols))

    # Calculate proportion of samples where each miRNA is detected
    detection_prop <- colMeans(detection_matrix == 1, na.rm = TRUE)

    # Keep miRNAs detected in at least 50% of samples
    miRNAs_to_keep <- names(detection_prop[detection_prop >= 0.5])
    miRNAs_to_keep <- gsub("\\.detection$", "", miRNAs_to_keep)

    # Filter data to keep only well-detected miRNAs
    data01_filtered <- data01 %>%
      select(pat_id, disease, age, sex, all_of(miRNAs_to_keep))

    # Save filtered data
    cat("Filtered", n_mirnas_before - length(miRNAs_to_keep), "of", n_mirnas_before, "miRNAs based on detection matrix\n")
  } else {
    # No filtering, keep all miRNAs
    data01_filtered <- data01 %>%
      select(pat_id, disease, age, sex, all_of(miRNA_cols))
  }

  # Perform differential expression analysis
  # Set up design matrix
  design <- model.matrix(~ 0 + disease + sex + age, data = data01_filtered)
  colnames(design) <- gsub("disease", "", colnames(design))

  # Create contrast for disease vs control
  contrast <- makeContrasts(
    contrasts = paste0(all_combis$diseases[i], " - control"),
    levels = design
  )

  # Create DGEList object
  dge <- DGEList(
    counts = data01_filtered %>% select(all_of(miRNAs_to_keep)) %>% as.matrix(),
    genes = data.frame(miRNA = miRNAs_to_keep)
  )

  # Filter by expression
  keep <- filterByExpr(dge, design)
  dge <- dge[keep, ]

  # Normalize
  dge <- calcNormFactors(dge)

  # Fit model
  fit <- lmFit(dge, design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)

  # Get results
  tt <- topTable(fit2, n = Inf)

  # Add gene names and annotations
  results <- tt %>%
    rownames_to_column("miRNA") %>%
    mutate(
      miRNA_formatted = convert_mir_name_V(miRNA),
      significant = adj.P.Val < 0.05,
      regulation = case_when(
        logFC > 0 & adj.P.Val < 0.05 ~ "Upregulated",
        logFC < 0 & adj.P.Val < 0.05 ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )

  # Create volcano plot
  volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = regulation)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = paste("Volcano Plot -", toupper(all_combis$diseases[i])),
      x = "Log2 Fold Change",
      y = "-log10(FDR)"
    ) +
    theme_bw() +
    my_base_theme() +
    theme(legend.position = "bottom")

  # Label top differentially expressed miRNAs
  top_de <- results %>%
    filter(significant == TRUE) %>%
    arrange(adj.P.Val) %>%
    head(10)

  volcano_plot_labeled <- volcano_plot +
    geom_text_repel(
      data = top_de,
      aes(label = miRNA_formatted),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.2,
      segment.color = "grey50"
    )

  # Save results
  if (SAVE.files) {
    # Create output directory if it doesn't exist
    output_dir <- file.path(data_path_bestageing2022, "output", "de_results", all_combis$diseases[i])
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    # Save processed data
    saveRDS(
      data01_filtered,
      file = glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_{all_combis$diseases[i]}_data01.rds")
    )

    # Save DE results
    saveRDS(
      results,
      file = glue("{output_dir}/001c_de_results_{all_combis$diseases[i]}_DET_MATRIX.rds")
    )

    # Save volcano plot
    ggsave(
      file = glue("{output_dir}/001c_volcano_plot_{all_combis$diseases[i]}_DET_MATRIX.png"),
      plot = volcano_plot_labeled,
      width = 10,
      height = 8,
      dpi = 300
    )

    # Save top DE results as CSV
    write.csv(
      results %>% arrange(adj.P.Val) %>% head(100),
      file = glue("{output_dir}/001c_top_de_results_{all_combis$diseases[i]}_DET_MATRIX.csv"),
      row.names = FALSE
    )
  }

  # Print summary
  cat("\nDifferential Expression Analysis Summary:\n")
  cat("Disease:", toupper(all_combis$diseases[i]), "\n")
  cat("Total miRNAs analyzed:", nrow(results), "\n")
  cat("Significant miRNAs (FDR < 0.05):", sum(results$significant), "\n")
  cat("Upregulated:", sum(results$regulation == "Upregulated"), "\n")
  cat("Downregulated:", sum(results$regulation == "Downregulated"), "\n\n")

  # Print completion message
  print(paste0(
    "|||-----------------------DE analysis finished for disease: ", stringr::str_to_upper(all_combis$diseases[i]),
    " -----------------------|||"
  ))
} # END LOOP
