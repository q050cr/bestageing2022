library(tidyverse)
library(readxl)
library(glue)

data_path_bestageing2022 <- "/mnt/nas185/reich/rockerprojects/bestageing2022"
data_path_BestAgeing <- "/mnt/nas185/reich/BestAgeing"

# Load processed data
data_processed_acs <- readRDS(glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_acs_data01.rds"))
data_processed_cad <- readRDS(glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_cad_data01.rds"))

# Load metadata
meta_acs <- read_excel(glue("{data_path_BestAgeing}/data/pheno_acs.xlsx")) %>%
  mutate(
    BestAgeingCode,
    phenodat = "acs",
    ntprobnp = as.numeric(`ACS-Lab-ProBNP`),
    crea = as.numeric(`ACS-Lab-Cr`),
    .keep = "none"
  )

meta_cad <- read_excel(glue("{data_path_BestAgeing}/data/pheno_cad.xlsx")) %>%
  mutate(
    BestAgeingCode,
    phenodat = "cad",
    ntprobnp = as.numeric(`CAD-Lab-ProBNP`),
    crea = as.numeric(`CAD-Lab-Cr`),
    .keep = "none"
  )

meta_control <- read_excel(glue("{data_path_BestAgeing}/data/pheno_controls.xlsx")) %>%
  mutate(
    BestAgeingCode,
    phenodat = "control",
    ntprobnp = as.numeric(`Controls-Lab-ProBNP`),
    crea = as.numeric(`Controls-Lab-Cr`),
    .keep = "none"
  )

# Function to clean data
clean_ntprobnp_data <- function(data, disease_name) {
  cat("Initial data dimensions:", dim(data), "\n")

  # Remove controls with NT-proBNP > 450
  controls_high_bnp <- data$disease == "control" & !is.na(data$ntprobnp) & data$ntprobnp > 450
  n_removed_controls <- sum(controls_high_bnp, na.rm = TRUE)
  if (n_removed_controls > 0) {
    data <- data[!controls_high_bnp, ]
    cat("Removed", n_removed_controls, "controls with NT-proBNP > 450\n")
  }

  # Mean imputation for controls
  control_ntprobnp_mean <- mean(data$ntprobnp[data$disease == "control"], na.rm = TRUE)
  data$ntprobnp[data$disease == "control" & is.na(data$ntprobnp)] <- control_ntprobnp_mean
  attr(data, "control_imputation_value") <- control_ntprobnp_mean

  # Create GLM data (remove all with missing NT-proBNP)
  data_glm <- data[!is.na(data$ntprobnp), ]

  return(list(full_data = data, glm_data = data_glm))
}

# Merge and clean ACS
acs_merged <- data_processed_acs %>%
  left_join(bind_rows(meta_acs, meta_control), by = c("pat_id" = "BestAgeingCode")) %>%
  mutate(preprocessing = "acs_pipeline")
acs_cleaned <- clean_ntprobnp_data(acs_merged, "acs")

# Merge and clean CAD
cad_merged <- data_processed_cad %>%
  left_join(bind_rows(meta_cad, meta_control), by = c("pat_id" = "BestAgeingCode")) %>%
  mutate(preprocessing = "cad_pipeline")
cad_cleaned <- clean_ntprobnp_data(cad_merged, "cad")

# Save
saveRDS(acs_cleaned, glue("{data_path_bestageing2022}/revision2025/data/acs_cleaned_data.rds"))
saveRDS(cad_cleaned, glue("{data_path_bestageing2022}/revision2025/data/cad_cleaned_data.rds"))

cat("ACS and CAD cleaned data saved successfully\n")
cat("ACS GLM data dimensions:", dim(acs_cleaned$glm_data), "\n")
cat("CAD GLM data dimensions:", dim(cad_cleaned$glm_data), "\n")
