# Utility functions for the BestAgeing2022 project

# Load configuration
source(file.path(dirname(dirname(getwd())), "scripts/config/config.R"))
paths <- get_project_paths()

# ----- Data Processing Functions -----

#' Check if two dataframes are identical
#'
#' @param df1 First dataframe
#' @param df2 Second dataframe
#' @param verbose Print detailed differences if TRUE
#' @return Logical value indicating if dataframes are identical
check_identical_dataframes <- function(df1, df2, verbose = TRUE) {
  source(file.path(paths$scripts_path, "helper/check_identical_dataframes.R"))
  check_identical_dataframes_impl(df1, df2, verbose)
}

#' Save an R object with standardized naming
#'
#' @param object Object to save
#' @param name Base name for the file
#' @param type Type of analysis (e.g., "DE", "ML", "Survival")
#' @param extension File extension (default: "rds")
#' @param dir Directory to save in (default: output_path)
#' @return Path to the saved file
save_project_object <- function(object, name, type, extension = "rds", dir = paths$output_path) {
  timestamp <- format(Sys.time(), "%Y%m%d")
  filename <- paste0(timestamp, "_", name, "_", type, ".", extension)
  filepath <- file.path(dir, filename)

  ensure_dir_exists(dirname(filepath))

  if (extension == "rds") {
    saveRDS(object, filepath)
  } else if (extension == "csv") {
    write.csv(object, filepath, row.names = FALSE)
  } else {
    save(object, file = filepath)
  }

  message(paste("Saved:", filepath))
  return(filepath)
}

# ----- Plotting Helpers -----

#' Get consistent ggplot theme for project
#'
#' @param base_size Base font size
#' @param base_family Base font family
#' @return A ggplot theme object
get_project_theme <- function(base_size = 12, base_family = "") {
  # Load the original theme function
  source(file.path(paths$scripts_path, "helper/custom_ggplot_theme.R"))

  # Create theme with consistent styling
  theme_ba <- function() {
    theme_bw(base_size = base_size, base_family = base_family) +
      theme(
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom"
      )
  }

  return(theme_ba())
}

#' Create custom survival plot
#'
#' @param fit Survival fit object
#' @param ... Additional parameters for ggsurvplot
#' @return A survival plot
custom_survival_plot <- function(fit, ...) {
  source(file.path(paths$scripts_path, "helper/ggsurvplot_custom.R"))
  ggsurvplot_custom(fit, ...)
}

#' Create custom variable importance plot
#'
#' @param model Model object
#' @param top_n Number of top features to show
#' @return A VIP plot
custom_vip_plot <- function(model, top_n = 10) {
  source(file.path(paths$scripts_path, "helper/custom_vip_plot.R"))
  custom_vip_plot_impl(model, top_n)
}

# ----- Model Performance Functions -----

#' Create custom calibration plot
#'
#' @param obs Observed values
#' @param pred Predicted values
#' @param title Plot title
#' @return A calibration plot
custom_calibration_plot <- function(obs, pred, title = "Calibration Plot") {
  source(file.path(paths$scripts_path, "helper/custom_calibration_plot.R"))
  custom_calibration_plot_impl(obs, pred, title)
}

#' Create custom performance summary barplot
#'
#' @param perf_data Performance data
#' @param metrics Metrics to include
#' @return A performance summary barplot
custom_performance_summary_barplot <- function(perf_data, metrics = c("AUC", "Accuracy", "Sensitivity", "Specificity")) {
  source(file.path(paths$scripts_path, "helper/custom_performance_summary_barplot.R"))
  custom_performance_summary_barplot_impl(perf_data, metrics)
}

# ----- Data Import Functions -----

#' Import standardized project data
#'
#' @param data_type Type of data to import
#' @param version Version of the data (optional)
#' @return Imported data
import_project_data <- function(data_type, version = NULL) {
  valid_types <- c("clinical", "mirna", "metadata", "matched")

  if (!(data_type %in% valid_types)) {
    stop(paste("Invalid data type. Must be one of:", paste(valid_types, collapse = ", ")))
  }

  # Default patterns for different data types
  patterns <- list(
    clinical = "Patientenliste.*\\.xlsx$",
    mirna = ".*miRNA.*\\.rds$",
    metadata = ".*annotation.*\\.xlsx$",
    matched = ".*matched.*\\.rds$"
  )

  # Find the latest file matching the pattern if version is not specified
  if (is.null(version)) {
    files <- list.files(paths$data_path,
      pattern = patterns[[data_type]],
      full.names = TRUE, recursive = TRUE
    )

    if (length(files) == 0) {
      stop(paste("No files found matching the pattern for data type:", data_type))
    }

    # Sort by modification time and get the most recent
    file_info <- file.info(files)
    file_info$name <- rownames(file_info)
    file_info <- file_info[order(file_info$mtime, decreasing = TRUE), ]
    file_path <- file_info$name[1]
  } else {
    # Use specified version
    if (data_type == "clinical") {
      file_path <- file.path(paths$data_path, paste0(version, "_Patientenliste_BestAgeing_Master_clean.xlsx"))
    } else if (data_type == "mirna") {
      file_path <- file.path(paths$data_path, paste0(version, "_XMELD_abfrage_best_ageing.rds"))
    } else if (data_type == "metadata") {
      file_path <- file.path(paths$data_path, paste0(version, "_annotation_chrstoph_reich.xlsx"))
    } else if (data_type == "matched") {
      file_path <- file.path(paths$data_path, paste0(version, "_matched_data.rds"))
    }
  }

  if (!file.exists(file_path)) {
    stop(paste("File does not exist:", file_path))
  }

  # Import based on file extension
  if (grepl("\\.rds$", file_path)) {
    data <- readRDS(file_path)
  } else if (grepl("\\.xlsx$", file_path)) {
    data <- readxl::read_excel(file_path)
  } else if (grepl("\\.csv$", file_path)) {
    data <- read.csv(file_path)
  } else {
    stop("Unsupported file format")
  }

  message(paste("Imported data from:", file_path))
  return(data)
}
