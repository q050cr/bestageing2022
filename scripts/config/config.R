# Global configuration file for BestAgeing2022 project
# This file centralizes all path definitions and common settings

# ----- System Configuration -----
# Get system name
system_name <- Sys.info()["nodename"]

# ----- Project Settings -----
PROJECT_NAME <- "bestageing2022"
SAVE_FILES <- TRUE
RUN_TESTS <- TRUE

# ----- Path Configuration -----
# Define function to get appropriate paths
get_project_paths <- function(mount_filesystem = TRUE) {
  # Default paths for cluster environment
  paths <- list(
    lib_path = "/mnt/users/reich/programs/R43/lib",
    project_path = "/mnt/users/reich/rockerprojects/bestageing2022",
    data_path_BestAgeing = "/mnt/users/reich/BestAgeing"
  )

  # Override paths for MacBook
  if (
    grepl(
      "MacBook-Pro-CR-2065|dhcp172-619.laptop-zim.uni-heidelberg.de",
      system_name
    )
  ) {
    if (mount_filesystem) {
      # Mounted filesystem paths
      paths$lib_path <- .libPaths()[1]
      paths$project_path <- "/Users/christophreich/Desktop/mount/rockerprojects/bestageing2022"
      paths$data_path_BestAgeing <- "/Users/christophreich/Desktop/mount/BestAgeing"
    } else {
      # Local paths
      paths$lib_path <- .libPaths()[1]
      paths$project_path <- "/Volumes/T7CR/data/bestageing2022"
      paths$data_path_BestAgeing <- "/Volumes/T7CR/data/BestAgeing"
    }
  }

  # Add derived paths
  paths$scripts_path <- file.path(paths$project_path, "scripts")
  paths$data_path <- file.path(paths$project_path, "data")
  paths$output_path <- file.path(paths$project_path, "output")
  paths$figures_path <- file.path(paths$project_path, "figures")

  return(paths)
}

# ----- Library Management -----
# Function to load required libraries
load_project_libraries <- function(required_packages) {
  paths <- get_project_paths()

  # Set library path
  .libPaths(paths$lib_path)

  # Check if packages are installed, install if needed
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message(paste("Installing package:", pkg))
      install.packages(pkg, lib = paths$lib_path)
      library(pkg, character.only = TRUE)
    }
  }
}

# ----- Common Functions -----
# Function to create directories if they don't exist
ensure_dir_exists <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message(paste("Created directory:", dir_path))
  }
}

# ----- Usage Example -----
# To use this configuration in your scripts:
# source("/mnt/users/reich/rockerprojects/bestageing2022/scripts/config/config.R")
# paths <- get_project_paths()
# load_project_libraries(c("tidyverse", "caret", "survival"))
