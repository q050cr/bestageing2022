

version$version.string

#  Load libraries
library(miRetrieve); packageVersion('miRetrieve')
library(magrittr) # Load magrittr for %>%
library(ggplot2) # Load ggplot2 for plotting
library(glue)
library(rcrossref)
library(dplyr) # Data wrangling
library(tidyr) # Data wrangling

saveFILE <- TRUE

##
# QUERIES ---------
##
diseases <- c("ACS", "CAD", "DCM", "HFrEF")

# Use the vectorized glue function to generate a vector of strings
path2files <- glue::glue_data(.x = data.frame(diseases = diseases),
                             "./data-literature/miRetrieve/{diseases}/2023-06-13-human-df_count_both.rds")
path2files.bm <- glue::glue_data(.x = data.frame(diseases = diseases),
                                 "./data-literature/miRetrieve/{diseases}/2023-06-13-human-disease_biomarker.rds")

for (data in seq_along(path2files)) {
  df_count_both <- readRDS(path2files[data])
  disease_biomarker <- readRDS(path2files.bm[data])  # also has biomarker_score
  
}

