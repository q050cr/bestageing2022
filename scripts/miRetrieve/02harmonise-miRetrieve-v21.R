

# Get system name
system_name <- Sys.info()["nodename"]
mount_filesystem <- TRUE
saveFILE <- TRUE

# Define library and data paths based on system
if (system_name == "MacBook-Pro-CR-2065.local") {
  lib_path <- .libPaths()[1]
  data_path_bestageing2022 <- "/Volumes/T7CR/data/bestageing2022"
  data_path_BestAgeing <- "/Volumes/T7CR/data/BestAgeing"
  if(mount_filesystem == TRUE) {
    data_path_bestageing2022 <- "/Users/christophreich/Desktop/mount/rockerprojects/bestageing2022"  # mount -t nfs 10.55.1.185:/data/users/reich/ ~/Desktop/mount/
    data_path_BestAgeing <- "/Users/christophreich/Desktop/mount/BestAgeing"
  }
} else {  # assuming cluster
  lib_path <- "/mnt/users/reich/programs/R43/lib" 
  data_path_bestageing2022 <- "/mnt/users/reich/rockerprojects/bestageing2022"
  data_path_BestAgeing <- "/mnt/users/reich/BestAgeing"
}

version$version.string

#  Load libraries
require(miRetrieve, lib.loc = lib_path); packageVersion('miRetrieve')
require(magrittr, lib.loc = lib_path) # Load magrittr for %>%
require(ggplot2, lib.loc = lib_path) # Load ggplot2 for plotting
require(glue, lib.loc = lib_path)
require(rcrossref, lib.loc = lib_path)
require(dplyr, lib.loc = lib_path) # Data wrangling
require(tidyr, lib.loc = lib_path) # Data wrangling
require(miRBaseConverter, lib.loc = lib_path)
require(janitor, lib.loc = lib_path)

# file oi created in "01miRetrieve_loop.R" and stored in "./data-literature/miRetrieve/ACS/2023-06-14-human-disease_biomarker.rds" 
# for biomarker relevant miRNA and converted to target v21

##
# Harmonize literature miRNAs to v21 (best version in bestageing mirna naming) -----------------------------
##
diseases <- c("ACS", "CAD", "DCM", "HFrEF")
subset <- "-human" # changed to "-human" only, could also be changed to "-humantrials" only, or "" for full pubmed search !!!

load(file = glue('{data_path_BestAgeing}/data/data.rda'))  # "UKL-HD" n=731
# get some miRNAs from data's colnames
miRNANames <- colnames(data[10:2558])
length(miRNANames)
version_mir = miRBaseConverter::checkMiRNAVersion(miRNANames, verbose = TRUE) 
result1 = miRBaseConverter::miRNA_NameToAccession(miRNANames,version = c("v21"))  # get Accession MIMAT (dentifier that define miRNA uniquely in miRBase)
result1[1:10,]

# Use the vectorized glue function to generate a vector of strings
path2files <- glue::glue_data(.x = data.frame(diseases = diseases),
                             "{data_path_bestageing2022}/data-literature/miRetrieve/{tolower(diseases)}/2023-07-27", {subset}, "-df_count_both.rds")  # changed from 2023-06-13
path2files.bm <- glue::glue_data(.x = data.frame(diseases = diseases),
                                 "{data_path_bestageing2022}/data-literature/miRetrieve/{tolower(diseases)}/2023-07-27", {subset}, "-disease_biomarker.rds")

for (mydata in seq_along(path2files)) {
  df_count_both <- readRDS(path2files[mydata])
  disease_biomarker <- readRDS(path2files.bm[mydata])  # also has biomarker_score
  
  ## fast conversion ----------------------------------------------------------
  convert.fast = miRBaseConverter::miRNAVersionConvert(paste0("hsa-", df_count_both$miRNA),targetVersion = "v21",exact = TRUE, verbose = TRUE) # default exact and verbose are TRUE
  convert.fast.bm = miRBaseConverter::miRNAVersionConvert(paste0("hsa-", unique(disease_biomarker$miRNA)),targetVersion = "v21",exact = TRUE, verbose = TRUE)  # also better than manual conversion  
  
  convert.fast <- convert.fast %>%
    separate_rows(TargetName, Accession, sep = "&")
  convert.fast.bm <- convert.fast.bm %>%
    separate_rows(TargetName, Accession, sep = "&")
  
  df_count_both_accession <- df_count_both %>% 
    mutate(miRNA = paste0("hsa-", miRNA)) %>% 
    left_join(convert.fast, by = c("miRNA"="OriginalName"))
  
  disease_biomarker_accession <- disease_biomarker %>% 
    mutate(miRNA = paste0("hsa-", miRNA)) %>% 
    left_join(convert.fast.bm, by = c("miRNA"="OriginalName"), relationship = "many-to-many")
  
  # SAVE
  if (saveFILE == TRUE) {
    # df_count_both_accession
    saveRDS(object = df_count_both_accession, 
            file = glue::glue("{data_path_bestageing2022}/data-literature/miRetrieve/{tolower(diseases[mydata])}/{Sys.Date()}{subset}-df_count_both_with_accession.rds"))
    #write.csv2(x = df_count_both_accession, file = glue::glue("{data_path_bestageing2022}/data-literature/miRetrieve/{diseases[mydata]}/{Sys.Date()}{subset}-df_count_both_with_accession.csv"))
    # disease_biomarker_with_accession
    saveRDS(object = disease_biomarker_accession, 
            file = glue::glue("{data_path_bestageing2022}/data-literature/miRetrieve/{tolower(diseases[mydata])}/{Sys.Date()}{subset}-disease_biomarker_with_accession.rds"))
    #write.csv2(x = disease_biomarker_accession, file = glue::glue("{data_path_bestageing2022}/data-literature/miRetrieve/{diseases[mydata]}/{Sys.Date()}{subset}-disease_biomarker_with_accession.csv"))
  }
  
  ## manually match to most frequent version -----------------------------------
  version_mir = miRBaseConverter::checkMiRNAVersion(paste0("hsa-", df_count_both$miRNA), verbose = TRUE) 
  result1 = miRBaseConverter::miRNA_NameToAccession(paste0("hsa-", df_count_both$miRNA),version = version_mir) 
  # match to second most frequent version
  index <- which(is.na(result1[2]))
  version2_mir = miRBaseConverter::checkMiRNAVersion(paste0("hsa-", df_count_both$miRNA[index]), verbose = TRUE) 
  result2 = miRBaseConverter::miRNA_NameToAccession(paste0("hsa-", df_count_both$miRNA[index]),version = version2_mir)
  # combine to result1
  result1$Accession[index] <- result2$Accession
  
  # miRBase Accessions to miRNA Names of the target version
  convert.manual = miRBaseConverter::miRNA_AccessionToName(result1[,2],targetVersion = "v21")
  
  print(glue("|||----------------------- Run finished for disease: ", {diseases[mydata]} , " -----------------------|||"))
}


