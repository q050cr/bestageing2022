

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

require(glue)
require(miRetrieve)



diseases <- c("ACS", "CAD", "DCM", "HFrEF")
diseases_tibble <- tibble(abbrev_disease = diseases,
                          disease_name = c("Acute Coronary Syndrome", "Coronary Artery Disease", "Dilatative Cardiomyopathy", "Heart Failure"))
path2dat <- paste0(data_path_bestageing2022, "/data-literature/pubmed/", list.files(path = glue("{data_path_bestageing2022}/data-literature/pubmed/")))
path2dat <- path2dat[grep("-2023-06-16-", path2dat)]   # use newest files from 2023-06-16


for(disease in seq_along(diseases)) {
  # subset research!
  subset <- "-human" # changed to "-human" only, could also be changed to "-humantrials" only, or "" for full pubmed search !!!, only works without date in next line!
  file <- glue::glue("{data_path_bestageing2022}/data-literature/pubmed-downloads/2023-06-14-pubmed-MicroRNA{subset}-{toupper(diseases[disease])}.txt")  
  
  ## load abstracts ----------------------------------------------------------
  # pubmed read only
  pubmed_read_only_df <- read_pubmed(file, topic = diseases[disease])
  print(glue("For {diseases[disease]} we yielded {length(unique(pubmed_read_only_df$PMID))} unique abstracts on PubMed."))
  print("-------------- ")
  
  no_research_subset_df <- read_pubmed(file, topic = diseases[disease])  %>% 
    # Subset for original articles
    #subset_research() %>% 
    extract_mir_df(threshold = 1,
                   extract_letters = TRUE)
  print(glue("Abstracts: For {diseases[disease]} without subset RESEARCH {length(unique(no_research_subset_df$PMID))} unique abstracts were found."))
  print(glue("miRNAs: For {diseases[disease]} without subset RESEARCH {length(unique(no_research_subset_df$miRNA))} unique miRNAs were found."))
  print("-------------- ")
 
  
  subset_df <- read_pubmed(file, topic = diseases[disease])  %>% 
    # Subset for original articles
    subset_research() %>% 
    extract_mir_df(threshold = 1,
                   extract_letters = TRUE)
  print(glue("Abstracts: For {diseases[disease]}  WITH subset RESEARCH {length(unique(subset_df$PMID))} unique abstracts were found."))
  print(glue("miRNAs:: For {diseases[disease]} WITH subset RESEARCH {length(unique(subset_df$miRNA))} unique miRNAs were found."))
  
  print("-------------- ")
  print("-------------- ")
}
