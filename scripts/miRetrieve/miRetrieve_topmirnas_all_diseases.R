


diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>%   # arranged automatically
  filter(analysis=="selected")

miRetrieve_alldiseases <- tibble(PMID = NA, Topic=NA, TargetName=NA, Accession=NA)

for(i in 1:length(all_combis$diseases)){
  path2biomarker <- glue("/mnt/users/reich/rockerprojects/bestageing2022/data-literature/miRetrieve/{all_combis$diseases[i]}/2023-07-12-human-disease_biomarker_with_accession.rds")
  path2biomarker_count <- glue("/mnt/users/reich/rockerprojects/bestageing2022/data-literature/miRetrieve/{all_combis$diseases[i]}/2023-07-12-human-df_count_both_with_accession.rds")
  human_disease_biomarker <- readRDS(file = path2biomarker)
  human_disease_biomarker_count <- readRDS(file = path2biomarker_count)
  
  # 1) most hits on miRetrieve
  human_disease_biomarker_count_tidy <- 
    human_disease_biomarker_count %>%
    drop_na(Accession) %>%
    arrange(desc(miRetrieve)) %>%
    dplyr::slice(1:30)
  
  human_disease_biomarker <- human_disease_biomarker %>% 
    drop_na(Accession) %>% 
    arrange(desc(Biomarker_score)) %>% 
    left_join(human_disease_biomarker_count %>% 
                select(Accession, miRetrieve, PubMed), 
              by=c("Accession"="Accession")) %>% 
    # drop_na(miRetrieve) %>%   # some PMIDs not in most hits ;)
    mutate(miRetrieve = ifelse(is.na(miRetrieve), 0, miRetrieve) ) %>% 
    ## !! GROUP BY distinct miRNAs!!
    group_by(Accession) %>% 
    mutate(max_value = (Biomarker_score + miRetrieve)) %>% 
    arrange(desc(max_value)) %>% 
    slice(1) %>% 
    ungroup() %>% 
    arrange(desc(max_value)) %>% 
    drop_na(TargetName) %>% 
    slice(1:50)  # top 50 
  
  human_disease_biomarker$TargetName <-  make_clean_names(human_disease_biomarker$TargetName) %>% 
    str_replace(pattern = "mi_r", replacement = "mir")
  # sanity check if miRNA names are in sequenced df
  human_disease_biomarker_selection <- human_disease_biomarker %>% 
    filter(TargetName %in% names(data01)) %>% 
    select(PMID, Topic, TargetName, Accession)
  
  miRetrieve_alldiseases <- rbind(miRetrieve_alldiseases, human_disease_biomarker_selection)
  
}

miRetrieve_alldiseases <- miRetrieve_alldiseases[-1, ]
saveRDS(object = miRetrieve_alldiseases, file = "/mnt/users/reich/rockerprojects/bestageing2022/data-literature/miRetrieve/top50mirnas_all_diseases.rds")

length(unique(miRetrieve_alldiseases$Accession))



