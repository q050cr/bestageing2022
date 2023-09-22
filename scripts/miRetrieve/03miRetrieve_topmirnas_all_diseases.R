

# Get system name
system_name <- Sys.info()["nodename"]
mount_filesystem <- TRUE
SAVE.files <- TRUE

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

require(readxl, lib.loc = lib_path)
require(janitor, lib.loc = lib_path)
require(miRetrieve, lib.loc = lib_path)
require(glue, lib.loc = lib_path)
require(dplyr, lib.loc = lib_path)
require(tidyr, lib.loc = lib_path)
require(stringr, lib.loc = lib_path)
require(purrr, lib.loc = lib_path)
require(dplyr, lib.loc = lib_path)
require(ggplot2, lib.loc = lib_path)


# MIRNA DAT
model_data1 <- clean_names(readRDS(file = glue('{data_path_BestAgeing}/data_new/model_data1.RDS')))  # has also multiclass col + diagnoses
load(file = glue('{data_path_BestAgeing}/data/mirnas.rda'))  # "UKL-HD" n=765
load(file = glue('{data_path_BestAgeing}/data/data.rda'))  # "UKL-HD" n=731
all_mirnas <- clean_names(mirnas)
names(all_mirnas) <- str_replace(string = names(all_mirnas), pattern = "mi_r", replacement = "mir")
mirnas_disease <- clean_names(data)
names(mirnas_disease) <- str_replace(string = names(mirnas_disease), pattern = "mi_r", replacement = "mir")
rm(data, mirnas)


diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>%   # arranged automatically
  filter(analysis=="selected")

miRetrieve_alldiseases <- tibble(PMID = NA, Topic=NA, TargetName=NA, Accession=NA)
table_mirna_top50bm_score_alldiseases <- tibble(
  Topic = NA,
  TargetName = NA,
  miRNA = NA,
  Accession = NA,
  Biomarker_score = NA,
  miRetrieve = NA,
  max_value = NA,
  PMIDs = NA,
  upregulated_percent_study = NA,
  serum_plasma_tissue = NA,
  sample_size_studies = NA,
  prognostic = NA
)

miRetrieve_alldiseases_ALL <- tibble(PMID = NA, Topic=NA, TargetName=NA, Accession=NA)
table_mirna_ALL_bm_score_alldiseases <- tibble(
  Topic = NA,
  TargetName = NA,
  miRNA = NA,
  Accession = NA,
  Biomarker_score = NA,
  miRetrieve = NA,
  max_value = NA,
  PMIDs = NA,
  upregulated_percent_study = NA,
  serum_plasma_tissue = NA,
  sample_size_studies = NA,
  prognostic = NA
)

for(i in 1:length(all_combis$diseases)){
  # LOOP start ----------------------------------------------------------------
  # files created in "harmonise-miRetrieve-v21.R"
  path2biomarker <- glue("{data_path_bestageing2022}/data-literature/miRetrieve/{all_combis$diseases[i]}/2023-07-27-human-disease_biomarker_with_accession.rds")  # changed from 2023-07-12
  path2biomarker_count <- glue("{data_path_bestageing2022}/data-literature/miRetrieve/{all_combis$diseases[i]}/2023-07-27-human-df_count_both_with_accession.rds")
  path2biomarker_with_article_count <-  glue("{data_path_bestageing2022}/data-literature/miRetrieve/{all_combis$diseases[i]}/2023-07-27-human-disease_biomarker_with_article_count.rds")
  human_disease_biomarker <- readRDS(file = path2biomarker)
  human_disease_biomarker_count <- readRDS(file = path2biomarker_count)
  human_disease_biomarker_with_article_count <- readRDS(file = path2biomarker_with_article_count)
  
  # 1) most hits on miRetrieve
  human_disease_biomarker_count_tidy <- 
    human_disease_biomarker_count %>%
    drop_na(Accession) %>%
    arrange(desc(miRetrieve)) %>%
    dplyr::slice(1:30)
  
  
  # NO slicing top 50 -------------------------------------------------------
  
  human_disease_biomarker_ALL <- human_disease_biomarker %>% 
    drop_na(Accession) %>% 
    arrange(desc(Biomarker_score)) %>% 
    left_join(human_disease_biomarker_count %>% 
                select(Accession, miRetrieve, PubMed) %>% 
                drop_na(Accession) %>% 
                filter(!duplicated(Accession)), 
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
    drop_na(TargetName) #%>% 
    #slice(1:50)  # top 50 
  
  # change here
  human_disease_biomarker_ALL <- human_disease_biomarker_ALL %>% 
    select(-c(miRetrieve, PubMed, Biomarker_score, max_value)) %>% 
    left_join(human_disease_biomarker_with_article_count %>% mutate(miRNA=paste0("hsa-", miRNA)),
              by=c("miRNA"="miRNA", "Topic"="Topic"))
  
  # TABLE miRNAs ::: get all PMIDs with BM score above threshold (threshold = 5 according to '01miRetrieve_loop.R') corresponding to unique miRNAs  -----------------------------------------------
  no.unique.articles <- length(unique(human_disease_biomarker$PMID))
  
  mirna_pmids <- human_disease_biomarker %>% 
    group_by(miRNA) %>% 
    summarize(PMIDs = toString(unique(PMID)))
  
  table_mirna_ALL <- human_disease_biomarker_ALL %>% 
    left_join(mirna_pmids, by=c("miRNA"="miRNA")) %>% 
    select(Topic, TargetName, miRNA, Accession, Biomarker_score, miRetrieve, max_value, PMIDs)
  
  # gpt response load
  gpt_response <- as_tibble(read.csv(file = glue("{data_path_bestageing2022}/output/gpt_dataframe/{toupper(all_combis$diseases[i])}/2023-06-30-df-pubmed-gpt_response.csv"))) %>% 
    mutate(
      direction_upreg_downreg = ifelse(direction_upreg_downreg == "Not Given", NA, direction_upreg_downreg),
      serum_plasma_tissue = ifelse(serum_plasma_tissue == "Not Sure", NA, serum_plasma_tissue),
      mortality = ifelse(mortality == "Not Sure", NA, mortality),
      sample_size = ifelse(sample_size == 0, NA, sample_size)) %>% 
    mutate(miRNA = paste0("hsa-", miRNA))
  
  gpt_response_structured <- gpt_response %>% 
    group_by(miRNA) %>% 
    summarize(upregulated_percent_study = round(mean(direction_upreg_downreg == "Upregulated", na.rm=TRUE), 2),
              serum_plasma_tissue = round(mean(serum_plasma_tissue == "Serum" | serum_plasma_tissue == "Plasma", na.rm=TRUE), 2),
              sample_size_studies = toString(sample_size), 
              prognostic = round(mean(mortality == "Yes", na.rm=TRUE), 2))
  # combine
  table_mirna_ALL <- table_mirna_ALL %>% 
    left_join(gpt_response_structured, by = c("miRNA"="miRNA")) %>% 
    mutate(
      TargetName = make_clean_names(TargetName) %>% 
        str_replace(pattern = "mi_r", replacement = "mir")
    ) %>% 
    filter(TargetName %in% names(mirnas_disease)) %>% 
    arrange(desc(miRetrieve))
  
  # prepare for matching to BestAgeing dat -------------------------------------
  human_disease_biomarker_ALL$TargetName <-  make_clean_names(human_disease_biomarker_ALL$TargetName) %>% 
    str_replace(pattern = "mi_r", replacement = "mir")
  # sanity check if miRNA names are in sequenced df
  human_disease_biomarker_selection_ALL <- human_disease_biomarker_ALL %>% 
    filter(TargetName %in% names(mirnas_disease)) %>% 
    select(PMID, Topic, TargetName, Accession)
  
  miRetrieve_alldiseases_ALL <- rbind(miRetrieve_alldiseases_ALL, human_disease_biomarker_selection_ALL)  # old
  table_mirna_ALL_bm_score_alldiseases <- rbind(table_mirna_ALL_bm_score_alldiseases, table_mirna_ALL)  # new
  

  # TOP 50 per disease ------------------------------------------------------

  human_disease_biomarker_top50 <- human_disease_biomarker %>% 
    drop_na(Accession) %>% 
    arrange(desc(Biomarker_score)) %>% 
    left_join(human_disease_biomarker_count %>% 
                select(Accession, miRetrieve, PubMed) %>% 
                drop_na(Accession) %>% 
                filter(!duplicated(Accession)), 
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
  
  # change here
  human_disease_biomarker_top50 <- human_disease_biomarker_top50 %>% 
    select(-c(miRetrieve, PubMed, Biomarker_score, max_value)) %>% 
    left_join(human_disease_biomarker_with_article_count %>% mutate(miRNA=paste0("hsa-", miRNA)),
              by=c("miRNA"="miRNA", "Topic"="Topic"))
  
  # TABLE miRNAs ::: get all PMIDs with BM score above threshold (threshold = 5 according to '01miRetrieve_loop.R') corresponding to unique miRNAs  -----------------------------------------------
  no.unique.articles <- length(unique(human_disease_biomarker$PMID))
  
  mirna_pmids <- human_disease_biomarker %>% 
    group_by(miRNA) %>% 
    summarize(PMIDs = toString(unique(PMID)))
  
  table_mirna_top50 <- human_disease_biomarker_top50 %>% 
    left_join(mirna_pmids, by=c("miRNA"="miRNA")) %>% 
    select(Topic, TargetName, miRNA, Accession, Biomarker_score, miRetrieve, max_value, PMIDs)
  
  # gpt response load
  gpt_response <- as_tibble(read.csv(file = glue("{data_path_bestageing2022}/output/gpt_dataframe/{toupper(all_combis$diseases[i])}/2023-06-30-df-pubmed-gpt_response.csv"))) %>% 
    mutate(
      direction_upreg_downreg = ifelse(direction_upreg_downreg == "Not Given", NA, direction_upreg_downreg),
      serum_plasma_tissue = ifelse(serum_plasma_tissue == "Not Sure", NA, serum_plasma_tissue),
      mortality = ifelse(mortality == "Not Sure", NA, mortality),
      sample_size = ifelse(sample_size == 0, NA, sample_size)) %>% 
    mutate(miRNA = paste0("hsa-", miRNA))
  
  gpt_response_structured <- gpt_response %>% 
    group_by(miRNA) %>% 
    summarize(upregulated_percent_study = round(mean(direction_upreg_downreg == "Upregulated", na.rm=TRUE), 2),
              serum_plasma_tissue = round(mean(serum_plasma_tissue == "Serum" | serum_plasma_tissue == "Plasma", na.rm=TRUE), 2),
              sample_size_studies = toString(sample_size), 
              prognostic = round(mean(mortality == "Yes", na.rm=TRUE), 2))
  # combine
  table_mirna_top50 <- table_mirna_top50 %>% 
    left_join(gpt_response_structured, by = c("miRNA"="miRNA")) %>% 
    mutate(
      TargetName = make_clean_names(TargetName) %>% 
        str_replace(pattern = "mi_r", replacement = "mir")
    ) %>% 
    filter(TargetName %in% names(mirnas_disease)) %>% 
    arrange(desc(miRetrieve))
  
  # prepare for matching to BestAgeing dat -------------------------------------
  human_disease_biomarker_top50$TargetName <-  make_clean_names(human_disease_biomarker_top50$TargetName) %>% 
    str_replace(pattern = "mi_r", replacement = "mir")
  # sanity check if miRNA names are in sequenced df
  human_disease_biomarker_selection <- human_disease_biomarker_top50 %>% 
    filter(TargetName %in% names(mirnas_disease)) %>% 
    select(PMID, Topic, TargetName, Accession)
  
  miRetrieve_alldiseases <- rbind(miRetrieve_alldiseases, human_disease_biomarker_selection)  # old
  table_mirna_top50bm_score_alldiseases <- rbind(table_mirna_top50bm_score_alldiseases, table_mirna_top50)  # new
  

}





# top 50 ------------------------------------------------------------------

miRetrieve_alldiseases <- miRetrieve_alldiseases[-1, ]
table_mirna_top50bm_score_alldiseases <- table_mirna_top50bm_score_alldiseases[-1, ]
# save
if (SAVE.files ==TRUE) {
  saveRDS(object = miRetrieve_alldiseases, file = glue("{data_path_bestageing2022}/data-literature/miRetrieve/{Sys.Date()}_top50mirnas_all_diseases.rds"))  # old
  saveRDS(object = table_mirna_top50bm_score_alldiseases, file = glue("{data_path_bestageing2022}/data-literature/miRetrieve/{Sys.Date()}_top50mirnas_all_diseases_pmids_gpt.rds"))  # new
}

length(unique(miRetrieve_alldiseases$Accession))


# all miRetrieve above bm threshold ---------------------------------------

miRetrieve_alldiseases_ALL <- miRetrieve_alldiseases_ALL[-1, ]
table_mirna_ALL_bm_score_alldiseases <- table_mirna_ALL_bm_score_alldiseases[-1, ]
# save
if (SAVE.files ==TRUE) {
  saveRDS(object = miRetrieve_alldiseases_ALL, file = glue("{data_path_bestageing2022}/data-literature/miRetrieve/{Sys.Date()}_ALL_bm_mirnas_all_diseases.rds"))  # old
  saveRDS(object = table_mirna_ALL_bm_score_alldiseases, file = glue("{data_path_bestageing2022}/data-literature/miRetrieve/{Sys.Date()}_ALL_bm_mirnas_all_diseases_pmids_gpt.rds"))  # new
}


# keywords biomarker from miRetrieve
miRetrieve::biomarker_keywords

