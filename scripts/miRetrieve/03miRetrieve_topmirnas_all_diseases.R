

require(readxl, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(janitor, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(miRetrieve, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(glue, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(gt, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(dplyr, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(tidyr, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(stringr, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(purrr, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(dplyr, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(ggplot2, lib.loc = "/mnt/users/reich/programs/R43/lib")


diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>%   # arranged automatically
  filter(analysis=="selected")

miRetrieve_alldiseases <- tibble(PMID = NA, Topic=NA, TargetName=NA, Accession=NA)

for(i in 1:length(all_combis$diseases)){
  # LOOP start ----------------------------------------------------------------
  # files created in "harmonise-miRetrieve-v21.R"
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
  
  human_disease_biomarker_top50 <- human_disease_biomarker %>% 
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
  
  # TABLE miRNAs ::: get all PMIDs with BM score above threshold (threshold = 5 according to '01miRetrieve_loop.R') corresponding to unique miRNAs  -----------------------------------------------
  no.unique.articles <- length(unique(human_disease_biomarker$PMID))
  
  mirna_pmids <- human_disease_biomarker %>% 
    group_by(miRNA) %>% 
    summarize(PMIDs_aboveBMscore5 = toString(unique(PMID)))
  
  table_mirna <- human_disease_biomarker_top50 %>% 
    left_join(mirna_pmids, by=c("miRNA"="miRNA")) %>% 
    select(Topic, TargetName, miRNA, Accession, Biomarker_score, miRetrieve, PubMed, PMIDs_aboveBMscore5)
  
  # gpt response load
  gpt_response <- as_tibble(read.csv(file = glue("./output/gpt_dataframe/{toupper(all_combis$diseases[i])}/2023-06-30-df-pubmed-gpt_response.csv"))) %>% 
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
              prognostic = round(mean(mortality == "Yes", na.rm=TRUE), 2)) -> test
  table_mirna %>% 
    left_join(gpt_response_structured, by = c("miRNA"="miRNA")) -> test
  View(test)
  
  # prepare for matching to BestAgeing dat -------------------------------------
  human_disease_biomarker_top50$TargetName <-  make_clean_names(human_disease_biomarker_top50$TargetName) %>% 
    str_replace(pattern = "mi_r", replacement = "mir")
  # sanity check if miRNA names are in sequenced df
  human_disease_biomarker_selection <- human_disease_biomarker_top50 %>% 
    filter(TargetName %in% names(data01)) %>% 
    select(PMID, Topic, TargetName, Accession)
  
  miRetrieve_alldiseases <- rbind(miRetrieve_alldiseases, human_disease_biomarker_selection)
  
}

miRetrieve_alldiseases <- miRetrieve_alldiseases[-1, ]
saveRDS(object = miRetrieve_alldiseases, file = "/mnt/users/reich/rockerprojects/bestageing2022/data-literature/miRetrieve/top50mirnas_all_diseases.rds")

length(unique(miRetrieve_alldiseases$Accession))


# keywords biomarker from miRetrieve
miRetrieve::biomarker_keywords

