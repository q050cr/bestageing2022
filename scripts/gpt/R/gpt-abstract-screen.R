
# the query was made on personal Mac M1 with API key stored in environment vars

# CONSIDER
## Default rate limits  https://platform.openai.com/docs/guides/rate-limits
## https://github.com/openai/openai-cookbook/blob/main/examples/How_to_handle_rate_limits.ipynb
## 3,500 RPM (requests per minute)
## 90,000 TPM (tokens per minute)


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

#  Load libraries
require(miRetrieve, lib.loc = lib_path); packageVersion('miRetrieve')
require(glue, lib.loc = lib_path)
require(magrittr, lib.loc = lib_path) # Load magrittr for %>%
require(ggplot2, lib.loc = lib_path) # Load ggplot2 for plotting
require(svglite, lib.loc = lib_path)
require(easyPubMed, lib.loc = lib_path) # Load easyPubMed, access to PubMed
require(rcrossref, lib.loc = lib_path)
require(dplyr, lib.loc = lib_path) # Data wrangling
require(tidyr, lib.loc = lib_path) # Data wrangling
require(patchwork, lib.loc = lib_path) # Group graphs
#easyPubMed does not work
require(rentrez, lib.loc = lib_path)
require(pubmedR, lib.loc = lib_path)



library(openai)
library(chatgpt)
library(httr)
library(jsonlite)
library(reticulate)  # for 'tiktoken' library to count tokens
library(tidyverse)
library(glue)
## comment out if this script is used in the future!! 
# use_python("/Library/Frameworks/Python.framework/Versions/3.11/bin/python3", required = T)
# tiktoken <- import("tiktoken")
# encoding <- tiktoken$encoding_for_model("gpt-3.5-turbo")  # model
# myapikey <- Sys.getenv("OPENAI_API_KEY")
# gptmodel <- openai::list_models()$data$root[39]   # gpt-3.5-turbo-0613

source(file = "./scripts/gpt/gpt-fns.R")


# load data -----------------------------------------------------

diseases <- c("ACS", "CAD", "DCM", "HFrEF")
diseases_tibble <- tibble(abbrev_disease = diseases,
                          disease_name = c("Acute Coronary Syndrome", "Coronary Artery Disease", "Dilatative Cardiomyopathy", "Heart Failure"))
path2dat <- paste0(data_path_bestageing2022, "/data-literature/pubmed/", list.files(path = glue("{data_path_bestageing2022}/data-literature/pubmed/")))
path2dat <- path2dat[grep("-2023-06-16-", path2dat)]   # use newest files from 2023-06-16

cost_input <- numeric()
processed_input_tokens <- numeric()
for(disease in seq_along(diseases)) {
  
  
  my_abstracts <- as_tibble(read.csv2(file = path2dat[disease]))
  print(glue("For {diseases[disease]} {length(unique(my_abstracts$PMID))} unique abstracts were found."))
  print(glue("For {diseases[disease]} {length(unique(my_abstracts$miRNA))} unique miRNAs were found."))
  # check api costs
  #calculate_input_token_cost() # now done in helper-fn in "gpt-fns.R"
  
}

length(unique(my_abstracts$PMID))
length(unique(my_abstracts$miRNA))

grouped_abstracts <- my_abstracts %>% 
  group_by(PMID) %>% 
  summarize(miRNA = list(unique(miRNA)))

# initialize column where we save our query tibbles later
grouped_abstracts$df_temp_abstract_x <- vector("list", nrow(grouped_abstracts))
grouped_abstracts$answer <-  character(length = nrow(grouped_abstracts))
grouped_abstracts$usage_total_tokens <- integer(length = nrow(grouped_abstracts))
grouped_abstracts$query_cost <- numeric(length = nrow(grouped_abstracts))
grouped_abstracts$time <- numeric(length = nrow(grouped_abstracts))


# get disease name for query
abstract_disease <- diseases_tibble$disease_name[disease]
for (unique_abstract in 1:nrow(grouped_abstracts)) {
  # get abstract: CAVE dim(my_abstracts)[1] != dim(grouped_abstracts)[1], that's why we are using the first index with "min"
  selected_abstract_index <- min(which(my_abstracts$PMID %in% grouped_abstracts$PMID[unique_abstract]))
  selected_abstract <- my_abstracts$Abstract[selected_abstract_index]
  
  # get mirna set of interest
  abstract_mirna_list <- unlist(grouped_abstracts$miRNA[unique_abstract])
  abstract_mirna_list_prompt <- paste0(abstract_mirna_list, collapse = ", ")
  system_prompt <- glue('[no prose, no explanation, no notes] [output only json] You are a Data Scientist analyzing microRNA data. You read medical abstracts and return answers 
                        with the column names "related_topic" for Q1, "direction_upreg_downreg" for Q2, "primary_literature" for Q3, 
                        "serum_plasma_tissue" for Q4, "mortality" for Q5, "measurement_type" for Q6, and "sample_size" for Q7.
                        
                        You provide {length(abstract_mirna_list)} separate json structured data for each of the following {length(abstract_mirna_list)} microRNAs in this ordered list: MY_LIST=[{abstract_mirna_list_prompt}]
                        
                        Q1: Does the abstract deal with {abstract_disease} in humans? (Yes/No)
                        Q2: Is the given miRNA from MY_LIST up- (increased) or downregulated (decreased) in {abstract_disease}? (Upregulated/Downregulated)
                        Q3: Is the article an original research and not a secondary literature (review articles, systematic reviews, meta-analysis)? (Yes/No)
                        Q4: Was the miRNA from MY_LIST measured in serum, plasma, heart tissue or somewhere else? (Serum/Plasma/Tissue/Else)
                        Q5: Is miRNA from MY_LIST in the abstract clearly associated with mortality/death? (Yes, No, Not given)
                        Q6: How was  miRNA from MY_LIST measured?
                        Q7: What was the human sample size in the study? (number)')
  user_prompt <- glue('Investigate the following Abstract: "{selected_abstract}"')
  
  # call to API
  resultsFromQuery <- retry_chatGPT(system_prompt = system_prompt, user_prompt = user_prompt, retries = 3)
  
  # resultsFromQuery <- chatGPT(system_prompt = system_prompt, user_prompt = user_prompt, cost_per_1000_input_tokens = 0.0015,
  #                             modelName = "gpt-3.5-turbo", temperature = 0.7)
  ## store results
  #grouped_abstracts$df_temp_abstract_x[[unique_abstract]] <- resultsFromQuery$df_temp_abstract_x
  grouped_abstracts$usage_total_tokens[unique_abstract] <- resultsFromQuery$usage_total_tokens
  grouped_abstracts$query_cost[unique_abstract] <- resultsFromQuery$query_cost
  grouped_abstracts$time[unique_abstract] <- resultsFromQuery$time
  grouped_abstracts$answer[unique_abstract] <- resultsFromQuery$answer
  
  print(glue::glue("|||-----------------------Run finished for unique abstract: {unique_abstract}/{nrow(grouped_abstracts)}-----------------------|||"))
}
sum_cost <- glue("{round(sum(grouped_abstracts$query_cost),3)} $")

saveRDS(object = grouped_abstracts, file = glue("./data-literature/pubmed-gpt/raw_results/{diseases[disease]}-{Sys.Date()}-human-gptpubmed_query-unique-articles{nrow(grouped_abstracts)}.rds"))



# saveRDS(object = grouped_abstracts, file = glue("./data-literature/pubmed-gpt/query_results/{diseases[disease]}-{Sys.Date()}-human-gptpubmed_query-unique-articles{nrow(grouped_abstracts)}.rds"))
