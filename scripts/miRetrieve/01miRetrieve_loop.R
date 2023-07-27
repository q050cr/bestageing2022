

## 
# file:///Users/christophreich/Downloads/lqab117_supplemental_files/Supplementary%20File%205%20Second%20Revision.html

#version$version.string

#  Load libraries
require(miRetrieve, lib.loc = "/mnt/users/reich/programs/R43/lib"); packageVersion('miRetrieve')
require(glue, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(magrittr, lib.loc = "/mnt/users/reich/programs/R43/lib") # Load magrittr for %>%
require(ggplot2, lib.loc = "/mnt/users/reich/programs/R43/lib") # Load ggplot2 for plotting
require(svglite, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(easyPubMed, lib.loc = "/mnt/users/reich/programs/R43/lib") # Load easyPubMed, access to PubMed
require(rcrossref, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(dplyr, lib.loc = "/mnt/users/reich/programs/R43/lib") # Data wrangling
require(tidyr, lib.loc = "/mnt/users/reich/programs/R43/lib") # Data wrangling
require(patchwork, lib.loc = "/mnt/users/reich/programs/R43/lib") # Group graphs
#easyPubMed does not work
require(rentrez, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(pubmedR, lib.loc = "/mnt/users/reich/programs/R43/lib")

source(file = "/mnt/users/reich/rockerprojects/bestageing2022/scripts/helper/custom_ggplot_theme.R")

saveFILE <- TRUE

##
# QUERIES ---------
##
diseases <- c("ACS", "CAD", "DCM", "HFrEF")

query_list <- list()
# initial query until 2022
query_list[["ACS"]] <- '(("MicroRNAs"[Mesh] OR "miRNAs"[Title/Abstract] OR "microRNA"[Title/Abstract] OR "miRNA"[Title/Abstract]) 
AND ("Acute Coronary Syndrome"[Mesh] OR "Acute Coronary Syndrome" OR "Myocardial Infarction"[Mesh])) 
AND 2000:2022[DP]'

# updated queries 2023, added diagnosis and biomarker terms in June 2023
query_list[["ACS"]] <- '(("MicroRNAs"[Mesh] OR "miRNAs"[Title/Abstract] OR "micro RNA"[Title/Abstract]) 
AND ("Acute Coronary Syndrome"[Mesh] OR "Acute Coronary Syndrome" OR "Myocardial Infarction"[Mesh]) 
AND ("Biomarkers"[Mesh] OR "Biomarkers"[Title/Abstract])
AND ("Diagnosis"[Mesh] OR "Diagnostic"[Title/Abstract] OR "Diagnostic Techniques and Procedures"[Mesh] OR "Diagnosis"[Title/Abstract])
AND "Humans"[Mesh] AND 2000:2023[DP] 
NOT ("Animals"[Mesh] NOT "Humans"[Mesh] OR "Cell Line"[Mesh] OR "Animal Experimentation"[Mesh] 
OR "In Vitro Techniques"[Mesh] OR "Cell Culture Techniques"[Mesh] OR "animal model"[Title/Abstract] OR "laboratory study"[Title/Abstract] OR "cell line"[Title/Abstract]))'

query_list[["CAD"]] <- '(("MicroRNAs"[Mesh] OR "miRNAs"[Title/Abstract] OR "microRNA"[Title/Abstract] OR "miRNA"[Title/Abstract]) 
AND ("Coronary Artery Disease"[Mesh] OR "Coronary Artery Disease" )) 
AND ("Biomarkers"[Mesh] OR "Biomarkers"[Title/Abstract])
AND ("Diagnosis"[Mesh] OR "Diagnostic"[Title/Abstract] OR "Diagnostic Techniques and Procedures"[Mesh] OR "Diagnosis"[Title/Abstract])
AND "Humans"[Mesh] AND 2000:2023[DP] 
NOT ("Animals"[Mesh] NOT "Humans"[Mesh] OR "Cell Line"[Mesh] OR "Animal Experimentation"[Mesh] 
OR "In Vitro Techniques"[Mesh] OR "Cell Culture Techniques"[Mesh] OR "animal model"[Title/Abstract] OR "laboratory study"[Title/Abstract] OR "cell line"[Title/Abstract]))'

query_list[["DCM"]] <- '(("MicroRNAs"[Mesh] OR "miRNAs"[Title/Abstract] OR "microRNA"[Title/Abstract] OR "miRNA"[Title/Abstract]) 
AND ("Cardiomyopathy, Dilated"[Mesh] OR "Cardiomyopathy, Dilated" OR "dilated cardiomyopathy")) 
AND ("Biomarkers"[Mesh] OR "Biomarkers"[Title/Abstract])
AND ("Diagnosis"[Mesh] OR "Diagnostic"[Title/Abstract] OR "Diagnostic Techniques and Procedures"[Mesh] OR "Diagnosis"[Title/Abstract])
AND "Humans"[Mesh] AND 2000:2023[DP] 
NOT ("Animals"[Mesh] NOT "Humans"[Mesh] OR "Cell Line"[Mesh] OR "Animal Experimentation"[Mesh] 
OR "In Vitro Techniques"[Mesh] OR "Cell Culture Techniques"[Mesh] OR "animal model"[Title/Abstract] OR "laboratory study"[Title/Abstract] OR "cell line"[Title/Abstract]))'

query_list[["HFrEF"]] <- '(("MicroRNAs"[Mesh] OR "miRNAs"[Title/Abstract] OR "microRNA"[Title/Abstract] OR "miRNA"[Title/Abstract]) 
AND ("Heart Failure"[Mesh] OR "Heart Failure, Systolic"[Mesh] OR "Heart Failure, Systolic" OR "heart failure with reduced ejection fraction"))
AND ("Biomarkers"[Mesh] OR "Biomarkers"[Title/Abstract])
AND ("Diagnosis"[Mesh] OR "Diagnostic"[Title/Abstract] OR "Diagnostic Techniques and Procedures"[Mesh] OR "Diagnosis"[Title/Abstract])
AND "Humans"[Mesh] AND 2000:2023[DP] 
NOT ("Animals"[Mesh] NOT "Humans"[Mesh] OR "Cell Line"[Mesh] OR "Animal Experimentation"[Mesh] 
OR "In Vitro Techniques"[Mesh] OR "Cell Culture Techniques"[Mesh] OR "animal model"[Title/Abstract] OR "laboratory study"[Title/Abstract] OR "cell line"[Title/Abstract]))'

# adapted queries for only looking at clinical trials, meta-analysis and randomized controlled trials
add_string <- 'AND ("Clinical Trial"[Publication Type] OR "Meta-Analysis"[Publication Type] OR "Randomized Controlled Trial"[Publication Type])'
query_list_trials <- list()
query_list_trials[["ACS"]] <- glue::glue({query_list[["ACS"]]},{add_string})
query_list_trials[["CAD"]] <- glue::glue({query_list[["CAD"]]},{add_string})
query_list_trials[["DCM"]] <- glue::glue({query_list[["DCM"]]},{add_string})
query_list_trials[["HFrEF"]] <- glue::glue({query_list[["HFrEF"]]},{add_string})


# fn to obtain PubMed IDs from query "{miRNA} + {disease}" 
compare_count <- function(mir, mesh_query) {
  query <- paste0('(',mir, ') AND ', mesh_query)
  pubmed_ids <- entrez_search(db="pubmed", term=query, retmax=9999)$ids # retmax max 9999, only EDirect (cmd line tool can obtain more queries, but we could e.g. split time)
  return(pubmed_ids)
}

# fn figure 1

count_plot_fn <- function(df_count) {
  count_plot <- df_count %>%
    mutate(miRNA = forcats::fct_reorder(miRNA, miRetrieve)) %>% 
    pivot_longer(cols = c(miRetrieve, PubMed)) %>%
    mutate(name = forcats::fct_rev(name)) %>% 
    ggplot(aes(x = miRNA, y = value, fill = name)) +
    geom_col(position = "dodge") +
    coord_flip() +
    guides(fill = guide_legend(reverse = TRUE, title = "Method")) +
    # scale_fill_brewer(palette = "Dark2") +
    # theme_classic() +
    ylab("# of articles") + 
    scale_y_continuous(expand = c(0,0)) +
    theme_minimal(base_size = 16, base_family = 'Source Sans Pro')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    my_base_theme()
  
  return(count_plot)
}



##
# QUERY LOOP -----------------------------------------------------------------
##
docl <- list()  # create a list object for storing results

topn.mirna.plt <- 20 # how many mirnas should be plotted in barplot?

for (disease in seq_along(diseases)) {
  subset <- "-human" # changed to "-human" only, could also be changed to "-humantrials" only, or "" for full pubmed search !!!, only works without date in next line!
  file <- glue::glue("./data-literature/pubmed-downloads/2023-06-14-pubmed-MicroRNA{subset}-{toupper(diseases[disease])}.txt")  
  
  ## load abstracts ----------------------------------------------------------
  df <- read_pubmed(file, topic = diseases[disease]) %>% 
    # Subset for original articles
    subset_research() %>%
    # Extract miRNA names mentioned at least twice/abstract
    extract_mir_df(threshold = 2)
  
  # Number of abstracts
  print(glue::glue("{length(unique(df[['PMID']]))} unique abstracts with extracted miRNA names mentioned at least twice/ abstract for disease: {diseases[disease]}"))
  # Number of rows
  nrow(df)
  
  # receive dois
  # testid <- id_converter(x = df$PMID).  Get a PMID from a DOI, and vice versa.
  
  ## Compare extraction count to PubMed --------------------------------------
  # Extract all miRNAs
  df_comparison <- 
    read_pubmed(file,
                topic = diseases[disease]) %>% 
    extract_mir_df(threshold = 1,
                   extract_letters = TRUE)
  no.unique.articles <- length(unique(df_comparison$PMID))
  print(glue("---There are {no.unique.articles} unique articles concerning {diseases[disease]}---"))
  
  ## Pubmed Manual Inquiry and for GPT API ------------------------------------------------
  df_comparison_manual_pubmed <- df_comparison %>% 
    # Subset for original articles
    subset_research() %>%
    # filter for Journal Articles only (no "Comment"); implemented before the above function, should do the same...
    filter(purrr::map_lgl(Type, ~ any(. == "Journal Article")) ) %>% 
    mutate(unrelated_topic01=NA, direction_upreg_downreg=NA, english="yes", 
           related_topic=NA, direction_upreg_downreg=NA, primary_literature=NA, 
           serum_plasma_tissue=NA, mortality=NA, measurement_type=NA, sample_size=NA,
           comment=NA) %>% 
    arrange(PMID)
  
  df_comparison_manual_pubmed$Type <- df_comparison_manual_pubmed$Type %>% 
    purrr::map_chr(., ~ paste(., collapse = ", "))
  
  #if (saveFILE == TRUE) {
  #  write.csv2(x = df_comparison_manual_pubmed, 
  #             file = glue::glue("./data-literature/pubmed/{diseases[disease]}-{Sys.Date()}{subset}-pubmed_check-unique-articles{no.unique.articles}.csv"))
  #}
  
  # Top 30 miRNAS letters/ no letters --------------------------------------------------------
  # Count miRNAs Top 30 
  df_count_letter <- 
    df_comparison %>% 
    count_mir() %>% 
    # Filter for miRNAs with lettered suffixes
    dplyr::filter(stringr::str_detect(miRNA, "\\d[a-z]")) %>% 
    dplyr::slice(1:30)
  
  # Filter for miRNAs with lettered suffixes
  # Extract miRNA numbers of miRNAs with lettered suffixes
  df_mir_letter <- 
    df_comparison %>% 
    count_mir() %>% 
    dplyr::filter(stringr::str_detect(miRNA, "\\d[a-z]")) %>% 
    filter(Mentioned_n >= 10) %>%                               # filter mentioned >= 10 !!
    pull(miRNA) %>% 
    stringr::str_extract_all("miR-\\d+|let-7") %>% 
    unlist()
  
  # Filter for miRNAs without lettered suffixes Top 30
  df_count_no_letter <- 
    df_comparison %>% 
    count_mir() %>% 
    dplyr::filter(stringr::str_detect(miRNA, "\\d$")) %>% 
    filter(!miRNA %in% df_mir_letter) %>% 
    dplyr::slice(1:30)
  
  # Obtain PubMed IDs from query "{miRNA} + {disease}" ----------------------------
  # takes time for querying entrez
  count_pubmed_vec <- purrr::map(df_count_no_letter$miRNA, ~compare_count(mir = .x, mesh_query = query_list[[diseases[disease]]])) %>% 
    purrr::set_names(df_count_no_letter$miRNA)
  
  count_pubmed_vec_letter <- purrr::map(df_count_letter$miRNA, ~compare_count(mir = .x, mesh_query = query_list[[diseases[disease]]])) %>% 
    purrr::set_names(df_count_letter$miRNA)
  
  ## Figure 1: Compare Top 30 miRetrieve with PubMeds miRNA count -------------------
  
  ## a) no letter
  # Obtain number of PubMed results
  length_pmid <- purrr::map(count_pubmed_vec, length) %>% 
    unlist()
  concatenated_list_no_letter <- purrr::map(count_pubmed_vec, ~paste(.x, collapse = ", ")) %>% 
    unlist()
  # Add number of PubMed results to dataframe "df_count_no_letter"
  df_count_no_letter[3] <- length_pmid
  df_count_no_letter$PMIDs <- concatenated_list_no_letter
  names(df_count_no_letter)[3] <- "PubMed"
  names(df_count_no_letter)[2] <- "miRetrieve"
  # Plot # of miRNAs with miRetrieve vs. PubMed (no letters)
  count_plot1 <- count_plot_fn(df_count = df_count_no_letter)
  
  ## b) with letter
  # Obtain number of PubMed results
  length_pmid_letter <- purrr::map(count_pubmed_vec_letter, length) %>% 
    unlist()
  concatenated_list_letter <- purrr::map(count_pubmed_vec_letter, ~paste(.x, collapse = ", ")) %>% 
    unlist()
  # Add number of PubMed results to dataframe "df_count_letter"
  df_count_letter[3] <- length_pmid_letter
  df_count_letter$PMIDs <- concatenated_list_letter
  names(df_count_letter)[3] <- "PubMed"
  names(df_count_letter)[2] <- "miRetrieve" 
  
  # combine with and without letter suffix
  df_count_both <- rbind(df_count_letter, df_count_no_letter) %>%  
    arrange(desc(PubMed))
  
  # Plot # of miRNAs with miRetrieve vs. PubMed (no letters)
  count_plot1_letter <- count_plot_fn(df_count = df_count_letter)
  
  # Combine plots with {patchwork} (Figure 1)
  combined <- count_plot1 + count_plot1_letter & theme(legend.position = "bottom")
  
  combined <- combined + 
    plot_layout(guides = "collect") + 
    plot_annotation(tag_levels = "A")
  print(combined)
  
  # only plot top 30miRNAs mixed letters FINAL PLOT 1
  df_count_mixed <- bind_rows(df_count_letter, df_count_no_letter) %>% arrange(desc(miRetrieve)) %>% 
    slice(1:30) %>% 
    mutate(miRNA = paste0("hsa-", miRNA))
  count_plot_mixed <- count_plot_fn(df_count=df_count_mixed)
  print(count_plot_mixed)
  
  if (saveFILE == TRUE) {
    saveRDS(object = df_comparison,  # all miRNAs in df 
            file = glue::glue("./data-literature/miRetrieve/{tolower(diseases[disease])}/{Sys.Date()}{subset}-df_comparison.rds"))
    write.csv2(x = df_comparison %>% select(-Type),  # is a list column...
               file = glue::glue("./data-literature/miRetrieve/{tolower(diseases[disease])}/{Sys.Date()}{subset}-df_comparison.csv"))
    saveRDS(object = df_count_both,  # count top 30 mirnas with letter suffix and top 30 without 
            file = glue::glue("./data-literature/miRetrieve/{tolower(diseases[disease])}/{Sys.Date()}{subset}-df_count_both.rds"))
    write.csv2(x = df_count_both,  # count top 30 mirnas with letter suffix and top 30 without 
            file = glue::glue("./data-literature/miRetrieve/{tolower(diseases[disease])}/{Sys.Date()}{subset}-df_count_both.csv"))
    
    filename.fig1 <- glue::glue("./output/plots/miRetrieve/{diseases[disease]}/{Sys.Date()}{subset}-figure1-compareTop30")
    ggsave(filename = paste0(filename.fig1, ".svg"), plot =  combined, 
           width = 10, height = 6, 
           units = "in"  # default
    )
    filename.fig1_mixed <- glue::glue("./output/plots/miRetrieve/{diseases[disease]}/{Sys.Date()}{subset}-figure1-compareTop30_mixed_letter")
    ggsave(filename = paste0(filename.fig1_mixed, ".svg"), plot =  count_plot_mixed, 
           width = 8, height = 6, 
           units = "in"  # default
    )
  }
  
  # Calculate score [miRetrieve count] / [PubMed count] (no letters)
  df_count_no_letter %>% 
    mutate(Score = miRetrieve / PubMed) %>%
    mutate(Score = ifelse(is.infinite(Score), 0.01, Score)) %>%
    summarise(mean(Score), sd(Score), sum(miRetrieve), sum(PubMed))
  
  # Calculate score [miRetrieve count] / [PubMed count] (with letters)
  df_count_letter %>% 
    mutate(Score = miRetrieve / PubMed) %>%
    mutate(Score = ifelse(is.infinite(Score), 0.01, Score)) %>%
    summarise(mean(Score), sd(Score), sum(miRetrieve), sum(PubMed))
  
  ## Figure 2: Top miRNAs in {Disease} ---------------------------------------
  ## equals the `count_plot_mixed` plot!!
  #most_frequently_miRNA <- plot_mir_count(df_comparison,
  #                                        top = topn.mirna.plt,  # defined above how many miRNAs should be plotted
  #                                        title = glue::glue("Most frequently mentioned miRNAs in {diseases[disease]}")) +
  #  theme_minimal(base_size = 16, base_family = 'Source Sans Pro') +
  #  scale_fill_manual(values = thematic::okabe_ito(6)) +
  #  my_base_theme()
  
  #if (saveFILE == TRUE) {
  #  filename.fig2 <- glue::glue("./output/plots/miRetrieve/{diseases[disease]}/{Sys.Date()}{subset}-figure2-mostfrequent")
  #  ggsave(filename = paste0(filename.fig2, ".svg"), plot =  most_frequently_miRNA, 
  #         width = 8, height = 8, 
  #         units = "in"  # default
  #  )
  #}
  #
  ## Figure 7: Potential miRNA biomarker in diseases --------------------
  # pinpointed abstracts containing at least five marker 
  #   words for biomarker(bio-marker, biological marker, 
  #   biomarker, body fluid,bodyfluid, circulating, diagnostic, 
  #   exosomal, exosomes,extracellular vesicles, plasma, serum, urinary, urine). 
  
  # Potential miRNA biomarker in ACS
  disease_biomarker <- calculate_score_biomarker(df_comparison,
                                                 biomarker_keywords[-c(5,6,9,10)],  # no urine, body fluid
                                                 # circulating, biomarker, bio-marker, extracellular vesicles, exosomes, exosomal, diagnostic, biological marker, biomarker, bio-marker, biomarker, bio-marker, serum, plasma
                                                 threshold = 5,
                                                 discard = TRUE)
  
  # Plot top 7 miRNA biomarker in DCM
  biomarker_plot <- plot_mir_count(disease_biomarker, 
                                   top = topn.mirna.plt,
                                   colour = thematic::okabe_ito(6)[2],
                                   title = glue::glue("Potential biomarker miRNAs in {diseases[disease]}")) +
    theme_minimal(base_size = 16, base_family = 'Source Sans Pro') +
    labs(title = NULL) +
    #scale_fill_manual(values = thematic::okabe_ito(6)[3]) +
    #scale_colour_manual(values = thematic::okabe_ito(6)[3]) +
    #geom_col(alpha=0.5)+
    my_base_theme()
  print(biomarker_plot)
  # modify plot to really only include 20 bars (problem with DCM loads of miRNAs have only 1 abstract)
  biomarker_plot <- biomarker_plot$data %>% 
    arrange(desc(n)) %>% 
    mutate(miR=paste0("hsa-",miR)) %>% 
    mutate(miR = forcats::fct_reorder(miR,n)) %>% 
    slice(1:topn.mirna.plt) %>% 
    ggplot(aes(x=miR, y = n)) + 
    geom_col(fill = thematic::okabe_ito(6)[2]) +
    coord_flip() + xlab("miRNA") + ylab("Mentioned in # of abstracts") + 
    theme_minimal(base_size = 16, base_family = 'Source Sans Pro') +
    my_base_theme()
  
  if (saveFILE == TRUE) {
    # RDS will be used in "02harmonise-miRetrieve-v21.R" script
    saveRDS(object = disease_biomarker, 
            file = glue::glue("./data-literature/miRetrieve/{tolower(diseases[disease])}/{Sys.Date()}{subset}-disease_biomarker.rds"))
    #write.csv2(x = disease_biomarker, 
    #           file = glue::glue("./data-literature/miRetrieve/{tolower(diseases[disease])}/{Sys.Date()}{subset}-disease_biomarker.csv"))
    
    filename.fig7 <- glue::glue("./output/plots/miRetrieve/{diseases[disease]}/{Sys.Date()}{subset}-figure7-mostfrequent_potential_miRNA_biomarker")
    ggsave(filename = paste0(filename.fig7, ".svg"), plot =  biomarker_plot, 
           width = 6, height = 6, 
           units = "in"  # default
    )
  }
  print(glue::glue("|||-----------------------Run finished for disease: {diseases[disease]}-----------------------|||"))
}



# not included in loop ----------------------------------------------------

# Figure 3: Top terms {disease} -------------------------------------------
# Generate new stop word list specifically for diseases to 
# improve overview
disease_stop <- generate_stopwords(c("atherosclerosis",
                                     "atherosclerotic",
                                     "heart failure",
                                     "catheterization",
                                     "human",
                                     "blood",
                                     "tissue",
                                     "iPS", "iPSC",
                                     "troponin",
                                     "nt-proBNP",
                                     "biomarker",
                                     "patients",
                                     "pathogenesis",
                                     "umbilical",
                                     "huvecs",
                                     "umbilical vein",
                                     "human umbilical",
                                     "cardiovascular disease",
                                     "atherosclerosis however",
                                     "atherosclerosis as",
                                     "development atherosclerosis",
                                     "pathogenesis atherosclerosis",
                                     "oxidized low",
                                     "vein endothelial",
                                     "cardiovascular diseases",
                                     "cells huvecs",
                                     "necrosis factor",
                                     "cells vsmcs",
                                     "cells ecs",
                                     "demonstrated mir",
                                     "low density",
                                     "atherosclerotic plaque",
                                     "artery disease",
                                     "muscle cells",
                                     "tumor necrosis",
                                     "coronary artery",
                                     "dependent manner",
                                     "vascular endothelial",
                                     "patients coronary",
                                     "serum levels",
                                     "vascular smooth",
                                     "α"),
                                   combine_with = stopwords_miretrieve)

# Plot single word terms for miR-133
mir133_1 <- plot_mir_terms(df_acs, "miR-133", 
                           stopwords = disease_stop,
                           top = 15,
                           title = "Top single terms for miR-133 in ACS")

# Plot 2-grams for miR-133
mir133_2 <-  plot_mir_terms(df_acs, "miR-133", token = "ngrams", n = 2,
                            stopwords = disease_stop,
                            top = 14,
                            title = "Top 2-grams for miR-133 in ACS")

# Plot single word terms for miR-21
mir21_1 <-  plot_mir_terms(df_acs, "miR-21",
                           stopwords = disease_stop,
                           top = 11,
                           title = "Top single terms for miR-21 in ACS") 

# Plot 2-grams for miR-21
mir21_2 <- plot_mir_terms(df_acs, "miR-21", token = "ngrams", n = 2,
                          stopwords = disease_stop,
                          top = 13,
                          title = "Top 2-grams for miR-21 in acs")

# Combine plots (Figure 3)
(mir133_1 | mir133_2) /
  (mir21_1 | mir21_2) + plot_annotation(tag_levels = 'A')



# Figure 4: Targets in {disease} ------------------------------------------

# Add miRTarBase targets to df
df_targets_acs_mir <- join_mirtarbase(df_acs)

# Plot targets of miR-155 and miR-21 (Figure 4)
plot_target_mir_scatter(df_targets_acs_mir,
                        mir = c("miR-133", "miR-21"),
                        alpha = 1,
                        title = "Targets of miR-21 and miR-133 in ACS")



# Disease Associations ----------------------------------------------------

# Top terms in DCM
plot_mir_count(df_dcm,
               title = "Most frequently mentioned miRNAs in DCM")


# Figure 6: Compare miR associations between Diseases ---------------------

# Generate stop words for miR-21 comparison to improve overview
comparison_stop <- generate_stopwords(c("dependent manner",
                                        "vascular endothelial",
                                        "peripheral blood",
                                        "serum mir",
                                        "serum levels",
                                        "tensin homolog",
                                        "phosphatase tensin",
                                        "α"),
                                      combine_with = stopwords_miretrieve)

# Combine
df_acs_dcm <- combine_df(df_acs, df_dcm)

# Compare shared terms, single word
comp_mir21_1 <- compare_mir_terms(df_acs_dcm, "miR-21", top = 7,
                                  title = "Comparison of miR-21 single term association",
                                  stopwords = comparison_stop) +
  ggplot2::scale_fill_manual(values = c("#F5793A", "#85C0F9"), name = "Topic")

# Compare shared terms, 2-gram
comp_mir21_2 <- compare_mir_terms(df_acs_dcm, "miR-21", token = "ngrams", n = 2, top = 9,
                                  title = "Comparison of miR-21 2-gram association",
                                  stopwords = comparison_stop) +
  ggplot2::scale_fill_manual(values = c("#F5793A", "#85C0F9"), name = "Disease")

# Combine plot (Figure 6)
plot_combined <- (comp_mir21_1 | comp_mir21_2) & theme(legend.position = "bottom")

plot_combined + 
  plot_annotation(tag_levels = 'A')  + 
  plot_layout(guides = "collect")


