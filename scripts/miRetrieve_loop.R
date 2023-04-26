
## 
# file:///Users/christophreich/Downloads/lqab117_supplemental_files/Supplementary%20File%205%20Second%20Revision.html

version$version.string

#  Load libraries
library(miRetrieve); packageVersion('miRetrieve')
library(magrittr) # Load magrittr for %>%
library(ggplot2) # Load ggplot2 for plotting
library(easyPubMed) # Load easyPubMed, access to PubMed
library(dplyr) # Data wrangling
library(tidyr) # Data wrangling
library(patchwork) # Group graphs
#easyPubMed does not work
library(rentrez)
library(pubmedR)

saveFILE <- TRUE


diseases <- c("ACS", "CAD", "DCM", "HFrEF")

query_list <- list()
query_list[["ACS"]] <- '(("MicroRNAs"[Mesh] OR "miRNAs"[Mesh] OR "microRNA"[Title/Abstract] OR "miRNA"[Title/Abstract]) 
AND ("Acute Coronary Syndrome"[Mesh] OR "Acute Coronary Syndrome" OR "Myocardial Infarction"[Mesh])) 
AND 2000:2022[DP]'
query_list[["CAD"]] <- '(("MicroRNAs"[Mesh] OR "miRNAs"[Mesh] OR "microRNA"[Title/Abstract] OR "miRNA"[Title/Abstract]) 
AND ("Coronary Artery Disease"[Mesh] OR "Coronary Artery Disease" )) 
AND 2000:2022[DP]'
query_list[["DCM"]] <- '(("MicroRNAs"[Mesh] OR "miRNAs"[Mesh] OR "microRNA"[Title/Abstract] OR "miRNA"[Title/Abstract]) 
AND ("Cardiomyopathy, Dilated"[Mesh] OR "Cardiomyopathy, Dilated" OR "dilated cardiomyopathy")) 
AND 2000:2022[DP]'
query_list[["HFrEF"]] <- '(("MicroRNAs"[Mesh] OR "miRNAs"[Mesh] OR "microRNA"[Title/Abstract] OR "miRNA"[Title/Abstract]) 
AND ("Heart Failure"[Mesh] OR "Heart Failure, Systolic"[Mesh] OR "Heart Failure, Systolic" OR "heart failure with reduced ejection fraction")) 
AND 2000:2022[DP]'

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
    guides(fill = "none") +
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    ylab("# of articles") + 
    scale_y_continuous(expand = c(0,0))
  return(count_plot)
}

docl <- list()  # create a list object for storing results

for (disease in seq_along(diseases)) {
  file <- glue::glue("./data-literature/pubmed-MicroRNA-{diseases[disease]}.txt")
  
  # load abstracts ----------------------------------------------------------
  df <- read_pubmed(file, topic = diseases[disease]) %>% 
    # Subset for original articles
    subset_research() %>%
    # Extract miRNA names mentioned at least twice/abstract
    extract_mir_df(threshold = 2)
  
  # Number of abstracts
  length(unique(df[["PMID"]]))
  # Number of rows
  nrow(df)
  
  # Compare extraction count to PubMed --------------------------------------
  # Extract all miRNAs
  df_comparison <- 
    read_pubmed(file,
                topic = diseases[disease]) %>% 
    extract_mir_df(threshold = 1,
                   extract_letters = TRUE)
  
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
    filter(Mentioned_n >= 10) %>% 
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
  count_pubmed_vec <- purrr::map(df_count_no_letter$miRNA, ~compare_count(mir = .x, mesh_query = query_list[[diseases[disease]]])) %>% 
    purrr::set_names(df_count_no_letter$miRNA)
  
  count_pubmed_vec_letter <- purrr::map(df_count_letter$miRNA, ~compare_count(mir = .x, mesh_query = query_list[[diseases[disease]]])) %>% 
    purrr::set_names(df_count_letter$miRNA)
  
  # Figure 1: Compare miRetrieve with PubMeds miRNA count -------------------
  
  ## a) no letter
  # Obtain number of PubMed results
  length_pmid <- purrr::map(count_pubmed_vec, length) %>% 
    unlist()
  
  # Add number of PubMed results to dataframe "df_count_no_letter"
  df_count_no_letter[3] <- length_pmid
  names(df_count_no_letter)[3] <- "PubMed"
  names(df_count_no_letter)[2] <- "miRetrieve"
  # Plot # of miRNAs with miRetrieve vs. PubMed (no letters)
  count_plot1 <- count_plot_fn(df_count = df_count_no_letter)
  
  ## b) with letter
  # Obtain number of PubMed results
  length_pmid_letter <- purrr::map(count_pubmed_vec_letter, length) %>% 
    unlist()
  
  # Add number of PubMed results to dataframe "df_count_letter"
  df_count_letter[3] <- length_pmid_letter
  names(df_count_letter)[3] <- "PubMed"
  names(df_count_letter)[2] <- "miRetrieve" 
  # Plot # of miRNAs with miRetrieve vs. PubMed (no letters)
  count_plot1_letter <- count_plot_fn(df_count = df_count_letter)
  
  # Combine plots (Figure 1)
  combined <- count_plot + count_letter_plot & theme(legend.position = "bottom")
  
  combined <- combined + 
    plot_layout(guides = "collect") + 
    plot_annotation(tag_levels = "A")
  
  if (saveFILE == TRUE) {
    ## arranged plots are toooo small 
    filename.fig1 <- glue::glue("./output/plots/miRetrieve/{diseases[disease]}/{Sys.Date()}-figure1-compare")
    ggsave(filename = paste0(filename.fig1, ".svg"), plot =  combined, 
           width = 8, height = 8, 
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
  
  # Figure 2: Top miRNAs in {Disease} ---------------------------------------
  
  most_frequently_miRNA <- plot_mir_count(df,
                                          title = glue::glue("Most frequently mentioned miRNAs in {diseases[disease]}"))
  
  if (saveFILE == TRUE) {
    ## arranged plots are toooo small 
    filename.fig2 <- glue::glue("./output/plots/miRetrieve/{diseases[disease]}/{Sys.Date()}-figure2-mostfrequent")
    ggsave(filename = paste0(filename.fig2, ".svg"), plot =  most_frequently_miRNA, 
           width = 8, height = 8, 
           units = "in"  # default
    )
  }
  
  # Figure 7: Potential miRNA biomarker in diseases --------------------
  # pinpointed abstracts containing at least five marker 
  #   words for biomarker(bio-marker, biological marker, 
  #   biomarker, body fluid,bodyfluid, circulating, diagnostic, 
  #   exosomal, exosomes,extracellular vesicles, plasma, serum, urinary, urine). 
  
  # Potential miRNA biomarker in ACS
  disease_biomarker <- calculate_score_biomarker(df,
                                                 biomarker_keywords[-c(5,6,9,10)],  # no urine, body fluid
                                                 threshold = 5,
                                                 discard = TRUE)
  # Plot top 7 miRNA biomarker in DCM
  biomarker_plot <- plot_mir_count(disease_biomarker, 
                                   top = 8,
                                   title = glue::glue("Potential biomarker miRNAs in {diseases[disease]}"))
  if (saveFILE == TRUE) {
    ## arranged plots are toooo small 
    filename.fig7 <- glue::glue("./output/plots/miRetrieve/{diseases[disease]}/{Sys.Date()}-figure7-mostfrequent")
    ggsave(filename = paste0(filename.fig7, ".svg"), plot =  biomarker_plot, 
           width = 8, height = 8, 
           units = "in"  # default
    )
  }
}



# not included in loop ----------------------------------------------------

# Figure 3: Top terms {disease} -------------------------------------------
# Generate new stop word list specifically for atherosclerosis to 
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




# also included in loop ---------------------------------------------------

# Figure 7: Potential miRNA biomarker in both diseases --------------------
# pinpointed abstracts containing at least five marker 
#   words for biomarker(bio-marker, biological marker, 
#   biomarker, body fluid,bodyfluid, circulating, diagnostic, 
#   exosomal, exosomes,extracellular vesicles, plasma, serum, urinary, urine). 

# Potential miRNA biomarker in ACS
acs_biomarker <- calculate_score_biomarker(df_acs,
                                           biomarker_keywords[-c(5,6,9,10)],  # no urine, body fluid
                                           threshold = 5,
                                           discard = TRUE)


# Plot top 7 miRNA biomarker in ACS
bm_acs_plot <- plot_mir_count(acs_biomarker, 
                              top = 8,
                              title = "Potential biomarker miRNAs in ACS")

# Potential miRNA biomarker in DCM
dcm_biomarker <- calculate_score_biomarker(df_dcm,
                                           threshold = 5,
                                           discard = TRUE)

# Plot top 7 miRNA biomarker in DCM
bm_dcm_plot <- plot_mir_count(dcm_biomarker, 
                              top = 8,
                              title = "Potential biomarker miRNAs in DCM")

# Combine plot (Figure 7)
(bm_acs_plot) + (bm_dcm_plot) + plot_annotation(tag_levels = 'A')



