


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
  .libPaths("/mnt/users/reich/programs/R43/lib")
  lib_path <- "/mnt/users/reich/programs/R43/lib" 
  data_path_bestageing2022 <- "/mnt/users/reich/rockerprojects/bestageing2022"
  data_path_BestAgeing <- "/mnt/users/reich/BestAgeing"
}


library(rentrez, lib.loc = lib_path)
library(pubmedR, lib.loc = lib_path)
library(dplyr, lib.loc = lib_path)
library(tidyr, lib.loc = lib_path)
library(ggplot2, lib.loc = lib_path)
library(ggpubr, lib.loc = lib_path)

source(file = glue("{data_path_bestageing2022}/scripts/helper/custom_ggplot_theme.R"))

api_key <- Sys.getenv("PUBMED_API")

diseases <- c("dcm","acs", "cad", "hfref")

## search MeSH terms -------------------------------------------------------
# https://www.ncbi.nlm.nih.gov/mesh

# precise query -------------------------------------------------------------------
query <- '(("MicroRNAs"[Mesh] OR "miRNAs"[Mesh] OR "microRNA"[Title/Abstract] OR "miRNA"[Title/Abstract]) 
AND ("Cardiovascular Diseases"[Mesh] OR "Cardiovascular Diseases"[Title/Abstract] OR "Heart Diseases"[Mesh] OR "Heart Diseases"[Title/Abstract])) 
AND 2000:2022[DP]'
# HELP: write a query: https://pubmed.ncbi.nlm.nih.gov/help/#search-tags or use ChatGPT ;)

query_acs <- '(("MicroRNAs"[Mesh] OR "miRNAs"[Mesh] OR "microRNA"[Title/Abstract] OR "miRNA"[Title/Abstract]) 
AND ("Acute Coronary Syndrome"[Mesh] OR "Acute Coronary Syndrome" OR "Myocardial Infarction"[Mesh])) 
AND 2000:2022[DP]'

query_cad <- '(("MicroRNAs"[Mesh] OR "miRNAs"[Mesh] OR "microRNA"[Title/Abstract] OR "miRNA"[Title/Abstract]) 
AND ("Coronary Artery Disease"[Mesh] OR "Coronary Artery Disease" )) 
AND 2000:2022[DP]'

query_hfref <- '(("MicroRNAs"[Mesh] OR "miRNAs"[Mesh] OR "microRNA"[Title/Abstract] OR "miRNA"[Title/Abstract]) 
AND ("Heart Failure"[Mesh] OR "Heart Failure, Systolic"[Mesh] OR "Heart Failure, Systolic" OR "heart failure with reduced ejection fraction")) 
AND 2000:2022[DP]'

query_dcm <- '(("MicroRNAs"[Mesh] OR "miRNAs"[Mesh] OR "microRNA"[Title/Abstract] OR "miRNA"[Title/Abstract]) 
AND ("Cardiomyopathy, Dilated"[Mesh] OR "Cardiomyopathy, Dilated" OR "dilated cardiomyopathy")) 
AND 2000:2022[DP]'

# run queries for each CV disease separately
query_loop_df <- tibble(diseases = diseases, query = c(query_dcm,query_acs, query_cad, query_hfref ))

### using {rentrez} ---------------------------------------------------------
# https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html#web_history
entrez_dbs()
entrez_db_summary("pubmed")
entrez_db_searchable("pubmed")  #  find out which search terms can be used

r_search <- entrez_search(db="pubmed", term=query, retmax=0) # retmax max 9999, only EDirect (cmd line tool can obtain more queries, but we could e.g. split time)
r_search

# advanced counting
search_year <- function(year, term){
  query <- paste(term, "AND (", year, "[PDAT])")
  entrez_search(db="pubmed", term=query, retmax=0)$count
}

year <- 2000:2022
# run function searching for research in ALL CV research
papers <- sapply(year, search_year, term=query, USE.NAMES=FALSE)
paper_count_df <- tibble(year = year, papers = papers)
year_annotation <- c(2005, 2010, 2015, 2020)

ggplot(mapping=aes(x=year, y=papers), data = paper_count_df %>% filter(year>=2005 & year <=2020))+
  geom_point(alpha =0.3)+
  geom_line(linewidth=1)+
  #geom_smooth(method="auto", se=FALSE, fullrange=FALSE, level=0.95) +  # plot smoothing line
  labs(title = "The Rise of Cardiovascular miRNA Research",
       #subtitle = subtitle.custom,
       caption = '("MicroRNAs"[Mesh] OR "miRNAs"[Mesh] OR "microRNA"[Title/Abstract] OR "miRNA"[Title/Abstract]) 
AND ("Cardiovascular Diseases"[Mesh] OR "Cardiovascular Diseases"[Title/Abstract] OR 
"Heart Diseases"[Mesh] OR "Heart Diseases"[Title/Abstract]) AND 2000:2022[DP]') +
  scale_x_continuous(name = NULL,  
                     #labels = comma, 
                     breaks = year_annotation, #c( 10000, n[7:length(n)]),
                     #limits = c(10000, 300000)
                     )+
  ylab("Count of Papers")+
  theme_minimal(base_size = 16, base_family = 'Arial') +
  labs(title = NULL) +  # update title to be NULL
  scale_fill_manual(values = thematic::okabe_ito(6)) +
  scale_colour_manual(values = thematic::okabe_ito(6)) +
  my_base_theme() -> paper.count.plot1
paper.count.plot1

paper_count_df_diseases <- tibble(year = year)
for (my_disease_query in 1:nrow(query_loop_df)){
  papers <- sapply(year, search_year, term=query_loop_df$query[my_disease_query] , USE.NAMES=FALSE)
  varname <- paste0("papers_",query_loop_df$diseases[my_disease_query])
  paper_count_df_diseases <- paper_count_df_diseases %>% 
    mutate( !!varname := papers )
  print(paste0("|||-----------------------Run finished for disease: ", toupper(query_loop_df$diseases[my_disease_query]), " -----------------------|||"))
}
paper_count_df_diseases1 <- paper_count_df_diseases  # copy bc creating takes some time
## now plot
paper_count_df_diseases1 %>% 
  pivot_longer(cols = !year, names_to = "Disease") %>% 
  mutate(Disease = factor(Disease, labels = c("ACS", "CAD", "DCM", "HF"))) %>% 
  filter(year>=2005 & year <=2020) %>% 
  # plot
  ggplot(mapping=aes(x=year, y=value, color=Disease)) +
  geom_point(alpha =0.3) + 
  geom_line(linewidth=1) +
  labs(#title = "Papers published in distinct phenotypes",
       subtitle = NULL,  #"Papers published in Distinct Phenotypes"
       ) +
  xlab("Year")+
  ylab("Count of Papers")+
  theme_minimal(base_size = 16, base_family = 'Arial') +
  scale_fill_manual(values = thematic::okabe_ito(6)) +
  scale_colour_manual(values = thematic::okabe_ito(6)) +
  my_base_theme() -> paper.count.plot2
paper.count.plot2

advanced.counting.plot <- ggarrange(paper.count.plot1, paper.count.plot2,
                  ncol=1, nrow=2, 
                  #heights = c(2, 1.5, 1.5),
                  legend = "bottom", common.legend = TRUE, align = "hv"
)
advanced.counting.plot

saveRDS(advanced.counting.plot, file = "./output/pubmed/paper_counting_plot.rds")
# advanced.counting.plot <- readRDS("./output/pubmed/paper_counting_plot.rds")

ggsave(filename = glue("{data_path_bestageing2022}/output/pubmed/{Sys.Date()}_paper_counting_plot.svg"), plot = advanced.counting.plot,
       width = 10, 
       height = 10)
  
# ### using {pubmedR}
# # https://cran.r-project.org/web/packages/pubmedR/vignettes/A_Brief_Example.html 
# 
# # Check the effectiveness of the query
# res <- pmQueryTotalCount(query = query, api_key = api_key, limit = 9998) # cannot be larger than 9998. For PubMed, ESearch can only retrieve the first 9,999 records matching the query.
# # To obtain more than 9,999 PubMed records, consider using EDirect https://www.ncbi.nlm.nih.gov/books/NBK25499/ (https://dataguide.nlm.nih.gov/edirect/documentation.html)
# 
# res$total_count
# 
# # Download the collection of document metadata
# D <- pmApiRequest(query = query, limit = res$total_count, api_key = api_key)  # 
# 
# #  Convert the download object into a “readable” and and “usable” format
# M <- pmApi2df(D)
# str(M)
