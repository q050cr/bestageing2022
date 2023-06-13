
library(rentrez)
library(pubmedR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

api_key <- "a585126d33c778791be5fb8ce8fb3d066408"

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

ggplot(mapping=aes(x=year, y=papers), data = paper_count_df %>% filter(year>=2005))+
  geom_point(alpha =0.3)+
  geom_line()+
  #geom_smooth(method="auto", se=FALSE, fullrange=FALSE, level=0.95) +  # plot smoothing line
  labs(title = "The Rise of Cardiovascular miRNA Research",
       #subtitle = subtitle.custom,
       caption = '("MicroRNAs"[Mesh] OR "miRNAs"[Mesh] OR "microRNA"[Title/Abstract] OR "miRNA"[Title/Abstract]) 
AND ("Cardiovascular Diseases"[Mesh] OR "Cardiovascular Diseases"[Title/Abstract] OR 
"Heart Diseases"[Mesh] OR "Heart Diseases"[Title/Abstract]) AND 2000:2022[DP]') +
  scale_x_continuous(name = "Year",  
                     #labels = comma, 
                     breaks = year_annotation, #c( 10000, n[7:length(n)]),
                     #limits = c(10000, 300000)
                     )+
  ylab("Count of Papers")+
  ggthemes::theme_few() +
  ggthemes::scale_color_few() -> paper.count.plot1

paper_count_df_diseases <- tibble(year = year)
for (my_disease_query in 1:nrow(query_loop_df)){
  papers <- sapply(year, search_year, term=query_loop_df$query[my_disease_query] , USE.NAMES=FALSE)
  varname <- paste0("papers_",query_loop_df$diseases[my_disease_query])
  paper_count_df_diseases <- paper_count_df_diseases %>% 
    mutate( !!varname := papers )
  print(paste0("|||-----------------------Run finished for disease: ", toupper(query_loop_df$diseases[my_disease_query]), " -----------------------|||"))
}

## now plot
paper_count_df_diseases %>% 
  pivot_longer(cols = !year, names_to = "Disease") %>% 
  mutate(Disease = factor(Disease, labels = c("ACS", "CAD", "DCM", "HF"))) %>% 
  ggplot(mapping=aes(x=year, y=value, color=Disease)) +
  geom_point(alpha =0.3) + 
  geom_line() +
  labs(#title = "Papers published in distinct phenotypes",
       subtitle = "Papers published in Distinct Phenotypes") +
  xlab("Year")+
  ylab("Count of Papers")+
  ggthemes::theme_few() +
  ggthemes::scale_color_few() -> paper.count.plot2

advanced.counting.plot <- ggarrange(paper.count.plot1, paper.count.plot2,
                  ncol=1, nrow=2, 
                  #heights = c(2, 1.5, 1.5),
                  legend = "bottom", common.legend = TRUE, align = "hv"
)

saveRDS(advanced.counting.plot, file = "./output/pubmed/paper_counting_plot.rds")
# advanced.counting.plot <- readRDS("./output/pubmed/paper_counting_plot.rds")

ggsave(filename = "./output/pubmed/paper_counting_plot.png", plot = advanced.counting.plot,
       width = 14, 
       height = 14)
  
# ### using {pubmedR} ---------------------------------------------------------
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
