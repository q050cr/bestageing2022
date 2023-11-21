
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

# libs
require(readxl, lib.loc = lib_path)
require(janitor, lib.loc = lib_path)
require(glue, lib.loc = lib_path)
require(gt, lib.loc = lib_path)
require(dplyr, lib.loc = lib_path)
require(tidyr, lib.loc = lib_path)
require(stringr, lib.loc = lib_path)
require(purrr, lib.loc = lib_path)
require(dplyr, lib.loc = lib_path)
require(e1071, lib.loc = lib_path)
require(ggplot2, lib.loc = lib_path)
require(ggridges, lib.loc = lib_path)
require(pROC, lib.loc = lib_path)
require(htmltools, lib.loc = lib_path)
require(webshot, lib.loc = lib_path)

source(file = glue("{data_path_bestageing2022}/scripts/helper/custom_ggplot_theme.R"))

# vars
diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>% 
  filter(analysis=="full")



# created in 03miRetrieve_topmirnas_all_diseases.R 
table_mirna_top50bm_score_alldiseases <- readRDS(file = glue("{data_path_bestageing2022}/data-literature/miRetrieve/2023-07-27_top50mirnas_all_diseases_pmids_gpt.rds"))  # new

# init
prepare_gt_dat <- tibble(Topic=NA, TargetName=NA, Accession=NA, Biomarker_score=NA, miRetrieve=NA, padj=NA, padj.glm=NA, sign_indicator=NA, boxplots=NA, rocaucs=NA, PMIDs=NA)

for(mydisease in 1:nrow(all_combis[all_combis[["analysis"]] == "full", ])){
  # created in model_de_analysis.R
  data01 <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/data01.rds"))
  de_results <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/de_results.rds"))
  results_logmedians <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/results_logmedians.rds"))
  
  # need to specify plot function here because of change in diseases
  ggboxplot_table <- function(plotdata, miRNA_name, plot_disease){
    ggplot(data=plotdata, aes(x= {{plot_disease}}, y= .data[[miRNA_name]] )) +
      geom_boxplot(outlier.shape = NA, alpha=0.7) +
      geom_jitter(aes(color={{plot_disease}}), position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), size=0.1, alpha=0.4) +
      # prevent geom_text() from searching for the disease variable, you can override the fill aesthetic by setting fill = NULL within the aes() function 
      labs(
        x = NULL, 
        y = NULL,  #expression(paste("log"[2], " expression")), 
        title = NULL) +
      #scale_y_continuous(breaks = seq(0, ceiling(max(miRNA_name)), by = 2), #limits = c(0, 1), expand = c(0, 0)) +
      # theme_classic()+
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = thematic::okabe_ito(6), labels = c(toupper(all_combis$diseases[mydisease]), "Control")) + guides(fill=guide_legend(title=NULL), color="none") +
      scale_color_manual(values = thematic::okabe_ito(6)) +
      #scale_color_brewer(palette = "Set1") + # Choose a color palette
      #scale_fill_brewer(palette = "Set1", labels = c(toupper(all_combis$diseases[mydisease]), "Control")) + guides(fill=guide_legend(title=NULL), color="none")+
      my_base_theme() +
      theme(axis.text = element_blank(), axis.ticks = element_blank())  # remove axis text and ticks
  }
  
  # univariate AUCs
  ggrocplot_table <- function(plotdata, miRNA_name) {
    roc_obj <- suppressMessages(roc(data01[["disease"]], data01[[miRNA_name]]))
    roc_df <- data.frame(
      TPR = roc_obj$sensitivities,
      FPR = 1 - roc_obj$specificities
    )
    
    # Create a ggplot2 ROC curve
    ggplot(roc_df, aes(x = FPR, y = TPR)) +
      geom_line() +
      labs(x = NULL, y = NULL, title = NULL) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      annotate("text", x = .6, y = .25, label = paste("AUC =", round(auc(roc_obj), 2)), size=16) +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      my_base_theme() +
      coord_equal() +
      theme(axis.text = element_blank(), axis.ticks = element_blank())
  }
  
  filtered4disease <- table_mirna_top50bm_score_alldiseases %>% 
    filter(toupper(Topic) == toupper(all_combis$diseases[mydisease])) %>% 
    left_join(results_logmedians, by=c("TargetName"="miR")) %>% 
    arrange(desc(auc)) 
  
  chosen_mirnas <- filtered4disease$TargetName
  # plot_list_boxplots <- purrr::map(.x = chosen_mirnas, .f = ~ ggboxplot_table(plotdata = data01, miRNA_name = .x, plot_disease = disease) )
  plot_list_rocaucs <- purrr::map(.x = chosen_mirnas, .f = ~ ggrocplot_table(plotdata = data01, miRNA_name = .x) )
  
  # with mutating the plots to existing df always computer hanged
  prepare_gt_dat_tmp <- tibble(TargetName = chosen_mirnas,
                               #boxplots = plot_list_boxplots,   # looads of space; size in MB: object.size(prepare_gt_dat_tmp) / 1024^2
                               rocaucs = plot_list_rocaucs) %>% 
    left_join(filtered4disease %>% select(Topic, TargetName, Accession, Biomarker_score, miRetrieve, max_value, padj, padj.glm, sign_indicator, auc, PMIDs), 
              by= c("TargetName"="TargetName"))
  
  prepare_gt_dat <- bind_rows(prepare_gt_dat, prepare_gt_dat_tmp)
}


# gt table beginning -----------------------------------------------------------

# Define custom order
custom_order <- c("Acute Coronary Syndrome", "Coronary Artery Disease", "Dilated Cardiomyopathy", "Heart Failure With Reduced Ejection Fraction")

prepare_gt_dat <- prepare_gt_dat[-1, ] 
# for gt to work the column in the tibble can be NA, but the plots should be stored in a list (which we can index from the old tibble)
prepare_gt_dat1 <- prepare_gt_dat %>% select(-boxplots)
prepare_gt_dat1 <- prepare_gt_dat1 %>% mutate(rocaucs=NA)

gt_plot_list <- prepare_gt_dat %>% 
  mutate(
    Topic = case_when(
      Topic == "ACS" ~ "Acute Coronary Syndrome",
      Topic == "CAD" ~ "Coronary Artery Disease",
      Topic == "DCM" ~ "Dilated Cardiomyopathy",
      Topic == "HFrEF" ~ "Heart Failure With Reduced Ejection Fraction",
    ),
    Topic = factor(Topic, levels = (custom_order))  # consistent with the behavior of ggplot2, where factor levels are plotted in reverse order to match the standard layout of a Cartesian coordinate system
  ) %>% 
  group_by(Topic) %>% 
  # THIS ALLOWS THE ORDERING OF THE SPANNERS! 
  arrange(Topic, desc(auc), desc(max_value)) %>% 
  slice(1:5) %>%   # if sliced list of plots MUST be adapted !!done here
  select(Topic, TargetName, rocaucs, auc)
  
  

prepare_gt_dat1 %>% 
  mutate(
    TargetName = str_replace_all(TargetName, "_", "-"), # replace underscores with hyphens
    padj = round(padj, 3), # round to 3 decimal places
    padj = ifelse(padj == 0, "≤ 0.001", as.character(padj)),
    padj.glm = round(padj.glm, 3), # round to 3 decimal places
    padj.glm = ifelse(padj.glm == 0, "≤ 0.001", as.character(padj.glm))
  ) %>% 
  mutate(
    Topic = case_when(
      Topic == "ACS" ~ "Acute Coronary Syndrome",
      Topic == "CAD" ~ "Coronary Artery Disease",
      Topic == "DCM" ~ "Dilated Cardiomyopathy",
      Topic == "HFrEF" ~ "Heart Failure With Reduced Ejection Fraction",
    ),
    Topic = factor(Topic, levels = (custom_order))  # consistent with the behavior of ggplot2, where factor levels are plotted in reverse order to match the standard layout of a Cartesian coordinate system
  ) %>% 

  mutate(
    PMIDs = paste0("<a href='http://www.ncbi.nlm.nih.gov/pubmed/",
                   PMIDs,
                   "' style='font-size:60%;'>",
                   PMIDs,
                   "</a>")
  ) %>% 
  group_by(Topic) %>% 
  # THIS ALLOWS THE ORDERING OF THE SPANNERS! 
  arrange(Topic, padj, desc(max_value)) %>%  # arrange(Topic, desc(auc), desc(max_value))
  slice(1:5) %>%   # only top 5
  relocate(
    max_value, .after = Accession
  ) %>% 
  select(-c(auc, Biomarker_score, miRetrieve, padj.glm)) %>% 
  
  ## start gt -------------------------------------------------------------
  gt(groupname_col = "Topic") %>% 
  opt_table_font(font = "Arial") %>%
  tab_header(title = md("**Literature miRNAs from miRetrieve**")) %>%
  fmt_markdown(columns = c(PMIDs)) %>%
  ## groupings spanner label
  #tab_spanner(
  #  label = "Topic", 
  #  columns = c(Topic)
  #) %>% 
  cols_label(
    #Topic = "Topic",
    TargetName = "microRNA",
    Accession = "Accession",
    max_value = "Biomarker Score",
    #miRetrieve = "miRetrieve",
    padj = "p-adj (Holm)",
    #padj.glm = "p-adj-glm (BH)",
    sign_indicator = "Significance Indicator",
    rocaucs = "ROC AUCs", 
    PMIDs = "PMIDs"
  ) %>%
  tab_spanner(
    label = "Performance",
    columns = c(
      max_value, padj,sign_indicator, rocaucs
    )
  ) %>% 
  #tab_spanner_delim(delim = "_") %>%
  data_color(
    columns = c(max_value),
    fn = scales::col_numeric(
      palette = c(thematic::okabe_ito(6)[1], thematic::okabe_ito(6)[2]),
      domain = NULL
    )
  ) %>% 
  text_transform(
    locations = cells_body(columns = rocaucs),
    fn = function(x) {
      gt_plot_list$rocaucs |>   # cave if ordering changed!!
        ggplot_image(height = px(50))  # had 100 before
      #prepare_gt_dat$rocaucs |>
      #  ggplot_image(height = px(100))
    }
  ) %>% 
  #fmt_number(
  #  # columns = 3:22,
  #  decimals = 3, 
  #  use_seps=TRUE
  #  ) %>% 
  tab_source_note(
    source_note = md("**Biomarker Score** as calculated by *{miRetrieve}*. P values 
                     were adjusted using the Benjamini-Hochberg procedure. 
                     The variable p.adj.glm represents p values from t-tests 
                     adjusted for age and sex.")
  ) %>% 
  tab_footnote(
    footnote = "PMID = PubMed Unique Identifier",
    locations = cells_column_labels(columns = PMIDs)
  ) %>% 
  opt_stylize(style=1, color="gray") %>% 
  tab_options(
    table.width = pct(100),
    footnotes.multiline = FALSE,
    footnotes.marks = letters,
    data_row.padding = px(0)
    ) %>% 
  # Align text to center
  cols_align(
    align = c("center"),
    columns = c(everything())
  ) %>% 
  # Adjust width of cols
  cols_width(
    TargetName ~ px(100),
    Accession ~ px(100),
    rocaucs ~ px(100), 
    sign_indicator ~ px(100),
    starts_with("p.adj") ~ px(40),
    PMIDs ~ px(150),
    everything() ~ px(50)
  ) %>% 
  tab_options(
    table.font.size = px(10L)
  ) -> gt_literature_miRetrieve

gt_literature_miRetrieve

if(SAVE.files == TRUE) {
  filename_gt_miretrieve <- glue("{data_path_bestageing2022}/output/tables/literature_miRetrieve/{Sys.Date()}_table_miRetrieve_miRNA_top5")
  gtsave(data = gt_literature_miRetrieve, glue("{filename_gt_miretrieve}.html")) 
  # gtsave(data = gt_literature_miRetrieve, glue("{filename_gt_miretrieve}.docx")) 
}



# Assuming you have a gt table named 'gt_table'
gt_literature_miRetrieve %>%
  as_raw_html() %>%
  html_print() %>%
  webshot(glue("{filename_gt_miretrieve}.pdf"))  # Should end with .png, .pdf, or .jpeg

gt_literature_miRetrieve %>%
  as_raw_html() %>%
  html_print() %>%
  webshot(glue("{filename_gt_miretrieve}.png"))

browsable(as.tags(gt_literature_miRetrieve %>% as_raw_html()))


# descriptive top 5 -----------------------------------------------------------

prepare_gt_dat1 %>% 
  mutate(
    TargetName = str_replace_all(TargetName, "_", "-"), # replace underscores with hyphens
    padj = round(padj, 3), # round to 3 decimal places
    padj = ifelse(padj == 0, "≤ 0.001", as.character(padj)),
    padj.glm = round(padj.glm, 3), # round to 3 decimal places
    padj.glm = ifelse(padj.glm == 0, "≤ 0.001", as.character(padj.glm))
  ) %>% 
  mutate(
    Topic = case_when(
      Topic == "ACS" ~ "Acute Coronary Syndrome",
      Topic == "CAD" ~ "Coronary Artery Disease",
      Topic == "DCM" ~ "Dilated Cardiomyopathy",
      Topic == "HFrEF" ~ "Heart Failure With Reduced Ejection Fraction",
    ),
    Topic = factor(Topic, levels = (custom_order))  # consistent with the behavior of ggplot2, where factor levels are plotted in reverse order to match the standard layout of a Cartesian coordinate system
  ) %>% 
  
  mutate(
    PMIDs = paste0("<a href='http://www.ncbi.nlm.nih.gov/pubmed/",
                   PMIDs,
                   "' style='font-size:60%;'>",
                   PMIDs,
                   "</a>")
  ) %>% 
  group_by(Topic) %>% 
  # THIS ALLOWS THE ORDERING OF THE SPANNERS! 
  arrange(Topic, desc(auc), desc(max_value)) %>% 
  #slice(1:5) %>%   # only top 5
  relocate(
    max_value, .after = Accession
  ) -> describe_mirna_lit

describe_mirna_lit %>% 
  select(sign_indicator) %>% 
  table() %>% 
  as_tibble() -> sign_count_df
sign_count_df

## barplot significance -----------------

ggplot(sign_count_df, aes(x = Topic, y = n, fill = sign_indicator)) +
  geom_bar(stat = "identity", position = "dodge", alpha=0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c("ACS", "CAD", "DCM", "HFrEF")) +
  labs(
    x = NULL,
    y = "Count",
    title = NULL,
    subtitle = NULL, 
    fill = NULL
  )+
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_fill_manual(values = thematic::okabe_ito(6)) +
  scale_color_manual(values = thematic::okabe_ito(6))  -> sign_indicator_barplot

sign_indicator_barplot

# AUC distribution
describe_mirna_lit %>% 
  select(auc) -> auc_distribution_topic_df

# Custom labels for the legend
custom_labels <- c(
  "Acute Coronary Syndrome" = "ACS",
  "Coronary Artery Disease" = "CAD",
  "Dilated Cardiomyopathy" = "DCM",
  "Heart Failure With Reduced Ejection Fraction" = "HFrEF"
)

## densities auc -----------------

ggplot(auc_distribution_topic_df, aes(x = auc, fill = Topic)) +
  geom_density(alpha = 0.5) +
  labs(
    x = "AUC",
    y = "Density",
    title = "Density Plot of AUC Values",
    subtitle = "Grouped by Topic"
  ) +
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_fill_manual(values = thematic::okabe_ito(6), labels = custom_labels) +
  scale_color_manual(values = thematic::okabe_ito(6)) 

ggplot(auc_distribution_topic_df, aes(x = auc, y = Topic, fill = Topic)) +
  geom_density_ridges(alpha = 0.7) +
  labs(
    x = "Univariate AUC",
    y = NULL,
    title = NULL,  # "Density Ridges Plot of AUC Values"
    subtitle = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_fill_manual(values = thematic::okabe_ito(6), labels = custom_labels) +
  scale_color_manual(values = thematic::okabe_ito(6)) +
  theme(axis.text.y = element_blank()) -> density_ridges_plot

density_ridges_plot

## arranged final plot with venn!!  -----------------

arranged_plot <- ggarrange(sign_indicator_barplot, density_ridges_plot, nrow = 2, align = "hv", labels = c("A", "B"))

# venn plot high from 001_model_de_analysis.R
venn_plot_high <- readRDS(glue("{data_path_bestageing2022}/output/plots/venn_dia/venn_de_mirnas_shared.rds"))
arranged_plot <- ggarrange(arranged_plot, venn_plot_high, nrow = 1, ncol = 2, align = "hv", labels = c("", "C"))
arranged_plot

if(SAVE.files == TRUE) {
  filename_plots_top50_pval_auc <- glue("{data_path_bestageing2022}/output/plots/miRetrieve/stats/top50_miRetrieve_stats_sign_aucs.svg")
  ggsave(filename = filename_plots_top50_pval_auc, plot =  arranged_plot, 
         width = 10, height = 8, 
         units = "in"  # default
  )
}
