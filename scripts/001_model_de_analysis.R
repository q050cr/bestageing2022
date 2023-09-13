
### INFO ----------------------------------------------------------------------
# Differential miRNA Expression Analysis
# this script is sourced from `scripts/render_param_reports.R`
# selection provided by `all_combis$diseases` and `all_combis$analysis`

# script creates plots: "fig01vogel2013", "fig02vogel2013"


# Get system name
system_name <- Sys.info()["nodename"]
mount_filesystem <- TRUE
SAVE.files <- TRUE
runTests <- TRUE

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

# dependencies ---------------------------------------------------------------
library(readxl, lib.loc = lib_path)
library(janitor, lib.loc = lib_path)
library(glue, lib.loc = lib_path)
library(gt, lib.loc = lib_path)
library(dplyr, lib.loc = lib_path)
library(tidyr, lib.loc = lib_path)
library(stringr, lib.loc = lib_path)
library(purrr, lib.loc = lib_path)
library(dplyr, lib.loc = lib_path)
library(ggplot2, lib.loc = lib_path)
library(gplots, lib.loc = lib_path)  # for Venn plot/ data intersections
library(gridExtra, lib.loc = lib_path)  # for adding marginal density plots
library(cowplot, lib.loc = lib_path)  # for adding marginal density plots
library(RColorBrewer, lib.loc = lib_path)
library(ggdist, lib.loc = lib_path)
library(gghalves, lib.loc = lib_path)
library(ggrepel, lib.loc = lib_path)
library(ggvenn, lib.loc = lib_path)
library(rstatix, lib.loc = lib_path)
library(ggthemes, lib.loc = lib_path)
library(ggpubr, lib.loc = lib_path)
library(pROC, lib.loc = lib_path)
library(reshape2, lib.loc = lib_path)
conflicted::conflict_prefer("expand", "tidyr")
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(janitor::make_clean_names)
conflicted::conflicts_prefer(pROC::roc)

source(file = glue("{data_path_bestageing2022}/scripts/helper/custom_ggplot_theme.R"))



diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) 


diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>% 
  filter(analysis=="full")


# load data ---------------------------------------------------------------

# MIRNA DAT
model_data1 <- clean_names(readRDS(file = glue('{data_path_BestAgeing}/data_new/model_data1.RDS')))  # has also multiclass col + diagnoses
load(file = glue('{data_path_BestAgeing}/data/mirnas.rda'))  # "UKL-HD" n=765
load(file = glue('{data_path_BestAgeing}/data/data.rda'))  # "UKL-HD" n=731
all_mirnas <- clean_names(mirnas)
names(all_mirnas) <- str_replace(string = names(all_mirnas), pattern = "mi_r", replacement = "mir")
mirnas_disease <- clean_names(data)
names(mirnas_disease) <- str_replace(string = names(mirnas_disease), pattern = "mi_r", replacement = "mir")
rm(data, mirnas)

# load-mirnas-from_research
# create vector of described mirnas
load(glue("{data_path_BestAgeing}/data_research/fromR/researchMiRNAAccession.rda"))
# check if all mirnas are named the same
researchMiRNAAccession$miRNAName_v21 <-  make_clean_names(researchMiRNAAccession$miRNAName_v21) %>% 
  str_replace(pattern = "mi_r", replacement = "mir")
# researchMiRNAAccession$miRNAName_v21[(!researchMiRNAAccession$miRNAName_v21 %in% colnames(all_mirnas))] 
## "hsa_mir_106a_5p" --> should be --> "hsa_mir_106b_5p"
## all_mirnas[,(str_detect(string = colnames(all_mirnas), pattern = "106"))]
## researchMiRNAAccession[which(!researchMiRNAAccession$miRNAName_v21 %in% colnames(all_mirnas)),]
researchMiRNAAccession$miRNAName_v21[researchMiRNAAccession$miRNAName_v21 == "hsa_mir_106a_5p"] <- "hsa_mir_106b_5p"

###
# load-metadat --
## DIAGNOSES DAT
load(glue("{data_path_BestAgeing}/data/diagnoses_df.rda"))

## SURVIVAL DAT
survival_dat <- clean_names(readRDS(glue("{data_path_bestageing2022}/data/202211908_XMELD_abfrage_best_ageing.rds"))) # %>% 
# original path "/mnt/users/reich/XMeldPortal_neu/meldeportal-tools-meldeportalclient-9.3/Rout/202211908_XMELD_abfrage_best_ageing.rds"

## metadata from DB
# https://www.bestageing.org/Pages/Login.aspx?ReturnUrl=%2f&AspxAutoDetectCookieSupport=1
load(glue("{data_path_BestAgeing}/data/clean_all_meta.rda"))  # created in "scripts/_prepare_metadata.R"
clean_all_meta <- clean_all_meta %>% 
  mutate(age = ifelse(age < 18, NA, age))  # wrong age remove
# cath data? "hkdb"

## load all original metadat xlsx files again to make sure that also overlapped 
#patients (e.g. dcm+cad) are in each group
control_ids <- read_excel(glue("{data_path_BestAgeing}/data/pheno_controls.xlsx")) %>% 
  dplyr::pull(BestAgeingCode)

# "UKL-HD-00318" both in Control and CAD dataset, looked it up (HK Nr 1289-2015): KHK ohne hg Stenosen, LV gut --> assign to CAD only
control_ids <- control_ids[control_ids != "UKL-HD-00318"]

###
# We apply both parametric t-tests and nonparametric U-tests ----------------

###
# START LOOP for specified combinations -------------------------------------
###

p2_list <- list()

for (i in 1:nrow(all_combis[all_combis[["analysis"]] == "full", ])) {  # only consider all miRNAs for DE (not research miRNAs only)
  ## reassign disease since only full analysis here
  disease <- all_combis[all_combis[["analysis"]] == "full", ]$diseases[i]
  
  ## DATA FIRST
  # create parameter specific {disease}_ids
  filename <- paste0(data_path_BestAgeing , "/data/pheno_", disease, ".xlsx")
  disease_vector <- paste0(disease, "_ids")
  # assign vector to string variable above    (CAVE: to evaluate -> `eval(as.symbol(X))`)
  assign(x = disease_vector,
         value = read_excel(filename) %>% 
           dplyr::pull(BestAgeingCode))
  # df with pat_id and disease identifier
  disease_ident_df <- tibble(pat_id = c(eval(as.symbol(disease_vector)),control_ids), 
                             disease = c(rep(disease, length(eval(as.symbol(disease_vector)))),
                                         rep("control", length(control_ids)))
  )
  
  # prepare-combine-data
  data01 <- all_mirnas %>% 
    # filter control patients and from {disease}_ids
    filter(pat_id %in% !!sym(disease_vector) |
             pat_id %in% control_ids
    ) %>% 
    # bind age and gender from metadata
    left_join(clean_all_meta %>% select(patID, age, sex), 
              by = c("pat_id"= "patID") ) %>% 
    left_join(disease_ident_df, by = "pat_id") %>% 
    relocate(c(disease, age, sex), .after = pat_id) %>% 
    mutate(disease = factor(disease),
           sex = factor(sex) ) %>% 
    as_tibble()
  
  # if (SAVE.files ==TRUE) {
  #   filename.data01 <- glue("{data_path_bestageing2022}/output/de_results/{disease}/data01.rds")
  #   saveRDS(object = data01, file = filename.data01)
  # }
  
  ## SUBSET of literature miRNAs
  # selected_mirna_dat <- paste0("data01_", length(researchMiRNAAccession$miRNAName_v21), "mirnas")  
  # assign(x = selected_mirna_dat, 
  #        value = data01 %>% 
  #          select(pat_id, disease, age, sex, researchMiRNAAccession$miRNAName_v21))
  
  ## run t-tests & log reg model  ----------------------------------------------------
  # create space to store results for ALL miRNAs
  pval.t.test<-rep(NA,ncol(all_mirnas)-1)
  pval.u.test<-rep(NA,ncol(all_mirnas)-1)
  pval.glm <- rep(NA,ncol(all_mirnas)-1)
  #average.difference <- rep(NA,ncol(all_mirnas)-1)
  log2FoldChange <- rep(NA,ncol(all_mirnas)-1)
  median.cont <- rep(NA,ncol(all_mirnas)-1)
  median.case <- rep(NA,ncol(all_mirnas)-1)
  mean.cont <- rep(NA,ncol(all_mirnas)-1)
  mean.case <- rep(NA,ncol(all_mirnas)-1)
  empse.cont <- rep(NA,ncol(all_mirnas)-1)
  empse.case <- rep(NA,ncol(all_mirnas)-1)
  aucs <- rep(NA,ncol(all_mirnas)-1)
  auc_glm <- rep(NA,ncol(all_mirnas)-1)
  name.mir <- rep(NA,ncol(all_mirnas)-1)
  # indexing
  cont.index <- data01$disease == "control"
  case.index <- data01$disease == disease
  # run tests
  
  if (runTests == TRUE){  #takes time
    for(miRNA in 1:(ncol(all_mirnas)-1)) {
      # update index since first colnames are [1] "pat_id"  "disease"  "age"   "sex"  "hsa_let_7a_3p" 
      miRNA_col <- miRNA+4
      cont <- data01[cont.index, miRNA_col] %>% as_vector()
      case <- data01[case.index, miRNA_col] %>% as_vector()
      # median for log-median expression plot (Figure 1 Vogel2013)
      mean.cont[miRNA] <- mean(cont)
      mean.case[miRNA] <- mean(case)
      median.cont[miRNA] <- median(cont)
      median.case[miRNA] <- median(case)
      empse.cont[miRNA] <- sd(cont)/ sqrt(length(cont))
      empse.case[miRNA] <- sd(case)/ sqrt(length(case))
      aucs[miRNA] <- suppressMessages(pROC::roc(controls=cont, cases=case)$auc[[1]])
      name.mir[miRNA] <- names(data01[miRNA_col])
      
      # average difference and logfold
      mean.control <- mean(cont)
      mean.case <- mean(case)
      log2FoldChange[miRNA] <- mean.case - mean.control
      ##log2FoldChange[miRNA] <- (mean.case - mean.control)/mean.control   # we are already on the log2 scale (https://support.bioconductor.org/p/117881/)
      # pvals
      pval.t.test[miRNA] <-t.test(cont,case)$p.value 
      pval.u.test[miRNA] <- wilcox.test(as.numeric(cont), as.numeric(case), exact = FALSE)$p.value
      # glm
      name_miRNA <- colnames(data01)[miRNA_col]
      f <- as.formula(paste("disease ~ ", name_miRNA, " + sex + age"))
      logreg <- glm(formula = f, data = data01, family = binomial(link = "logit") )
      predictions <- predict(logreg, newdata = data01[ , c("age", "sex", name_miRNA)], type = "response")
      auc_glm[miRNA] <- suppressMessages(roc(data01$disease, predictions)$auc[[1]])
      
      pval.glm[miRNA] <- coef(summary(logreg))[2,4]
    }
    
    # GATHER Results -----------------------------------------------------------
    # we conducted 2549 t-tests and 2549 glm-models for each gene
    de.results <- tibble(miRNA = colnames(all_mirnas)[-1],  # all miRNAs without patID
                         # average.difference = average.difference,
                         log2FoldChange = log2FoldChange,
                         pval.t.test = pval.t.test,
                         pval.u.test = pval.u.test,
                         pval.glm = pval.glm,
                         aucs_univariate = aucs,
                         aucs_glm = auc_glm)
    
    # for figure 1 vogel 2013 
    results_logmedians <- tibble(miR =name.mir, auc=aucs, aucs_glm = auc_glm, pval.t.test = pval.t.test, pval.glm = pval.glm,
                                 logmedian.cont = median.cont, logmedian.case = median.case,
                                 logmean.cont = mean.cont, logmean.case = mean.case, empse.case, empse.cont) %>% 
      mutate(auc= ifelse(auc<0.5, 1-auc, auc)) %>% 
      # changed on 2023-08-02
      mutate(padj = p.adjust(pval.t.test, method = "holm", n = length(name.mir)),  # inflation with "BH", use Bonferroni-Holm
             padj.glm = p.adjust(pval.glm, method = "holm", n = length(name.mir))) %>% 
      mutate(sign_indicator = ifelse(padj < 0.05, "p.adj≤0.05", "n.s."))
    
    if (SAVE.files ==TRUE) {
      filename.de.tibble <- glue("{data_path_bestageing2022}/output/de_results/{disease}/de_results.rds")
      saveRDS(object = de.results, file = filename.de.tibble)
      filename.logmedians <- glue("{data_path_bestageing2022}/output/de_results/{disease}/results_logmedians.rds")
      saveRDS(object = results_logmedians, file = filename.logmedians)
    }
  }
  
  rm(de.results, results_logmedians)
  print(paste0("|||-----------------------Run finished for disease: ", toupper(disease), " -----------------------|||"))
}



rm(disease)  # errors in plot otherwise because colname and varname

# changed from "rockerprojects/bestageing2022/data-literature/miRetrieve/top50mirnas_all_diseases.rds"
miRetrieve_alldiseases <- readRDS(glue("{data_path_bestageing2022}/data-literature/miRetrieve/2023-07-27_top50mirnas_all_diseases_pmids_gpt.rds")) # created in "scripts/miRetrieve/miRetrieve_topmirnas_all_diseases.R"
length(unique(miRetrieve_alldiseases$Accession))

p2_listnew <- list()
for(mydisease in 1:nrow(all_combis[all_combis[["analysis"]] == "full", ])){
  data01 <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/data01.rds"))
  de_results <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/de_results.rds"))
  results_logmedians <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/results_logmedians.rds"))
  
  # figure 01 vogel plot ------------------------------------------------------
  p1.1 <- ggplot(results_logmedians, aes(x = logmedian.cont, y = logmedian.case)) +
    geom_point(alpha=0.3, shape=16) +
    #geom_count(color="black", size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
    labs(x = glue("controls\n[log-median expression]"), y = glue("{toupper(all_combis$diseases[mydisease])} patients\n[log-median expression]")) +
    #theme_classic()
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    my_base_theme()
  
  p1.2 <- ggplot() +
    geom_histogram(data = results_logmedians, aes(x = padj), fill = "white", color = "black", # bins = 500,
                   breaks = seq(0, 1, by = 0.05) ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.05), expand = c(0, 0)) +
    geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", size = 1) +
    labs(x = "t-test P-value\n(Bonferroni-Holm adjusted)", y = "Frequency") +
    #theme_classic()
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    my_base_theme()
  
  p1.3 <- ggplot() +
    geom_histogram(data = results_logmedians, aes(x = auc), fill = "white", color = "black", binwidth = 0.01) +
    scale_x_continuous(breaks = seq(0.5, 0.7, by = 0.1), limits = c(0.45, 0.8), expand = c(0, 0)) +
    geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", size = 1) +
    labs(x = "Univariate AUC", y = "Frequency") +
    #theme_classic()
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    my_base_theme()
  
  baseline_glm_model <- glm(formula = disease~age+sex, data = data01, family = binomial(link = "logit") )
  predictions <- predict(baseline_glm_model, newdata = data01[ , c("age", "sex")], type = "response")
  auc_glm_baseline <- suppressMessages(roc(data01$disease, predictions)$auc[[1]])
  p1.3b <- ggplot() +  # adjusted glm for age and sex
    geom_histogram(data = results_logmedians, aes(x = aucs_glm), fill = "white", color = "black", binwidth = 0.01) +
    scale_x_continuous(breaks = seq(0.6, 0.9, by = 0.1), limits = c(0.55, 0.85), expand = c(0, 0)) +
    geom_vline(xintercept = auc_glm_baseline, color = "red", linetype = "dashed", size = 1) +
    labs(x = "Univariate AUC", y = "Frequency") +
    #theme_classic()
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    my_base_theme()
  
  # Arrange the plots into two rows
  plot_row1 <- p1.1 + labs(title = "A")
  plot_row2 <- ggarrange(p1.2 + labs(title = "B"), p1.3 + labs(title = "C"), nrow = 1, ncol = 2, widths = c(1, 1))
  plot_row2b <- ggarrange(p1.2 + labs(title = "B"), p1.3b + labs(title = "C"), nrow = 1, ncol = 2, widths = c(1, 1))
  # Combine the plots into a single plot with two rows
  fig1vogel2013 <- ggarrange(plot_row1, plot_row2, nrow = 2, heights = c(2, 1))
  fig1vogel2013b <- ggarrange(plot_row1, plot_row2b, nrow = 2, heights = c(2, 1))
  print(fig1vogel2013)
  
  if (SAVE.files ==TRUE) {
    fig1vogel2013_path <- glue("{data_path_bestageing2022}/output/plots/fig01vogel2013/{Sys.Date()}_{toupper(all_combis$diseases[mydisease])}_linear_correlation.svg")
    fig1vogel2013_path_b <- glue("{data_path_bestageing2022}/output/plots/fig01vogel2013/{Sys.Date()}_{toupper(all_combis$diseases[mydisease])}_linear_correlation_b.svg")
    ggsave(filename = fig1vogel2013_path, plot = fig1vogel2013, 
           width = 8, height = 8, 
           units = "in"  # default
    )
    ggsave(filename = fig1vogel2013_path_b, plot = fig1vogel2013b, 
           width = 8, height = 8, 
           units = "in"  # default
    )
  }
  
  ## figure 02 vogel 2013. ----------------------------
  # Step 1: select miRNAs from miRetrieve and disease
  disease_miRetrieve_mirnas <- miRetrieve_alldiseases %>% 
    mutate(Topic = toupper(Topic)) %>% 
    filter(Topic == toupper(all_combis$diseases[mydisease])) %>% 
    pull(TargetName)
  no_mirnas <- 10
  topNmiRetrieve <- results_logmedians %>% 
    filter(miR %in% disease_miRetrieve_mirnas) %>% 
    arrange(padj, desc(auc)) %>% 
    slice(1:no_mirnas)
  
  # long format
  rainbow_plot_dat <- data01 %>% 
    select(disease, all_of(topNmiRetrieve$miR)) %>% 
    pivot_longer(cols=-disease, names_to = "feature", values_to = "value") 

  # Calculate maximum values for each feature
  df_max <- rainbow_plot_dat %>% 
    group_by(feature, disease) %>% 
    summarise(max_value = max(value, na.rm = TRUE)) %>% 
    arrange(desc(max_value)) %>% 
    ungroup() %>% 
    select(-disease) %>% 
    distinct(feature, .keep_all = TRUE) %>% 
    # label
    left_join(topNmiRetrieve %>% select(miR, sign_indicator, padj), by=c("feature"="miR"))  %>% 
    mutate(feature = gsub("_", "-", feature))
  
  df_max <- df_max %>% 
    mutate(sign_indicator_asterisks = case_when(
      padj < 0.001 ~ "***",
      padj>= 0.001 & padj < 0.01 ~ "**",
      padj>= 0.01 & padj < 0.05 ~ "*",
      .default = "n.s."
    ))
  
  rainbow_plot_dat <- rainbow_plot_dat %>% 
    mutate(feature = gsub("_", "-", feature))
  
  # Order the levels of the "feature" variable based on the median values
  feature_order <- rainbow_plot_dat %>%
    group_by(feature) %>%
    summarise(median_value = median(value)) %>%
    arrange(desc(median_value)) %>%
    pull(feature)
  rainbow_plot_dat$feature <- factor(rainbow_plot_dat$feature, levels = feature_order)
  rainbow_plot_dat$disease <- factor(rainbow_plot_dat$disease, levels = c(all_combis$diseases[mydisease], "control"))
  
  # grouped box plots
  p2 <- ggplot(rainbow_plot_dat, aes(x=feature, y=value, fill = disease)) +
    geom_boxplot(position=position_dodge(width=0.8), outlier.shape = NA, alpha=0.7)+    
    geom_jitter(aes(color=disease),position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), size=0.1, alpha=0.2) +
    # prevent geom_text() from searching for the disease variable, you can override the fill aesthetic by setting fill = NULL within the aes() function 
    geom_text(data = df_max, aes(y = max_value, x=feature, fill = NULL), label = df_max$sign_indicator_asterisks, vjust = -0.5, size=2)+  
    labs(x = "", y = expression(paste("log"[2], " expression")), title = "") +
    scale_y_continuous(breaks = seq(0, ceiling(max(rainbow_plot_dat$value)), by = 2), #limits = c(0, 1), expand = c(0, 0)
    )+
    # theme_classic()+
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6), labels = c(toupper(all_combis$diseases[mydisease]), "Control")) + guides(fill=guide_legend(title=NULL), color="none") +
    scale_color_manual(values = thematic::okabe_ito(6)) +
    #scale_color_brewer(palette = "Set1") + # Choose a color palette
    #scale_fill_brewer(palette = "Set1", labels = c(toupper(all_combis$diseases[mydisease]), "Control")) + guides(fill=guide_legend(title=NULL), color="none")+
    my_base_theme()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p2)
  
  p2_listnew[[mydisease]] <- p2
  if (SAVE.files ==TRUE) {
    fig2_disease_vogel2013_path <- glue("{data_path_bestageing2022}/output/plots/fig02vogel2013/{Sys.Date()}_{toupper(all_combis$diseases[mydisease])}_topdysregulatedmiRetrievemiR.svg")
    ggsave(filename = fig2_disease_vogel2013_path, plot = p2, 
           width = 8, height = 8, 
           units = "in"  # default
    )
  }
  #combine fig01 and fig02 for disease
  #plot_0102_combined <- ggarrange(fig1vogel2013b, p2 + labs(title = "D"), nrow = 2, ncol = 1, heights = c(1.5, 1))
  #plot_0102_combined
}

fig2vogel2013 <- ggarrange(plotlist = p2_listnew, ncol = 2, nrow = 2, 
                           labels = c("A", "B", "C", "D"))

# To print the plot
print(fig2vogel2013)
if (SAVE.files ==TRUE) {
  fig2vogel2013_path <- glue("{data_path_bestageing2022}/output/plots/fig02vogel2013/{Sys.Date()}_arranged_topdysregulatedmiRetrievemiR.svg")
  ggsave(filename = fig2vogel2013_path, plot = fig2vogel2013, 
         width = 12, height = 12, 
         units = "in"  # default
  )
}



# DE volcanos ----------------------------------------------------------------
# moved here from "main.Rmd"
df_de_mirna_topic <- tibble(
  miRNA = character(0),
  Topic = character(0)
)

for(mydisease in 1:nrow(all_combis[all_combis[["analysis"]] == "full", ])) {
  # de.results_old <- readRDS(dplyr::last(list.files(path = "/Users/christophreich/Desktop/mount/rockerprojects/bestageing2022/output/de_results", pattern = paste0("_de_results_DISEASE_", disease_vector[disease]), full.names = TRUE)))
  data01 <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/data01.rds"))
  de.results <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/de_results.rds"))
  results_logmedians <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/results_logmedians.rds"))
  
  ### FDR calculation
  de.results$padj <- p.adjust(de.results$pval.t.test, method = "holm", n = length(de.results$pval.t.test))  # "BH" p-val-inflation
  de.results$padj.glm <- p.adjust(de.results$pval.glm, method = "holm", n = length(de.results$pval.t.test))
  ## manual calc 
  # n.comparisons <- length(de.results$pval.t.test)
  # inv.rank<-n.comparisons-rank(de.results$pval.t.test)
  # adjusted.results.t.test<-de.results$pval.t.test*(n.comparisons)/(n.comparisons-inv.rank+1)
  
  de.results$diffexpressed <- "NO"
  # if log2Foldchange > 1.0 and pvalue < 0.05, set as "UP" 
  de.results$diffexpressed[de.results$log2FoldChange > 1.0 & 
                             de.results$padj < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  de.results$diffexpressed[de.results$log2FoldChange < -1.0 & 
                             de.results$padj < 0.05] <- "DOWN"
  
  log2folds_thresholds <- c(1, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
  for (i in 1:length(log2folds_thresholds)) {
    if ( sum(de.results$diffexpressed != "NO") < 15 ) {
      de.results$diffexpressed[de.results$log2FoldChange > log2folds_thresholds[i] & 
                                 de.results$padj < 0.05] <- "UP"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      de.results$diffexpressed[de.results$log2FoldChange < - log2folds_thresholds[i] & 
                                 de.results$padj < 0.05] <- "DOWN"
      log2folds_thresh <- log2folds_thresholds[i]
    } else{ 
      log2folds_thresh <- log2folds_thresholds[i]
      break
    }
  }
  
  ## label miRNAs 
  # Step 1: select miRNAs from miRetrieve and disease
  disease_miRetrieve_mirnas <- miRetrieve_alldiseases %>% 
    mutate(Topic = toupper(Topic)) %>% 
    #filter(Topic == toupper(all_combis$diseases[mydisease])) %>% 
    pull(TargetName)
  
  # init
  de.results$delabel <- NA
  de.results$delabel[de.results$diffexpressed != "NO"] <- 
    de.results$miRNA[de.results$diffexpressed != "NO"]
  # ONLY label miRNAs that are also known from LITERATURE 
  ## old literature miRNAs 
  # de.results$delabel <- ifelse(de.results$delabel %in% researchMiRNAAccession$miRNAName_v21,  de.results$delabel, NA)
  # NEW from miRetrieve
  de.results$delabel <- ifelse(de.results$delabel %in% disease_miRetrieve_mirnas,  de.results$delabel, NA)
  
  ## color all literature miRNAs ----------------------------------------
  #de.results$research_mirna <- ifelse(de.results$miRNA %in% researchMiRNAAccession$miRNAName_v21,  "Literature miRNA", "Not found in Literature")
  de.results$research_mirna <- ifelse(de.results$miRNA %in% disease_miRetrieve_mirnas,  "Literature miRNA", "Not found in miRetrieve Top Hits")
  
  ## # label miRNA 106b-5p (special interest)
  ## freys.mirnas <- c("hsa_mir_106b_3p", "hsa_mir_106b_5p")
  ## de.results$freysmirna <- ifelse(de.results$miRNA %in% freys.mirnas, de.results$miRNA, NA)
  ## de.results$freysmirna.col <- ifelse(de.results$miRNA %in% freys.mirnas, de.results$miRNA, "not selected")
  
  ## PLOT ##
  # custom volcano -----------------------------------------------------------
  subtitle.custom <- paste0("log2foldchange = \u00b1 ", log2folds_thresh, " for significant miRNAs (p_adj < 0.05)")
  
  de.results %>% 
    mutate(delabel = gsub("_", "-", delabel)) %>% 
    ggplot(aes(x=log2FoldChange, y=-log10(pval.t.test), col=diffexpressed, label=delabel)) + 
    geom_point(alpha=ifelse(de.results$diffexpressed != "NO", 1, 0.2), size=1, shape=16) + 
    geom_vline(xintercept=c(-log2folds_thresh, log2folds_thresh), alpha=0.3) +
    geom_hline(
      yintercept=-log10(0.05/length(pval.t.test)),  # bonferroni only correction
      alpha=0.3
    ) + 
    labs(
      title = NULL, # paste0("Volcano plot of DE miRNAs in ", disease_vector[mydisease]),
      subtitle = NULL, #subtitle.custom
    ) +
    xlab("log2 fold change") + 
    ylab("-log10 (p-value)") + 
    theme(legend.position = "none", 
          plot.title = element_text(size = rel(1.5), hjust = 0.5), 
          axis.title = element_text(size = rel(1.25)))+
    #ggthemes::scale_color_few(name="Differentially expressed")+
    ggrepel::geom_text_repel(show.legend = FALSE, size=5)+
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    scale_color_manual(values = thematic::okabe_ito(6), name="Differentially expressed")+
    my_base_theme() -> volcano_1
  print(volcano_1)
  
  if (SAVE.files ==TRUE ){
    filename.volcano <- glue("{data_path_bestageing2022}/output/plots/de_analysis/{all_combis$diseases[mydisease]}/{Sys.Date()}_de_volcano.svg")
    ggsave(filename = filename.volcano, plot = volcano_1, 
           width = 9, height = 5, 
           units = "in"  # default
    )
  }
  
  # Color all literature miRNAs
  de.results %>% 
    mutate(delabel = gsub("_", "-", delabel)) %>% 
  ggplot(aes(x=log2FoldChange, y=-log10(pval.t.test), col=research_mirna, label=delabel)) + 
    geom_point(alpha=ifelse(de.results$research_mirna == "Literature miRNA", 1, 0.2), size=1, shape=16) + 
    geom_vline(xintercept=c(-log2folds_thresh, log2folds_thresh), alpha=0.3) +
    geom_hline(
      yintercept=-log10(0.05/length(pval.t.test)),  # bonferroni only correction
      alpha=0.3
    ) + 
    labs(
      title = NULL, # paste0("Volcano plot of DE miRNAs in ", disease_vector[mydisease]),
      subtitle = NULL, # subtitle.custom
    ) +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") + 
    theme(legend.position = "none", 
          plot.title = element_text(size = rel(1.5), hjust = 0.5), 
          axis.title = element_text(size = rel(1.25)))+
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    scale_color_manual(values = thematic::okabe_ito(6), name=NULL)+
    ggrepel::geom_text_repel(show.legend = FALSE, size=5) -> volcano_2
  print(volcano_2)
  
  if (SAVE.files ==TRUE ){
    filename.volcano2 <-glue("{data_path_bestageing2022}/output/plots/de_analysis/{all_combis$diseases[mydisease]}/{Sys.Date()}_de_volcano_LITERATmiRNAs_ANNOT.svg")
    ggsave(filename = filename.volcano2, plot = volcano_2, 
           width = 9, height = 5, 
           units = "in"  # default
    )
  }
  
  # append list for venn diagram (intersection of de mirnas in different diseases)
  df_de_mirna_topic_tmp <- de.results %>% 
    filter(padj.glm < 0.05) %>% 
    select(miRNA) %>% 
    mutate(Topic = all_combis$diseases[mydisease])
  df_de_mirna_topic <- df_de_mirna_topic %>% 
    bind_rows(df_de_mirna_topic_tmp)
}

# Venn diagramm 4 miRetrieve -------------------------------------------------------------

list_vann <- list(
  ACS = miRetrieve_alldiseases$TargetName[miRetrieve_alldiseases$Topic == "ACS"],
  CAD = miRetrieve_alldiseases$TargetName[miRetrieve_alldiseases$Topic == "CAD"],
  DCM = miRetrieve_alldiseases$TargetName[miRetrieve_alldiseases$Topic == "DCM"],
  HFrEF = miRetrieve_alldiseases$TargetName[miRetrieve_alldiseases$Topic == "HFrEF"]
)

# Create the Venn diagram
venn_plot <- ggvenn(list_vann) +
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_fill_manual(values = thematic::okabe_ito(6)) +
  scale_color_manual(values = thematic::okabe_ito(6)) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )
venn_plot_high <- venn_plot + coord_fixed(ratio = 2) 
venn_plot

if (SAVE.files ==TRUE ){
  filename.venn_miRetrieve <-glue("{data_path_bestageing2022}/output/plots/venn_dia/venn_LITERATmiRNAs_ANNOT.svg")
  ggsave(filename = filename.venn_miRetrieve, plot = venn_plot, 
         width = 8, height = 8, 
         units = "in"  # default
  )
  filename.venn_miRetrieve_dimension_high <-glue("{data_path_bestageing2022}/output/plots/venn_dia/venn_LITERATmiRNAs_ANNOT_high.svg")
  ggsave(filename = filename.venn_miRetrieve_dimension_high, plot = venn_plot_high, 
         width = 8, height = 8, 
         units = "in"  # default
  )
  saveRDS(object = venn_plot_high, file = glue("{data_path_bestageing2022}/output/plots/venn_dia/venn_de_mirnas_shared.rds"))
  
}

# Extracting shared names
shared_elements_all_diseases <- Reduce(intersect, list_vann)

# Displaying the shared names in a sophisticated way
if (length(shared_elements) > 0) {
  cat("Shared elements between the sets are:\n")
  cat(paste(shared_elements, collapse = ", "), "\n")
} else {
  cat("No shared elements between the sets.\n")
}

# Calculating intersections
intersections <- gplots::venn(list_vann)
list_of_intersections <- attr(intersections, "intersections")

# Displaying the intersections
intersections_diseases <- tibble(
  Intersection = character(),
  miRNAs = character(), 
  count_intersections = integer()
)
for (name in seq_along(list_of_intersections)) {
  intersection_name <- names(list_of_intersections[name])
  miRNAs_joined <- paste(gsub("_", "-", list_of_intersections[[name]]), collapse = ", ")
  no_mirnas <- length(list_of_intersections[[name]])
  
  cat(intersection_name, ":\n", miRNAs_joined, "\n\n")
  
  intersections_diseases <- intersections_diseases %>% 
    add_row(Intersection = intersection_name, miRNAs = miRNAs_joined, count_intersections = no_mirnas)
}

intersections_diseases <- intersections_diseases %>% arrange(desc(count_intersections))
write.csv2(x = intersections_diseases, file = glue("{data_path_bestageing2022}/output/tables/venn_dia/top50perdisease_venn.csv"))

length(unique(miRetrieve_alldiseases$TargetName))

# Venn plot DE-mirnas shared ----------------------------------------------

list_venn_de_mirnas <- list(
  ACS = df_de_mirna_topic$miRNA[df_de_mirna_topic$Topic == "acs"],
  CAD = df_de_mirna_topic$miRNA[df_de_mirna_topic$Topic == "cad"],
  DCM = df_de_mirna_topic$miRNA[df_de_mirna_topic$Topic == "dcm"],
  HFrEF = df_de_mirna_topic$miRNA[df_de_mirna_topic$Topic == "hfref"]
)

# Create the Venn diagram
venn_plot_de_mirnas <- ggvenn(list_venn_de_mirnas) +
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_fill_manual(values = thematic::okabe_ito(6)) +
  scale_color_manual(values = thematic::okabe_ito(6)) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )
venn_plot_de_mirnas

if (SAVE.files ==TRUE ){
  filename.venn_de_mirnas_shared <-glue("{data_path_bestageing2022}/output/plots/venn_dia/venn_de_mirnas_shared.svg")
  ggsave(filename = filename.venn_de_mirnas_shared, plot = venn_plot_de_mirnas, 
         width = 8, height = 8, 
         units = "in"  # default
  )
}


# Check p-values ------------------------------------------------------------

for(mydisease in 1:nrow(all_combis[all_combis[["analysis"]] == "full", ])) {
  # de.results_old <- readRDS(dplyr::last(list.files(path = "/Users/christophreich/Desktop/mount/rockerprojects/bestageing2022/output/de_results", pattern = paste0("_de_results_DISEASE_", disease_vector[disease]), full.names = TRUE)))
  data01 <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/data01.rds"))
  de.results <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/de_results.rds"))
  results_logmedians <- readRDS(file = glue("{data_path_bestageing2022}/output/de_results/{all_combis$diseases[mydisease]}/results_logmedians.rds"))
  
  ### FDR calculation
  de.results$padj <- p.adjust(de.results$pval.t.test, method = "holm", n = length(de.results$pval.t.test))  # "BH" p-val-inflation
  de.results$padj.glm <- p.adjust(de.results$pval.glm, method = "holm", n = length(de.results$pval.t.test))
  
  # Create a data frame for expected and observed p-values
  qq_df <- tibble(
    Expected = -log10((1:length(de.results$pval.t.test)) / length(de.results$pval.t.test)),
    Observed_t = -log10(sort(de.results$pval.t.test)),
    Observed_t_adj = -log10(sort(de.results$padj)),
    Observed_glm = -log10(sort(de.results$pval.glm)),
    Observed_glm_adj = -log10(sort(de.results$padj.glm)),
    )
  
  # Create the QQ-plot
  ggplot(qq_df, aes(Expected, Observed_t)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    xlab("Expected -log10 p-values") +
    ylab("Observed -log10 p-values")
  ggplot(qq_df, aes(Expected, Observed_t_adj)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    xlab("Expected -log10 p-values") +
    ylab("Observed -log10 p-values")
  ggplot(qq_df, aes(Expected, Observed_glm)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    xlab("Expected -log10 p-values") +
    ylab("Observed -log10 p-values")
  ggplot(qq_df, aes(Expected, Observed_glm_adj)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    xlab("Expected -log10 p-values") +
    ylab("Observed -log10 p-values")
  
}

# we experienced p-val inflation

# 1) Check experimental design and preprocessing steps

# Perform PCA -----------------------------
# prep dat
exprs_metadat <- all_mirnas %>% 
  # bind age and gender from metadata
  left_join(clean_all_meta %>% select(patID, disease, age, sex), 
            by = c("pat_id"= "patID") ) %>% 
  relocate(c(disease, age, sex), .after = pat_id) %>% 
  filter(disease != "pef") %>% 
  mutate(disease = factor(disease),
         sex = factor(sex) )

exprs <- exprs_metadat %>% select(-c(pat_id, disease, age, sex))
pca_res <- prcomp(exprs, center = TRUE, scale. = TRUE)

# Get the first two principal components
pca_df <- as_tibble(as.data.frame(pca_res$x[,1:5]))

# add metadata
exprs_metadat <- exprs_metadat %>% 
  bind_cols(pca_df)

# PC-plot
# Custom labels for the legend
custom_labels <- c(
  "acs" = "ACS",
  "cad" = "CAD",
  "control" = "Control",
  "dcm" = "DCM",
  "ref" = "HFrEF"
)

ggplot(exprs_metadat, aes(PC1, PC2, color = disease)) +
  geom_point(alpha=0.7, shape=16) +
  labs(title = NULL,  # "PCA Plot Colored by Disease",
       x = "PC1",
       y = "PC2") +
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_fill_manual(values = thematic::okabe_ito(6)) +
  scale_color_manual(values = thematic::okabe_ito(6), name=NULL, labels = custom_labels)+
  my_base_theme() -> plot_PC1_PC2
plot_PC1_PC2

ggplot(exprs_metadat, aes(PC1, PC3, color = disease)) +
  geom_point(alpha=0.7, shape=16) +
  labs(title = NULL,  # "PCA Plot Colored by Disease",
       x = "PC1",
       y = "PC3") +
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_fill_manual(values = thematic::okabe_ito(6)) +
  scale_color_manual(values = thematic::okabe_ito(6), name=NULL, labels = custom_labels)+
  my_base_theme() -> plot_PC1_PC3

ggplot(exprs_metadat, aes(PC1, PC4, color = disease)) +
  geom_point(alpha=0.7, shape=16) +
  labs(title = NULL,  # "PCA Plot Colored by Disease",
       x = "PC1",
       y = "PC4") +
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_fill_manual(values = thematic::okabe_ito(6)) +
  scale_color_manual(values = thematic::okabe_ito(6), name=NULL, labels = custom_labels)+
  my_base_theme() -> plot_PC1_PC4

ggplot(exprs_metadat, aes(PC1, PC5, color = disease)) +
  geom_point(alpha=0.7, shape=16) +
  labs(title = NULL,  # "PCA Plot Colored by Disease",
       x = "PC1",
       y = "PC5") +
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_fill_manual(values = thematic::okabe_ito(6)) +
  scale_color_manual(values = thematic::okabe_ito(6), name=NULL, labels = custom_labels)+
  my_base_theme() -> plot_PC1_PC5

# add marginal density plots
# Create the marginal density plots
x_density <- ggplot(exprs_metadat, aes(PC1, fill = disease)) +
  geom_density(alpha = 0.5) +
  guides(fill = "none") +
  xlab("") +
  ylab("")+
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_fill_manual(values = thematic::okabe_ito(6)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

y_density <- ggplot(exprs_metadat, aes(PC3, fill = disease)) +
  geom_density(alpha = 0.5) +
  coord_flip() +
  xlab("") +
  ylab("")+
  guides(fill = "none") +
  theme_minimal(base_size = 16, base_family = 'Arial')+
  scale_fill_manual(values = thematic::okabe_ito(6)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

# Combine the plots
grid.arrange(x_density, 
             NULL, 
             plot_PC1_PC2, 
             y_density, 
             ncol = 2, 
             nrow = 2, 
             widths = c(4, 1), 
             heights = c(1, 4))

# better plot alignment
combined_plot <- plot_grid(x_density, NULL, plot_PC1_PC2, y_density,
                           ncol = 2, nrow = 2, align = 'hv', axis = 'tblr',
                           rel_widths = c(4, 1.5), rel_heights = c(1, 4)
                           )
combined_plot
