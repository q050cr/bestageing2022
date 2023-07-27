
### INFO ----------------------------------------------------------------------
# Differential miRNA Expression Analysis
# this script is sourced from `scripts/render_param_reports.R`
# selection provided by `all_combis$diseases` and `all_combis$analysis`

# script creates plots: "fig01vogel2013", "fig02vogel2013"

.libPaths("/mnt/users/reich/programs/R43/lib")
# dependencies ---------------------------------------------------------------
library(readxl, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(janitor, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(glue, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(gt, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(dplyr, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(tidyr, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(stringr, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(purrr, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(dplyr, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(ggplot2, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(RColorBrewer, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(ggdist, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(gghalves, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(ggrepel, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(rstatix, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(ggthemes, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(ggpubr, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(pROC, lib.loc = "/mnt/users/reich/programs/R43/lib")
library(reshape2, lib.loc = "/mnt/users/reich/programs/R43/lib")
conflicted::conflict_prefer("expand", "tidyr")
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(janitor::make_clean_names)

source(file = "/mnt/users/reich/rockerprojects/bestageing2022/scripts/helper/custom_ggplot_theme.R")

SAVE.files <- TRUE
runTests <- FALSE

diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) 


diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>% 
  filter(analysis=="full")


# load data ---------------------------------------------------------------

# MIRNA DAT
model_data1 <- clean_names(readRDS(file = '/mnt/users/reich/BestAgeing/data_new/model_data1.RDS'))  # has also multiclass col + diagnoses

load(file = "/mnt/users/reich/BestAgeing/data/mirnas.rda")  # "UKL-HD" n=765
load(file = "/mnt/users/reich/BestAgeing/data/data.rda")  # "UKL-HD" n=731
all_mirnas <- clean_names(mirnas)
names(all_mirnas) <- str_replace(string = names(all_mirnas), pattern = "mi_r", replacement = "mir")
mirnas_disease <- clean_names(data)
names(mirnas_disease) <- str_replace(string = names(mirnas_disease), pattern = "mi_r", replacement = "mir")
rm(data, mirnas)

# load-mirnas-from_research
# create vector of described mirnas
load("/mnt/users/reich/BestAgeing/data_research/fromR/researchMiRNAAccession.rda")
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
load(file = "/mnt/users/reich/BestAgeing/data/diagnoses_df.rda")

## SURVIVAL DAT
survival_dat <- clean_names(readRDS(file = '/mnt/users/reich/rockerprojects/bestageing2022/data/202211908_XMELD_abfrage_best_ageing.rds'))# %>% 
# original path "/mnt/users/reich/XMeldPortal_neu/meldeportal-tools-meldeportalclient-9.3/Rout/202211908_XMELD_abfrage_best_ageing.rds"

## metadata from DB
# https://www.bestageing.org/Pages/Login.aspx?ReturnUrl=%2f&AspxAutoDetectCookieSupport=1
load(file = "/mnt/users/reich/BestAgeing/data/clean_all_meta.rda")  # created in "scripts/_prepare_metadata.R"
clean_all_meta <- clean_all_meta %>% 
  mutate(age = ifelse(age < 18, NA, age))  # wrong age remove
# cath data? "hkdb"

## load all original metadat xlsx files again to make sure that also overlapped 
#patients (e.g. dcm+cad) are in each group
control_ids <- read_excel("/mnt/users/reich/BestAgeing/data/pheno_controls.xlsx") %>% 
  dplyr::pull(BestAgeingCode)

# "UKL-HD-00318" both in Control and CAD dataset, looked it up (HK Nr 1289-2015): KHK ohne hg Stenosen, LV gut --> assign to CAD only
control_ids <- control_ids[control_ids != "UKL-HD-00318"]

miRetrieve_alldiseases <- readRDS(file = "/mnt/users/reich/rockerprojects/bestageing2022/data-literature/miRetrieve/top50mirnas_all_diseases.rds") # created in "scripts/miRetrieve/miRetrieve_topmirnas_all_diseases.R"
length(unique(miRetrieve_alldiseases$Accession))

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
  filename <- paste0("/mnt/users/reich/BestAgeing/data/pheno_", disease, ".xlsx")
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
  
  if (SAVE.files ==TRUE) {
    filename.data01 <- glue("/mnt/users/reich/rockerprojects/bestageing2022/output/de_results/{disease}/data01.rds")
    saveRDS(object = data01, file = filename.data01)
  }
  
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
      name.mir[miRNA] <- names(data01[, miRNA_col])
      
      # average difference and logfold
      mean.control <- mean(cont)
      mean.case <- mean(case)
      log2FoldChange[miRNA] <- mean.case - mean.control
      ##log2FoldChange[miRNA] <- (mean.case - mean.control)/mean.control   # we are already on the log2 scale (https://support.bioconductor.org/p/117881/)
      # pvals
      pval.t.test[miRNA]<-t.test(cont,case)$p.value 
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
      mutate(padj = p.adjust(pval.t.test, method = "BH", n = length(name.mir)),
             padj.glm = p.adjust(pval.glm, method = "BH", n = length(name.mir))) %>% 
      mutate(sign_indicator = ifelse(padj < 0.05, "p.adj≤0.05", "n.s."))
    
    if (SAVE.files ==TRUE) {
      filename.de.tibble <- glue("/mnt/users/reich/rockerprojects/bestageing2022/output/de_results/{disease}/de_results.rds")
      saveRDS(object = de.results, file = filename.de.tibble)
      filename.logmedians <- glue("/mnt/users/reich/rockerprojects/bestageing2022/output/de_results/{disease}/results_logmedians.rds")
      saveRDS(object = results_logmedians, file = filename.logmedians)
    }
  }
  
  rm(de.results, results_logmedians)
  print(paste0("|||-----------------------Run finished for disease: ", toupper(disease), " -----------------------|||"))
}



rm(disease)  # errors in plot otherwise because colname and varname

p2_listnew <- list()
for(mydisease in 1:nrow(all_combis[all_combis[["analysis"]] == "full", ])){
  data01 <- readRDS(file = glue("/mnt/users/reich/rockerprojects/bestageing2022/output/de_results/{all_combis$diseases[mydisease]}/data01.rds"))
  de_results <- readRDS(file = glue("/mnt/users/reich/rockerprojects/bestageing2022/output/de_results/{all_combis$diseases[mydisease]}/de_results.rds"))
  results_logmedians <- readRDS(file = glue("/mnt/users/reich/rockerprojects/bestageing2022/output/de_results/{all_combis$diseases[mydisease]}/results_logmedians.rds"))
  
  # figure 01 vogel plot ------------------------------------------------------
  p1.1 <- ggplot(results_logmedians, aes(x = logmedian.cont, y = logmedian.case)) +
    geom_point(alpha=0.3, shape=16) +
    #geom_count(color="black", size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
    labs(x = glue("controls\n[log-median expression]"), y = glue("{toupper(all_combis$diseases[mydisease])} patients\n[log-median expression]")) +
    #theme_classic()
    theme_minimal(base_size = 16, base_family = 'Source Sans Pro')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    my_base_theme()
  
  p1.2 <- ggplot() +
    geom_histogram(data = results_logmedians, aes(x = padj), fill = "white", color = "black", #binwidth = 0.05) +
                   breaks = seq(0, 1, by = 0.05) ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.05), expand = c(0, 0)) +
    geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", size = 1) +
    labs(x = "t-test P-value\n(Benjamini-Hochberg corrected)", y = "Frequency") +
    #theme_classic()
    theme_minimal(base_size = 16, base_family = 'Source Sans Pro')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    my_base_theme()
  
  p1.3 <- ggplot() +
    geom_histogram(data = results_logmedians, aes(x = auc), fill = "white", color = "black", binwidth = 0.01) +
    scale_x_continuous(breaks = seq(0.5, 0.7, by = 0.1), limits = c(0.45, 0.8), expand = c(0, 0)) +
    geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", size = 1) +
    labs(x = "Univariate AUC", y = "Frequency") +
    #theme_classic()
    theme_minimal(base_size = 16, base_family = 'Source Sans Pro')+
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
    theme_minimal(base_size = 16, base_family = 'Source Sans Pro')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    my_base_theme()
  
  # Arrange the plots into two rows
  plot_row1 <- p1.1 + labs(title = "A")
  plot_row2 <- ggarrange(p1.2 + labs(title = "B"), p1.3 + labs(title = "C"), nrow = 1, ncol = 2, widths = c(1, 1))
  plot_row2b <- ggarrange(p1.2 + labs(title = "B"), p1.3b + labs(title = "C"), nrow = 1, ncol = 2, widths = c(1, 1))
  # Combine the plots into a single plot with two rows
  fig1vogel2013 <- ggarrange(plot_row1, plot_row2, nrow = 2, heights = c(2, 1))
  fig1vogel2013b <- ggarrange(plot_row1, plot_row2b, nrow = 2, heights = c(2, 1))
  
  if (SAVE.files ==TRUE) {
    fig1vogel2013_path <- glue("/mnt/users/reich/rockerprojects/bestageing2022/output/plots/fig01vogel2013/{toupper(all_combis$diseases[mydisease])}_linear_correlation.svg")
    fig1vogel2013_path_b <- glue("/mnt/users/reich/rockerprojects/bestageing2022/output/plots/fig01vogel2013/{toupper(all_combis$diseases[mydisease])}_linear_correlation_b.svg")
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
    arrange(desc(auc)) %>% 
    slice(1:no_mirnas)
  
  # long format
  rainbow_plot_dat <- data01 %>% 
    select(disease, all_of(topNmiRetrieve$miR)) %>% 
    pivot_longer(cols=-disease, names_to = "feature", values_to = "value")
  # Order the levels of the "feature" variable based on the median values
  feature_order <- rainbow_plot_dat %>%
    group_by(feature) %>%
    summarise(median_value = median(value)) %>%
    arrange(desc(median_value)) %>%
    pull(feature)
  rainbow_plot_dat$feature <- factor(rainbow_plot_dat$feature, levels = feature_order)
  rainbow_plot_dat$disease <- factor(rainbow_plot_dat$disease, levels = c(all_combis$diseases[mydisease], "control"))
  # Calculate maximum values for each feature
  df_max <- rainbow_plot_dat %>% 
    group_by(feature, disease) %>% 
    summarise(max_value = max(value, na.rm = TRUE)) %>% 
    arrange(desc(max_value)) %>% 
    ungroup() %>% 
    select(-disease) %>% 
    distinct(feature, .keep_all = TRUE) %>% 
    # label
    left_join(topNmiRetrieve %>% select(miR, sign_indicator, padj), by=c("feature"="miR"))  
  
  df_max <- df_max %>% 
    mutate(sign_indicator_asterisks = case_when(
      padj < 0.001 ~ "***",
      padj>= 0.001 & padj < 0.01 ~ "**",
      padj>= 0.01 & padj < 0.05 ~ "*",
      .default = "n.s."
    ))
  
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
    theme_minimal(base_size = 16, base_family = 'Source Sans Pro')+
    scale_fill_manual(values = thematic::okabe_ito(6), labels = c(toupper(all_combis$diseases[mydisease]), "Control")) + guides(fill=guide_legend(title=NULL), color="none") +
    scale_color_manual(values = thematic::okabe_ito(6)) +
    #scale_color_brewer(palette = "Set1") + # Choose a color palette
    #scale_fill_brewer(palette = "Set1", labels = c(toupper(all_combis$diseases[mydisease]), "Control")) + guides(fill=guide_legend(title=NULL), color="none")+
    my_base_theme()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p2)
  p2_listnew[[mydisease]] <- p2
  if (SAVE.files ==TRUE) {
    fig2_disease_vogel2013_path <- glue("/mnt/users/reich/rockerprojects/bestageing2022/output/plots/fig02vogel2013/{toupper(all_combis$diseases[mydisease])}_topdysregulatedmiRetrievemiR.svg")
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
  fig2vogel2013_path <- glue("/mnt/users/reich/rockerprojects/bestageing2022/output/plots/fig02vogel2013/arranged_topdysregulatedmiRetrievemiR.svg")
  ggsave(filename = fig2vogel2013_path, plot = fig2vogel2013, 
         width = 12, height = 12, 
         units = "in"  # default
  )
}
