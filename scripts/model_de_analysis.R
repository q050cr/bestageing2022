
### INFO ----------------------------------------------------------------------
# Differential miRNA Expression Analysis
# this script is sourced from `scripts/render_param_reports.R`
# selection provided by `all_combis$diseases` and `all_combis$analysis`

# dependencies ---------------------------------------------------------------
library(readxl, lib.loc = "/mnt/users/reich/programs/R/lib")
library(janitor, lib.loc = "/mnt/users/reich/programs/R/lib")
library(dplyr, lib.loc = "/mnt/users/reich/programs/R/lib")
library(tidyr, lib.loc = "/mnt/users/reich/programs/R/lib")
library(stringr, lib.loc = "/mnt/users/reich/programs/R/lib")
library(purrr, lib.loc = "/mnt/users/reich/programs/R/lib")
library(ggplot2, lib.loc = "/mnt/users/reich/programs/R/lib")
library(ggrepel, lib.loc = "/mnt/users/reich/programs/R/lib")
library(ggthemes, lib.loc = "/mnt/users/reich/programs/R/lib")
conflicted::conflict_prefer("expand", "tidyr")
conflicted::conflicts_prefer(dplyr::filter)

if ( !exists("all_combis")  ) {  # check if variable name exists in env
  diseases <- c("dcm", "acs", "cad", "hfref")
  analysis <- c("selected", "full")
  all_combis <- tidyr::crossing(diseases, analysis) #%>% 
  #filter(analysis=="selected")
}

# load data ---------------------------------------------------------------

# MIRNA DAT
model_data1 <- clean_names(readRDS(file = '/mnt/users/reich/BestAgeing/data_new/model_data1.RDS'))  # has also multiclass col + diagnoses

load(file = "../../BestAgeing/data/mirnas.rda")  # "UKL-HD" n=765
load(file = "../../BestAgeing/data/data.rda")  # "UKL-HD" n=731
all_mirnas <- clean_names(mirnas)
names(all_mirnas) <- str_replace(string = names(all_mirnas), pattern = "mi_r", replacement = "mir")
mirnas_disease <- clean_names(data)
names(mirnas_disease) <- str_replace(string = names(mirnas_disease), pattern = "mi_r", replacement = "mir")
rm(data, mirnas)

# load-mirnas-from_research
# create vector of described mirnas
load("../../BestAgeing/data_research/fromR/researchMiRNAAccession.rda")
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
load(file = "../../BestAgeing/data/diagnoses_df.rda")

## SURVIVAL DAT
survival_dat <- clean_names(readRDS(file = 'data/202211908_XMELD_abfrage_best_ageing.rds'))# %>% 
# original path "../../XMeldPortal_neu/meldeportal-tools-meldeportalclient-9.3/Rout/202211908_XMELD_abfrage_best_ageing.rds"

## metadata from DB
# https://www.bestageing.org/Pages/Login.aspx?ReturnUrl=%2f&AspxAutoDetectCookieSupport=1
load(file = "../../BestAgeing/data/clean_all_meta.rda")  # created in "scripts/_prepare_metadata.R"
clean_all_meta <- clean_all_meta %>% 
  mutate(age = ifelse(age < 18, NA, age))  # wrong age remove
# cath data? "hkdb"

## load all original metadat xlsx files again to make sure that also overlapped 
#patients (e.g. dcm+cad) are in each group
control_ids <- read_excel("../../BestAgeing/data/pheno_controls.xlsx") %>% 
  dplyr::pull(BestAgeingCode)

# "UKL-HD-00318" both in Control and CAD dataset, looked it up (HK Nr 1289-2015): KHK ohne hg Stenosen, LV gut --> assign to CAD only
control_ids <- control_ids[control_ids != "UKL-HD-00318"]

###
# We apply both parametric t-tests and nonparametric U-tests ----------------

###
# START LOOP for specified combinations -------------------------------------
###

for (i in 1:nrow(all_combis[all_combis[["analysis"]] == "full", ])) {  # only consider all miRNAs for DE (not research miRNAs only)
  ## reassign disease since only full analysis here
  disease <- all_combis[all_combis[["analysis"]] == "full", ]$diseases[i]
  
  ## DATA FIRST
  # create parameter specific {disease}_ids
  filename <- paste0("../../BestAgeing/data/pheno_", disease, ".xlsx")
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
  # indexing
  cont.index <- data01$disease == "control"
  case.index <- data01$disease == disease
  # run tests
  for(miRNA in 1:(ncol(all_mirnas)-1)) {
    # update index since first colnames are [1] "pat_id"  "disease"  "age"   "sex"  "hsa_let_7a_3p" 
    miRNA_col <- miRNA+4
    cont <- data01[cont.index, miRNA_col] %>% as_vector()
    case <- data01[case.index, miRNA_col] %>% as_vector()
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
    pval.glm[miRNA] <- coef(summary(logreg))[2,4]
  }
  # we conducted 2549 t-tests and 2549 glm-models for each gene
  
  de.results <- tibble(miRNA = colnames(all_mirnas)[-1],  # all miRNAs without patID
                       # average.difference = average.difference,
                       log2FoldChange = log2FoldChange,
                       pval.t.test = pval.t.test,
                       pval.u.test = pval.u.test,
                       pval.glm = pval.glm)
  
  filename.de.tibble <-paste0("./output/de_results/", Sys.Date(), "_de_results_DISEASE_", toupper(disease), ".rds")
  saveRDS(object = de.results, file = filename.de.tibble)
  rm(de.results)
  
  print(paste0("|||-----------------------Run finished for disease: ", toupper(disease), " -----------------------|||"))
}


