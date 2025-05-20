


### INFO ----------------------------------------------------------------------
# Differential miRNA Expression Analysis
# this script is sourced from `scripts/render_param_reports.R`
# selection provided by `all_combis$diseases` and `all_combis$analysis`

# script creates plots: "fig01vogel2013", "fig02vogel2013"


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
library(RColorBrewer, lib.loc = lib_path)
library(ggdist, lib.loc = lib_path)
library(gghalves, lib.loc = lib_path)
library(ggrepel, lib.loc = lib_path)
library(rstatix, lib.loc = lib_path)
library(ggthemes, lib.loc = lib_path)
library(ggpubr, lib.loc = lib_path)
library(pROC, lib.loc = lib_path)
library(reshape2, lib.loc = lib_path)
conflicted::conflict_prefer("expand", "tidyr")
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(janitor::make_clean_names)


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

# DESCRIPTIVE -----------------------------------------------------------------

clean_all_meta %>% 
  filter(disease != "pef") %>% 
  filter(disease == "ref") %>% 
  select(age) %>% 
  summary()

# center recruitment

addmargins(table(sub("-[0-9]+$", "", clean_all_meta$patID), clean_all_meta$disease))


clean_all_meta %>% 
  filter(patID %in% all_mirnas$pat_id) %>% 
  mutate(
    abbrev_center_name = (sub("-[0-9]+$", "", patID)),
    full_center_name = case_when(
      abbrev_center_name == "AMC" ~  "Academisch Medisch Centrum bij de Universiteit van Amsterdam",
      abbrev_center_name == "GUF" ~  "JOHANN WOLFGANG GOETHE UNIVERSITAET FRANKFURT AM MAIN",
      abbrev_center_name == "INSERM" ~  "INSTITUT NATIONAL DE LA SANTE ET DE LA RECHERCHE MEDICALE",
      abbrev_center_name == "NSC" ~  "National Scientific Center Institute of cardiology n.a. M.D.Strazhesko",
      abbrev_center_name == "SERMAS" ~  "SERVICIO MADRILENO DE SALUD",
      abbrev_center_name == "SFNH" ~  "AZIENDA COMPLESSO OSPEDALIERO SAN FILIPPO NERI",
      abbrev_center_name == "UCSC" ~  "UNIVERSITA CATTOLICA DEL SACRO CUORE",
      abbrev_center_name == "UKL-HD" ~  "UNIVERSITAETSKLINIKUM HEIDELBERG", 
      abbrev_center_name == "UNIPD" ~  "UNIVERSITA DEGLI STUDI DI PADOVA", 
      abbrev_center_name == "UOA" ~  "NATIONAL AND KAPODISTRIAN UNIVERSITY OF ATHENS", 
      abbrev_center_name == "UU" ~  "UPPSALA UNIVERSITET"
    )
  ) %>% 
  select(full_center_name, disease) %>% 
  filter(disease!= "pef") %>% 
  table() %>% 
  addmargins() -> center_statistics

center_statistics

write.csv2(center_statistics, file = glue("{data_path_bestageing2022}/output/tables/centers/center_mirna_stats.csv") )

center_statistics_with_abbrev <- clean_all_meta %>% 
  filter(patID %in% all_mirnas$pat_id) %>% 
  mutate(
    abbrev_center_name = (sub("-[0-9]+$", "", patID)),
    full_center_name = case_when(
      abbrev_center_name == "AMC" ~  "Academisch Medisch Centrum bij de Universiteit van Amsterdam",
      abbrev_center_name == "GUF" ~  "JOHANN WOLFGANG GOETHE UNIVERSITAET FRANKFURT AM MAIN",
      abbrev_center_name == "INSERM" ~  "INSTITUT NATIONAL DE LA SANTE ET DE LA RECHERCHE MEDICALE",
      abbrev_center_name == "NSC" ~  "National Scientific Center Institute of cardiology n.a. M.D.Strazhesko",
      abbrev_center_name == "SERMAS" ~  "SERVICIO MADRILENO DE SALUD",
      abbrev_center_name == "SFNH" ~  "AZIENDA COMPLESSO OSPEDALIERO SAN FILIPPO NERI",
      abbrev_center_name == "UCSC" ~  "UNIVERSITA CATTOLICA DEL SACRO CUORE",
      abbrev_center_name == "UKL-HD" ~  "UNIVERSITAETSKLINIKUM HEIDELBERG", 
      abbrev_center_name == "UNIPD" ~  "UNIVERSITA DEGLI STUDI DI PADOVA", 
      abbrev_center_name == "UOA" ~  "NATIONAL AND KAPODISTRIAN UNIVERSITY OF ATHENS", 
      abbrev_center_name == "UU" ~  "UPPSALA UNIVERSITET"
    )
  ) %>% 
  select(abbrev_center_name, full_center_name) %>% 
  group_by(abbrev_center_name, full_center_name) %>% 
  summarize(n=n())

center_statistics_with_abbrev
  