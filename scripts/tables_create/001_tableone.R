

# create table one --------------------------------------------------------


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

require(readxl, lib.loc = lib_path)
require(openxlsx, lib.loc = lib_path)
require(glue, lib.loc = lib_path)
require(dplyr, lib.loc = lib_path)
require(skimr, lib.loc = lib_path)
require(tableone, lib.loc = lib_path)
require(gt, lib.loc = lib_path)
require(kableExtra, lib.loc = lib_path)

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


# TABLE ONE ---------------------------------------------------------------

myVars <- c("disease", "sex", "age", "weight", "height", "bmi", "smoke", 
            "dm", "hypertension", "dyslipo", "familiyhist", 
            "echo_ef", "echo_lv_dil_henry", "echo_lvedd", 
            "crea", "bnp", "troponinI", "troponinT", "hemo", "wbc", "cholesterol", 
            "hdl", "ldl")
catVars <- c("disease", "sex", "smoke", 
             "dm", "hypertension", "dyslipo", "familiyhist")

nonnormal_features <- c("bmi", "crea", "bnp", "troponinI", "troponinT", "hemo", "wbc")

#tableone::CreateTableOne(data=clean_all_meta, vars = myVars, factorVars = catVars)
#tableone <-  tableone::CreateTableOne(data=clean_all_meta, vars = myVars, factorVars = catVars)

# print(tableone, nonnormal = nonnormal_features, showAllLevels = TRUE, formatOptions = list(big.mark = ","))


clean_all_meta <- clean_all_meta %>% select(
  patID, disease, sex, 
  age, weight, height, bmi,
  smoke, dm, 
  hypertension, dyslipo, familiyhist,
  echo_ef, echo_lv_dil_henry, echo_lvedd, echo_lvesd, echo_lv_dil, echo_valv,
  crea, bnp, troponinI, troponinT, hemo, wbc, cholesterol, hdl, ldl
) %>%
  filter(disease != "pef") %>% 
  mutate(
    smoke = case_when(
      smoke == "Never" | smoke == "No" | smoke == "Past" ~ "No",
      smoke == "Present" | smoke == "Yes" ~ "Yes",
      .default = NA
    ),
    dm = case_when(
      dm == "Diab (cont diet)" | dm == "Diab (insul dep)" | dm == "Diab (oral med)" | dm == "Diab (unsp)" | dm == "New diagnosed" ~ "Present",
      dm == "No" | dm == "0" ~ "Not Present",
      .default = NA
    ),
    familiyhist = ifelse(familiyhist == "Yes", "Positive", "Negative"),
    dyslipo = case_when(
      dyslipo == "Yes" ~ "Yes",
      dyslipo == "No" ~ "No",
      .default = NA
    ),
    hypertension = case_when(
      hypertension == "Yes" ~ "Yes",
      hypertension == "No" ~ "No",
      .default = NA
    )
  )


clean_all_meta_tab01 <- clean_all_meta %>% 
  filter(patID %in% all_mirnas$pat_id) %>% 
  select(-patID)

mean_age_overall <- mean(clean_all_meta_tab01$age, na.rm=TRUE)
sd_age_overall <- sd(clean_all_meta_tab01$age, na.rm=TRUE)

mean_sex_overall <- mean(clean_all_meta_tab01$sex == "Male", na.rm=TRUE)
sum_men <- sum(clean_all_meta_tab01$sex == "Male", na.rm=TRUE)


# add strata
tableone_strata <- tableone::CreateTableOne(data=clean_all_meta_tab01, 
                                            vars = myVars[-1], # remove disease from vars
                                            factorVars = catVars[-1], 
                                            strata="disease"
)

kableone(tableone_strata) %>% 
  kable_classic_2(full_width = F)

# save
tab1Mat <- print(tableone_strata, nonnormal=nonnormal_features, quote=FALSE, noSpaces=TRUE, printToggle=FALSE)
tab1Mat

## Save to a CSV file and prepare it for final use ;)
write.csv2(tab1Mat, file = glue("{data_path_bestageing2022}/output/tables/tableone/{Sys.Date()}-table1strata_disease.csv"))
#openxlsx::write.xlsx(tab1Mat, file = "/mnt/users/reich/rockerprojects/bestageing2022/output/tables/tableone/table1strata_disease.xlsx")

# to gt

tableone_strata_gt <- as.data.frame(print(tableone_strata, nonnormal=nonnormal_features, quote=FALSE, noSpaces=TRUE))
tableone_strata_gt$variable <- rownames(tableone_strata_gt)
rownames(tableone_strata_gt) <- NULL
gt_table <- gt(tableone_strata_gt %>% dplyr::select(variable, everything()))
gt_table
