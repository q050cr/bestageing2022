

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

fmt_stars <- function(data, columns, rows = NULL) {
  
  # Capture expression in `rows`
  rows <- rlang::enquo(rows)
  
  # Pass `data`, `columns`, `rows`, and the formatting
  # functions as a function list to `fmt()`
  fmt(
    data = data,
    columns = columns,
    rows = !!rows,
    fns = list(
      default = function(x) {
        dplyr::case_when(
          x < 0.005 ~ "***",
          x < 0.01  ~ "**",
          x < 0.05  ~ "*",
          TRUE      ~ ""
        )
      }
    )
  )
}

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
      smoke == "Never" | smoke == "No"  ~ "No",
      smoke == "Past" ~ "Ex", 
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



# save patients with miRNA data only ----------------------------------
clean_all_meta_2023_table01 <- clean_all_meta %>% 
  filter(patID %in% all_mirnas$pat_id)

saveRDS(object = clean_all_meta_2023_table01, file = glue("{data_path_bestageing2022}/data/disease_identifier_table01/20240123clean_all_meta_2023_table01.rds"))

# add information! 2024
pheno_control <- read.xlsx(xlsxFile = glue("{data_path_BestAgeing}/data/pheno_orig/BestAgeing_Report_30_01_2019-Controls.xlsx")) %>% 
  mutate(`Controls-CliExam-BloodPressSys` = as.numeric(`Controls-CliExam-BloodPressSys`),
         `Controls-CliExam-BloodPressDias` = as.numeric(`Controls-CliExam-BloodPressDias`)) %>% 
  mutate(hypertension2 = case_when(`Controls-CliExam-BloodPressSys` > 160 ~ "Yes",
                                  `Controls-CliExam-BloodPressDias` > 100 ~ "Yes",
                                  .default = "No"), 
         diabetes2 = `Controls-ExcCrit-Diab`) %>% 
  as_tibble()

pheno_control$`Controls-ExcCrit-Diab` %>% table()

clean_all_meta_2023_table01 <- 
  clean_all_meta_2023_table01 %>% 
  left_join(
    pheno_control %>% select(BestAgeingCode, hypertension2, diabetes2),
    by = c("patID" = "BestAgeingCode")
  ) %>% 
  mutate(hypertension = case_when(
    hypertension2 == "Yes" ~ "Yes",
    hypertension2 == "No" ~ "No",
    .default = hypertension
  )) %>% 
  mutate(dm = case_when(
    diabetes2 == "Yes" ~ "Present",
    diabetes2 == "No" ~ "Not Present",
    .default = dm
  )) %>% 
  select(-hypertension2, -diabetes2)


## ukhd data add -----------------------------------------------------------
master_ukhd <- clean_names(read.xlsx(xlsxFile = glue("{data_path_BestAgeing}/patientID/20180212 _Patientenliste_BestAgeing_Master_clean_nonames.xlsx"))) %>% 
  as_tibble()

master_ukhd_clean <- master_ukhd %>% 
  mutate(lv_ef_hkt = as.numeric(lv_ef_hkt),
         lv_ef_echo = as.numeric(lv_ef_echo)) %>% 
  mutate(
    lvef_cat = case_when(
      lv_ef_hkt == 0 ~ "good",
      lv_ef_hkt == 1 ~ "slightly reduced",
      lv_ef_hkt == 2 ~ "moderately reduced",
      lv_ef_hkt == 3 ~ "severely reduced",
      .default = na_chr)
  ) %>% 
  mutate(
    lvef_cat = case_when(  # echo schlägt hk ;)
      lv_ef_echo == 0 ~ "good",
      lv_ef_echo == 1 ~ "slightly reduced",
      lv_ef_echo == 2 ~ "moderately reduced",
      lv_ef_echo == 3 ~ "severely reduced",
      .default = lvef_cat)
  ) %>% 
  mutate(
    aht = case_when(
      aht == "1" ~ "Yes",
      aht == "2" ~ "No",
      .default = na_chr
    )
  ) %>% 
  mutate(
    hyperlipidamie = case_when(
      hyperlipidamie == "1" ~ "Yes",
      hyperlipidamie == "2" ~ "No",
      .default = na_chr
    )
  ) %>% 
  mutate(
    dm = case_when(
      dm == "1" ~ "Yes",
      dm == "2" ~ "No",
      .default = na_chr
    ),
    smoking = case_when(
      smoking == "1" ~ "Yes",
      smoking == "2" ~ "No",
      smoking == "3" ~ "Ex",
      .default = na_chr
    )
  ) %>% 
  # medication
  mutate(
    aceORarb = case_when(
      arb_at1_antagonisten_sartane == "1" ~ "Yes",
      arb_at1_antagonisten_sartane == "2" ~ "No",
      .default = na_chr
    ),
    aceORarb = case_when(
      ace_h == "1" ~ "Yes",
      ace_h == "2" ~ "No",
      .default = aceORarb
    ),
    bbl = case_when(
      bbl == "1" ~ "Yes",
      bbl == "2" ~ "No",
      .default = na_chr
    ),
    calcium_antagonist = case_when(
      calcium_antagonist == "1" ~ "Yes",
      calcium_antagonist == "2" ~ "No",
      .default = na_chr
    ),
    ass = case_when(
      ass == "1" ~ "Yes",
      ass == "2" ~ "No",
      .default = na_chr
    ),
    P2Y12inhibitors = case_when(
      clopidogrel_plavix == "1" | brilique_ticagrelor == "1"  | prasugrel_efient == "1"  ~ "Yes",
      .default = "No"
    ),
    orale_antikoagulant = case_when(
      orale_antikoagulant == "1" ~ "Yes",
      orale_antikoagulant == "2" ~ "No",
      .default = na_chr
    ),
    statin = case_when(
      statin == "1" ~ "Yes",
      statin == "2" ~ "No",
      .default = na_chr
    ),
  ) %>% 
  select(best_ageing_code, lvef_cat, aht, hyperlipidamie, dm, aceORarb, bbl, calcium_antagonist, ass, P2Y12inhibitors, orale_antikoagulant, statin) 



## lab 2023 add ------------------------------------------------------------

lab2023 <- read_excel(path = glue("{data_path_bestageing2022}/data/elham2023/2023-11-15-elham_metadat_combined_LAB2023_b.xlsx")) %>% 
  filter(best_ageing_code %in% survival_dat1$pat_id) %>% 
  select(best_ageing_code, HSTNT, HSTNTHP, NTBNP, INR, BILI, KREA, LEUKO, CRP, HB, CHOL) %>% 
  mutate(
    HSTNT = case_when(is.na(HSTNT) ~ HSTNTHP,
                      .default = HSTNT)
  ) %>% 
  # drop duplicates
  distinct(best_ageing_code, .keep_all = TRUE) %>% 
  # convert to numeric and change "," -> "." for decimals
  mutate(across(!c(best_ageing_code), function(x) as.numeric(gsub(",", ".", x)) )) %>% 
  select(-HSTNTHP) # is second TNT, first TNT measured in HSTNT, kinetics

master_ukhd_clean <- master_ukhd_clean %>% 
  left_join(lab2023, by=c("best_ageing_code" = "best_ageing_code") ) %>% 
  rename(dm2 = dm)

# combine
clean_all_meta_2023_table01_new <- clean_all_meta_2023_table01 %>% 
  left_join(
    master_ukhd_clean %>% 
      dplyr::distinct(best_ageing_code, .keep_all = TRUE), 
    by=c("patID"="best_ageing_code")
  ) %>% 
  mutate(
    bnp = case_when(
      is.na(bnp) ~ NTBNP,
      .default = bnp
    ),
    troponinT = case_when(
      is.na(troponinT) ~ HSTNT,
      .default = troponinT
    ),
    crea = case_when(
      is.na(crea) ~ KREA,
      .default = crea
    ),
    hemo = case_when(
      is.na(hemo) ~ HB,
      .default = hemo
    ),
    wbc = case_when(
      is.na(wbc) ~ LEUKO,
      .default = wbc
    ),
    cholesterol = case_when(
      is.na(cholesterol) ~ CHOL,
      .default = cholesterol
    )
  ) %>% 
  select(-c(CHOL, LEUKO, HB, KREA, HSTNT, NTBNP)) %>% 
  mutate(
    dm = case_when(
      dm2 == "Yes" ~ "Present",
      dm2 == "No"  ~ "Not Present",
      .default = dm
    ),
    hypertension = case_when(
      aht == "Yes" ~ "Yes",
      aht == "No"  ~ "No",
      .default = hypertension
    ),
    dyslipo = case_when(
      hyperlipidamie == "Yes" ~ "Yes",
      hyperlipidamie == "No"  ~ "No",
      .default = dyslipo
    ),
  ) %>% 
  mutate(
    lvef_cat = case_when(
      echo_ef > 52 ~ "good",
      echo_ef <=52 & echo_ef > 40 ~ "slightly reduced",
      echo_ef <=40 & echo_ef > 30 ~ "moderately reduced",
      echo_ef <=30 ~ "severely reduced",
      .default = lvef_cat)
  ) %>% 
  select(-c(dm2, aht, hyperlipidamie))

clean_all_meta_2023_table01_new <- clean_all_meta_2023_table01_new %>% 
  mutate(lvef_cat = factor(lvef_cat, levels = c("good", "slightly reduced", "moderately reduced", "severely reduced" ))) %>% 
  mutate(
    disease = case_when(
      disease == "control" ~ "Control",
      disease == "acs" ~ "ACS",
      disease == "cad" ~ "CAD",
      disease == "dcm" ~ "DCM",
      disease == "ref" ~ "HFrEF",
    ),
    disease = factor(disease, levels = c("Control", "ACS", "CAD", "DCM", "HFrEF"))
  )

## old table 01 prep -----------------------------------------------------

clean_all_meta_tab01 <- clean_all_meta_2023_table01 %>% 
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
write.csv2(tab1Mat, file = glue("{data_path_bestageing2022}/output/tables/tableone/table1strata_disease2023.csv"))
#openxlsx::write.xlsx(tab1Mat, file = "/mnt/users/reich/rockerprojects/bestageing2022/output/tables/tableone/table1strata_disease.xlsx")

# to gt

tableone_strata_gt <- as.data.frame(print(tableone_strata, nonnormal=nonnormal_features, quote=FALSE, noSpaces=TRUE))
tableone_strata_gt$variable <- rownames(tableone_strata_gt)
rownames(tableone_strata_gt) <- NULL
gt_table <- gt(tableone_strata_gt %>% dplyr::select(variable, everything()))
gt_table


# new 2024 table01 --------------------------------------------------------

table01_dat <- clean_all_meta_2023_table01_new %>% 
  select(-patID) %>% 
  rename(
    Disease = disease,
    Sex = sex,
    Age = age,
    Weight = weight,
    Height = height,
    BMI = bmi,
    SmokingStatus = smoke,
    Diabetes = dm,
    Hypertension = hypertension,
    Dyslipidemia = dyslipo,
    FamilyHistory = familiyhist,
    Echocardiogram_EF = echo_ef,
    LVEF_cat = lvef_cat,
    Creatinine = crea,
    `nt-proBNP` = bnp,
    `Troponin-T` = troponinT,
    Hemoglobin = hemo,
    Leukocytes = wbc,
    TotalCholesterol = cholesterol,
    HDL_Cholesterol = hdl,
    LDL_Cholesterol = ldl,
    `ACE-Inhibitor or ARB` = aceORarb,
    `Beta Blocker` = bbl,
    `Calcium Antagonist` = calcium_antagonist,
    Aspirin = ass,
    `P2Y12 Inhibitors` = P2Y12inhibitors,
    `Oral Anticoagulants` = orale_antikoagulant,
    Statin = statin,
    INR = INR,
    Bilirubin = BILI
  )

myVars <- c("Disease", "Sex", "Age", "Weight", "Height", "BMI", "SmokingStatus", 
            "Diabetes", "Hypertension", "Dyslipidemia", "FamilyHistory", 
            "Echocardiogram_EF", "LVEF_cat",
            # LAB
            "Creatinine", "nt-proBNP", "Troponin-T", "Hemoglobin", "Leukocytes", "TotalCholesterol", 
            "HDL_Cholesterol", "LDL_Cholesterol", "INR", "Bilirubin", "CRP",
            # meds
            "ACE-Inhibitor or ARB", 
            "Beta Blocker", "Calcium Antagonist", "Aspirin", "P2Y12 Inhibitors", 
            "Oral Anticoagulants", "Statin"
            )
catVars <- c("Disease", "Sex", "SmokingStatus", 
             "Diabetes", "Hypertension", "Dyslipidemia", "FamilyHistory",
             "ACE-Inhibitor or ARB", "Beta Blocker", "Calcium Antagonist", "Aspirin", "P2Y12 Inhibitors", "Oral Anticoagulants", "Statin")

nonnormal_features <- c("Creatinine", "nt-proBNP", "Troponin-T")


# add strata
tableone_strata <- tableone::CreateTableOne(data=table01_dat, 
                                            vars = myVars[-1], # remove disease from vars
                                            factorVars = catVars[-1], 
                                            strata="Disease"
)

kableone(tableone_strata, nonnormal=nonnormal_features) %>% 
  kable_classic_2(full_width = F)

# save
tab1Mat_new <- print(tableone_strata, nonnormal=nonnormal_features, quote=FALSE, noSpaces=TRUE, printToggle=FALSE)
tab1Mat_new

## Save to a CSV file and prepare it for final use ;)
write.csv2(tab1Mat_new, file = glue("{data_path_bestageing2022}/output/tables/tableone/20240123-table1strata_disease.csv"))
#openxlsx::write.xlsx(tab1Mat_new, file = "/mnt/users/reich/rockerprojects/bestageing2022/output/tables/tableone/table1strata_disease.xlsx")

# to gt

tableone_strata_gt <- as.data.frame(print(tableone_strata, nonnormal=nonnormal_features, quote=FALSE, noSpaces=TRUE))
tableone_strata_gt$variable <- rownames(tableone_strata_gt)
rownames(tableone_strata_gt) <- NULL
gt_table <- gt(tableone_strata_gt %>% dplyr::select(variable, everything()))
gt_table

gt_table %>% 
  gtsave(glue("{data_path_bestageing2022}/output/tables/tableone/20240123-table1strata_disease.docx")) 

