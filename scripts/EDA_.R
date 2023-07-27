# ---------------------------------- #
# describe patient cohort # --------------------------------------------------
# ---------------------------------- #

# creates file `clean_all_meta` in "/mnt/users/reich/BestAgeing/output_new/EDA/clean_all_meta.rds"


require(readxl, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(openxlsx, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(dplyr, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(skimr, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(tableone, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(gt, lib.loc = "/mnt/users/reich/programs/R43/lib")
require(kableExtra, lib.loc = "/mnt/users/reich/programs/R43/lib")

source("/mnt/users/reich/BestAgeing/scripts/_prepare_metadata.R")

# get overview ------------------------------------------------------------
skimr::skim(clean_all_meta)
colnames(clean_all_meta)
# after inspecting NAs select reasonable cols to do further exploration

clean_all_meta <- clean_all_meta %>% select(
  disease, sex, 
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
    
#table_clean_all_meta <- skimr::skim(clean_all_meta)

#table_ref_meta <- skimr::skim(ref_meta)

#table_clean_all_meta$complete_rate <- round(table_clean_all_meta$complete_rate,2)
#table_clean_all_meta$numeric.mean <- round(table_clean_all_meta$numeric.mean,2)
#table_clean_all_meta$numeric.sd <- round(table_clean_all_meta$numeric.sd,2)
#
#table_clean_all_meta <- table_clean_all_meta %>% 
#  select(skim_variable, n_missing, complete_rate, numeric.mean, numeric.sd, numeric.hist)

#write_csv(table_clean_all_meta, "output_new/EDA/skimr_all_meta.csv")
# write_csv(table_ref_meta , "output_new/EDA/skimr_ref_meta.csv")



#tableone::CreateTableOne(data=clean_all_meta)

dput(names(clean_all_meta)) #INTERESTING

# saveRDS(clean_all_meta, file = "/mnt/users/reich/BestAgeing/output_new/EDA/clean_all_meta.rds")

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

# add strata
tableone_strata <- tableone::CreateTableOne(data=clean_all_meta, 
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
write.csv2(tab1Mat, file = "/mnt/users/reich/rockerprojects/bestageing2022/output/tables/tableone/table1strata_disease.csv")
#openxlsx::write.xlsx(tab1Mat, file = "/mnt/users/reich/rockerprojects/bestageing2022/output/tables/tableone/table1strata_disease.xlsx")

# to gt

tableone_strata_gt <- as.data.frame(print(tableone_strata, nonnormal=nonnormal_features, quote=FALSE, noSpaces=TRUE))
tableone_strata_gt$variable <- rownames(tableone_strata_gt)
rownames(tableone_strata_gt) <- NULL
gt_table <- gt(tableone_strata_gt %>% dplyr::select(variable, everything()))
gt_table
