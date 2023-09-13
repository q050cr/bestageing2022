# ---------------------------------- #
# describe patient cohort # --------------------------------------------------
# ---------------------------------- #

# creates file `clean_all_meta` in "/mnt/users/reich/BestAgeing/output_new/EDA/clean_all_meta.rds"

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

source(glue("{data_path_BestAgeing}/scripts/_prepare_metadata.R"))

# get overview ------------------------------------------------------------
skimr::skim(clean_all_meta)
colnames(clean_all_meta)
# after inspecting NAs select reasonable cols to do further exploration

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

# saveRDS(clean_all_meta, file = glue("{data_path_BestAgeing}/output_new/EDA/clean_all_meta.rds"))


