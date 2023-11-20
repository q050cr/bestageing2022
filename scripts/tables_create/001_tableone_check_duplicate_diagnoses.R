

load(glue("{data_path_BestAgeing}/data/diagnoses_df.rda"))

# load all_mirnas in another script

diagnoses_df <- diagnoses_df %>% 
  filter(patID %in% all_mirnas$pat_id) 

diagnoses_df %>% 
  select(disease) %>% 
  table()

diagnoses_df %>% 
  filter(patID %in% all_mirnas$pat_id) %>% 
  distinct(patID)

# Keeping rows that have duplicate patIDs
duplicates_df <- diagnoses_df %>%
  group_by(patID) %>%
  filter(n() > 1 ) %>%
  ungroup() %>% 
  arrange(patID)
duplicates_df

duplicates_df_dcm_ref <- diagnoses_df %>%
  group_by(patID) %>%
  filter(n() > 1 & "ref" %in% disease & "dcm" %in% disease) %>%  # for a given patID, both "ref" and "dcm" appear in the disease column in the group
  ungroup() %>% 
  arrange(patID)
duplicates_df_dcm_ref


# acs duplicates are combined with HFrEF!
duplicates_df_acs_duplicates <- diagnoses_df %>%
  group_by(patID) %>%
  filter(n() > 1 & "acs" %in% disease) %>%  # for a given patID, both "ref" and "dcm" appear in the disease column in the group
  ungroup() %>% 
  arrange(patID)
View(duplicates_df_acs_duplicates)

# STEMIS?
diagnoses_df %>% select(stemi_factor) %>% table()


# ref and ef > 45? cannot be due to inclusion criteria?
ref_subgroup_ef <- diagnoses_df %>% 
  filter(disease == "ref" & lvef > 45) %>% 
  arrange(patID)
View(ref_subgroup_ef)
