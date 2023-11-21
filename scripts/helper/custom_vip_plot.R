
library(ggplot2)
library(thematic)
library(dplyr)

# custom TMwR plot function p 292
ggplot_imp <- function(...) {
  obj <- list(...)
  metric_name <- attr(obj[[1]], "loss_name")
  metric_lab <- paste(metric_name, 
                      "after permutations\n(Global Feature Importance)")
  
  full_vip <- bind_rows(obj) %>%
    filter(variable != "_baseline_" & variable != "rapID" & variable != "source") %>% 
    mutate(variable = case_when(
      variable == "age"                 ~ "Age",
      variable == "t0_na_value"         ~ "Sodium",
      variable == "t0_inr_value"        ~ "INR",
      variable == "t0_hstnt_value"     ~ "hs-TroponinT",
      variable == "grace_score"        ~ "GRACE Score",
      variable == "vit_herzfrequenz"   ~ "Heart Rate",
      variable == "t0_crp_value"       ~ "CRP",
      variable == "t0_ck_value"        ~ "CK",
      variable == "t0_thrombo_value"   ~ "Thrombocytes",
      variable == "t0_gluc_value"      ~ "Glucose",
      variable == "t0_hb_value"        ~ "Hemoglobin",
      variable == "vit_rr_syst"        ~ "Systolic BP",
      variable == "t0_krea_value"      ~ "Creatinine",
      variable == "t0_leuko_value"     ~ "Leucocytes",
      variable == "symptombeginn"      ~ "Symptom Onset",
      variable == "h_khk"              ~ "CHD History",
      variable == "ekg_sinus_normal"   ~ "Normal Sinus ECG",
      variable == "aktiver_raucher"    ~ "Current Smoker",
      variable == "h_cholesterin"      ~ "Cholesterol History",
      variable == "sex_f1_m0"          ~ "Sex", # Assuming 1 is Female and 0 is Male
      variable == "h_familienana"      ~ "Family History",
      variable == "ekg_st_senkung"     ~ "ECG ST Depression",
      variable == "ekg_schrittmacher"  ~ "Pacemaker ECG",
      variable == "h_diabetes"         ~ "Diabetes History",
      variable == "h_lvdys_grad_BINARY"~ "LV Dysfunction",
      variable == "h_hypertonie"       ~ "Hypertension History",
      variable == "t0_ntbnp_value"     ~ "NTproBNP",
      variable == "t0_ckdepi_value"    ~ "GFR",
      .default = variable                 
    ))
  
  perm_vals <- full_vip %>% 
    filter(variable == "_full_model_") %>% 
    group_by(label) %>% 
    summarise(dropout_loss = mean(dropout_loss))
  
  p <- full_vip %>%
    filter(variable != "_full_model_") %>% 
    mutate(variable = fct_reorder(variable, dropout_loss)) %>%
    ggplot(aes(dropout_loss, variable)) 
  if(length(obj) > 1) {
    p <- p + 
      facet_wrap(vars(label)) +
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss, color = label),
                 size = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(aes(color = label, fill = label), alpha = 0.2)
  } else {
    p <- p + 
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss),
                 size = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(fill = thematic::okabe_ito(6)[1], alpha = 0.4)
    
  }
  p +
    theme(legend.position = "none",
          axis.text.y = element_text(size=4)) +
    labs(x = metric_lab, 
         y = NULL,  fill = NULL,  color = NULL) +
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6)) +
    scale_color_manual(values = thematic::okabe_ito(6)) +
    my_base_theme()
}
