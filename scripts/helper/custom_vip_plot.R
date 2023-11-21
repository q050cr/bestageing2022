
library(ggplot2)
library(thematic)
library(dplyr)
library(forcats)


convert_mir_name <- function(name) {
  # Replace 'mir' with 'miR'
  name <- gsub("mir", "miR", name)
  
  # Replace underscores with hyphens
  name <- gsub("_", "-", name)
  
  return(name)
}
convert_mir_name_V <- Vectorize(convert_mir_name)

# custom TMwR plot function p 292
ggplot_imp <- function(...) {
  obj <- list(...)
  metric_name <- attr(obj[[1]], "loss_name")
  metric_lab <- paste(metric_name, 
                      "after permutations\n(Global Feature Importance)")
  
  full_vip <- bind_rows(obj) %>%
    filter(variable != "_baseline_" & variable != "pat_id") %>% 
    mutate(variable = convert_mir_name_V(variable))
  
  perm_vals <- full_vip %>% 
    filter(variable == "-full-model-") %>% 
    group_by(label) %>% 
    summarise(dropout_loss = mean(dropout_loss))
  
  p <- full_vip %>%
    filter(variable != "-full-model-") %>% 
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
