
library(ggridges)
library(binom)

# easy calibration plot

calibration_plot_easy <- function(predicted_probabilities, observed_outcomes, outcome, conf.level=0.9, conf=TRUE, ...) {
  # the fn returns: 
  # fn_dat ... a dataframe of the model output probabilites and calibrated probabilities needed to create calibration plot
  # calib_plot ... the calibration plot
  # ridges_probs_plot ... the ridge plot of the predicted probabilities
  # ridges_probs_plot_recalib ... the ridge plot of the recalibrated probabilities
  # density_probs_plot ... the density plot of the predicted probabilities
  # density_probs_plot_recalib ... the density plot of the recalibrated probabilities
  
  fn_dat <- data.frame(predicted_probabilities, observed_outcomes)
  
  # bin the predicted probs
  fn_dat$bin <- cut(fn_dat$predicted_probabilities, breaks = seq(0, 1, 0.1), include.lowest=TRUE, right=FALSE)
  
  # Calculate Observed Frequencies and Mean Predicted Probabilities:
  calibration_data <- fn_dat %>%
    group_by(bin) %>%
    summarize(count = n(),
              events = sum(observed_outcomes==outcome),
              mean_predicted = mean(predicted_probabilities),
              observed_frequency = mean(observed_outcomes==outcome)) %>% 
    mutate(source = "original")
  
  # Compute 90% confidence intervals
  conf_int <- binom::binom.confint(calibration_data$events, calibration_data$count, conf.level = conf.level, methods = "logit")
  
  
  # create the calibration plot ------------------------------------------------
  if (conf==TRUE){
    calib_plot <- ggplot(calibration_data, aes(x = mean_predicted, y = observed_frequency)) +
      geom_point(size = 2.5, alpha=0.7, shape=16) + # , color = "blue") +
      geom_line(aes(group = 1), alpha=0.7) +
      # conf int
      geom_ribbon(aes(ymin = conf_int$lower, ymax = conf_int$upper), fill = thematic::okabe_ito(6)[3], alpha = 0.5) +
      geom_line(aes(x = mean_predicted, y = mean_predicted), color = "red", linetype = "dashed") +
      geom_rug(data = subset(fn_dat, observed_outcomes == outcome), aes(x = predicted_probabilities, y=NULL), sides = "b", color = "black") +
      #xlim(c(0, 1)) +
      ylim(c(0, 1)) +
      xlab("Mean Predicted Probability") +
      ylab("Observed Frequency") +
      ggtitle(NULL) +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = thematic::okabe_ito(6)) +
      my_base_theme()
  } else {
    calib_plot <- ggplot(calibration_data, aes(x = mean_predicted, y = observed_frequency)) +
      geom_point(size = 2.5, alpha=0.7, shape=16) + # , color = "blue") +
      geom_line(aes(group = 1), alpha=0.7) +
      # conf int
      # geom_ribbon(aes(ymin = conf_int$lower, ymax = conf_int$upper), fill = thematic::okabe_ito(6)[3], alpha = 0.5) +
      geom_line(aes(x = mean_predicted, y = mean_predicted), color = "red", linetype = "dashed") +
      geom_rug(data = subset(fn_dat, observed_outcomes == outcome), aes(x = predicted_probabilities, y=NULL), sides = "b", color = "black") +
      #xlim(c(0, 1)) +
      ylim(c(0, 1)) +
      xlab("Mean Predicted Probability") +
      ylab("Observed Frequency") +
      ggtitle(NULL) +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = thematic::okabe_ito(6)) +
      my_base_theme()
  } 
  
  # geom ridges
  ridges_probs_plot <- ggplot(fn_dat, aes(x = predicted_probabilities, y = observed_outcomes, fill = observed_outcomes)) +
    geom_density_ridges(alpha=0.7) +
    labs(title = NULL,
         x = "Predicted Probability",
         y = NULL,
         fill = NULL) +
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6), name=NULL) +
    theme(axis.text.y = element_blank()) +
    my_base_theme()
  
  # geom_density
  density_probs_plot <- ggplot(fn_dat, aes(x = predicted_probabilities, fill = observed_outcomes)) +
    geom_density(alpha = 0.7) +
    xlim(0, 1) +
    labs(title = NULL,
         x = "Predicted Probability",
         y = "Density",
         fill = "Class") +
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6), name=NULL) +
    my_base_theme()
  
  return(list(fn_dat %>% as_tibble(), calib_plot, ridges_probs_plot, density_probs_plot))
}

# calibration_plot_easy(predicted_probabilities = predictions_loop[[models]][[".pred_died"]], observed_outcomes = predictions_loop[[models]]$o_mortality, outcome="died")


# recalibration plot -----------------------------------------------------------

custom_calibration_plot <- function(predicted_probabilities, observed_outcomes, outcome, conf.level=0.9, conf=TRUE, recalibration = TRUE, ...) {
  # the fn returns: 
  # fn_dat ... a dataframe of the model output probabilites and calibrated probabilities needed to create calibration plot
  # calib_plot ... the calibration plot
  # ridges_probs_plot ... the ridge plot of the predicted probabilities
  # ridges_probs_plot_recalib ... the ridge plot of the recalibrated probabilities
  # density_probs_plot ... the density plot of the predicted probabilities
  # density_probs_plot_recalib ... the density plot of the recalibrated probabilities

  fn_dat <- data.frame(predicted_probabilities, observed_outcomes)
  
  # recalibrate with glm/ Platt scaling ----------------------------------------
  # https://www.tidyverse.org/blog/2022/11/model-calibration/ 
  # https://medium.com/@chandu.bathula16/machine-learning-concept-68-platts-scaling-b8245421739e
  recalibrated_model <- glm(observed_outcomes~predicted_probabilities, family = binomial(link = "logit"))
  preds_recalibrated <- predict(recalibrated_model, data.frame(predicted_probabilities=predicted_probabilities), type = "response")
  fn_dat$preds_recalibrated = preds_recalibrated
  
  # bin the predicted probs
  fn_dat$bin <- cut(fn_dat$predicted_probabilities, breaks = seq(0, 1, 0.1), include.lowest=TRUE, right=FALSE)
  fn_dat$bin_recalib <- cut(fn_dat$preds_recalibrated, breaks = seq(0, 1, 0.1), include.lowest=TRUE, right=FALSE)
  
  # Calculate Observed Frequencies and Mean Predicted Probabilities:
  calibration_data <- fn_dat %>%
    group_by(bin) %>%
    summarize(count = n(),
              events = sum(observed_outcomes==outcome),
              mean_predicted = mean(predicted_probabilities),
              observed_frequency = mean(observed_outcomes==outcome)) %>% 
    mutate(source = "original")
  
  calibration_data_recalibration <- fn_dat %>%
    group_by(bin_recalib) %>%
    summarize(count = n(),
              events = sum(observed_outcomes==outcome),
              mean_predicted = mean(preds_recalibrated),
              observed_frequency = mean(observed_outcomes==outcome)) %>% 
    mutate(source = "recalibrated")
  
  # Compute 90% confidence intervals
  conf_int <- binom::binom.confint(calibration_data$events, calibration_data$count, conf.level = conf.level, methods = "logit")
  conf_int_recalib <- binom::binom.confint(calibration_data_recalibration$events, calibration_data_recalibration$count, conf.level = conf.level, methods = "logit")
  conf_int_combined <- rbind(conf_int, conf_int_recalib)
  # combined
  combined_data <- bind_rows(calibration_data, calibration_data_recalibration)
  
  
  # create the calibration plot ------------------------------------------------
  if (conf==TRUE & recalibration != TRUE){
    calib_plot <- ggplot(calibration_data, aes(x = mean_predicted, y = observed_frequency)) +
      geom_point(size = 2.5, alpha=0.7, shape=16) + # , color = "blue") +
      geom_line(aes(group = 1), alpha=0.7) +
      # conf int
      geom_ribbon(aes(ymin = conf_int$lower, ymax = conf_int$upper), fill = thematic::okabe_ito(6)[3], alpha = 0.5) +
      geom_line(aes(x = mean_predicted, y = mean_predicted), color = "red", linetype = "dashed") +
      geom_rug(data = subset(fn_dat, observed_outcomes == outcome), aes(x = predicted_probabilities, y=NULL), sides = "b", color = "black") +
      #xlim(c(0, 1)) +
      ylim(c(0, 1)) +
      xlab("Mean Predicted Probability") +
      ylab("Observed Frequency") +
      ggtitle(NULL) +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = thematic::okabe_ito(6)) +
      my_base_theme()
  } else if (conf!=TRUE & recalibration != TRUE){
    calib_plot <- ggplot(calibration_data, aes(x = mean_predicted, y = observed_frequency)) +
      geom_point(size = 2.5, alpha=0.7, shape=16) + # , color = "blue") +
      geom_line(aes(group = 1), alpha=0.7) +
      # conf int
      # geom_ribbon(aes(ymin = conf_int$lower, ymax = conf_int$upper), fill = thematic::okabe_ito(6)[3], alpha = 0.5) +
      geom_line(aes(x = mean_predicted, y = mean_predicted), color = "red", linetype = "dashed") +
      geom_rug(data = subset(fn_dat, observed_outcomes == outcome), aes(x = predicted_probabilities, y=NULL), sides = "b", color = "black") +
      #xlim(c(0, 1)) +
      ylim(c(0, 1)) +
      xlab("Mean Predicted Probability") +
      ylab("Observed Frequency") +
      ggtitle(NULL) +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = thematic::okabe_ito(6)) +
      my_base_theme()
  } else if (conf==TRUE & recalibration == TRUE) {
    # a) method geom_smooth
    calib_plot <- ggplot(combined_data, aes(x = mean_predicted, y = observed_frequency)) +
      geom_point(size = 2.5, alpha=0.7, shape=16) + # , color = "blue") +
      #geom_line(aes(group = source), alpha=0.7) +
      # conf int
      geom_smooth(aes(color = source), method = NULL, se = TRUE) + 
      geom_line(aes(x = mean_predicted, y = mean_predicted), color = "red", linetype = "dashed") +
      geom_rug(data = subset(fn_dat, observed_outcomes == outcome), aes(x = predicted_probabilities, y=NULL), sides = "b", color = "black") +
      geom_rug(data = subset(fn_dat, observed_outcomes == outcome), aes(x = preds_recalibrated, y=NULL), sides = "t", color = "black") +
      #xlim(c(0, 1)) +
      ylim(c(0, 1)) +
      xlab("Mean Predicted Probability") +
      ylab("Observed Frequency") +
      ggtitle(NULL) +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_color_manual(values = thematic::okabe_ito(6), name=NULL) +
      scale_fill_manual(values = thematic::okabe_ito(6), name=NULL) +
      my_base_theme()
    
    # b) method geom_ribbon
    calib_plot <- ggplot(combined_data, aes(x = mean_predicted, y = observed_frequency)) +
      geom_point(size = 2.5, alpha=0.7, shape=16) + # , color = "blue") +
      geom_line(aes(group = source), alpha=0.7) +
      # conf int
      geom_ribbon(aes(ymin = conf_int_combined$lower, ymax = conf_int_combined$upper, group=source, fill = source),  alpha = 0.5) +
      geom_line(aes(x = mean_predicted, y = mean_predicted), color = "red", linetype = "dashed") +
      geom_rug(data = subset(fn_dat, observed_outcomes == outcome), aes(x = predicted_probabilities, y=NULL), sides = "b", color = "black") +
      geom_rug(data = subset(fn_dat, observed_outcomes == outcome), aes(x = preds_recalibrated, y=NULL), sides = "t", color = "black") +
      #xlim(c(0, 1)) +
      ylim(c(0, 1)) +
      xlab("Mean Predicted Probability") +
      ylab("Observed Frequency") +
      ggtitle(NULL) +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = thematic::okabe_ito(6), name=NULL) +
      my_base_theme()
  }
  
  # geom ridges
  ridges_probs_plot <- ggplot(fn_dat, aes(x = predicted_probabilities, y = observed_outcomes, fill = observed_outcomes)) +
    geom_density_ridges(alpha=0.7) +
    labs(title = NULL,
         x = "Predicted Probability",
         y = NULL,
         fill = NULL) +
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6), name=NULL) +
    theme(axis.text.y = element_blank()) +
    my_base_theme()
  ridges_probs_plot_recalib <- ggplot(fn_dat, aes(x = preds_recalibrated, y = observed_outcomes, fill = observed_outcomes)) +
    geom_density_ridges(alpha=0.7) +
    labs(title = NULL,
         x = "Predicted Probability",
         y = NULL,
         fill = NULL) +
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6), name=NULL) +
    theme(axis.text.y = element_blank()) +
    my_base_theme()
  
  # geom_density
  density_probs_plot <- ggplot(fn_dat, aes(x = predicted_probabilities, fill = observed_outcomes)) +
    geom_density(alpha = 0.7) +
    xlim(0, 1) +
    labs(title = NULL,
         x = "Predicted Probability",
         y = "Density",
         fill = "Class") +
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6), name=NULL) +
    my_base_theme()
  
  density_probs_plot_recalib <- ggplot(fn_dat, aes(x = preds_recalibrated, fill = observed_outcomes)) +
    geom_density(alpha = 0.7) +
    xlim(0, 1) +
    labs(title = NULL,
         x = "Predicted Probability",
         y = "Density",
         fill = "Class") +
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6), name=NULL) +
    my_base_theme()
  
  return(list(fn_dat %>% as_tibble(), calib_plot, ridges_probs_plot, ridges_probs_plot_recalib, density_probs_plot, density_probs_plot_recalib))
}

# custom_calibration_plot(predicted_probabilities = predictions_loop[[models]][[".pred_died"]], observed_outcomes = predictions_loop[[models]]$o_mortality, outcome="died",
#                         recalibration=TRUE)


