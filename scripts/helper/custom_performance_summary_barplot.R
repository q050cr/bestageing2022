


# performance bar plot function
# here only one model performance shown!

library(dplyr)
library(tidyr)
library(ggplot2)
library(epoxy)

my_performance_plot <- function(mydata, ... ){
  # Filter top 5 models by AUC for each disease
  top_models <- mydata %>%
    arrange(desc(AUC), desc(f1.score)) %>%
    slice_head(n = 5)
  
  # Reshape the data to long format
  long_data <- top_models %>%
    select(model, AUC, f1.score, sensitivity, specificity, npv, ppv, accuracy) %>%
    pivot_longer(cols = c(AUC, f1.score, sensitivity, specificity, npv, ppv, accuracy), 
                 names_to = "metric", values_to = "value") %>% 
    mutate(metric = ifelse(metric == "f1.score", "F1-score", metric),
           metric = ifelse(metric == "accuracy", "Accuracy", metric),
           metric = ifelse(metric == "sensitivity", "Sensitivity", metric),
           metric = ifelse(metric == "specificity", "Specificity", metric),
           metric = ifelse(metric == "npv", "NPV", metric),
           metric = ifelse(metric == "ppv", "PPV", metric))
  
  # Plot with ggplot2
  performance_barplot <- ggplot(long_data, aes(x = metric, y = value, fill = model,
                                               #group = model)
  )) +
    geom_bar(stat = "identity", position = "dodge", alpha=0.8) +
    #facet_wrap(~ metric, scales = "free_y") + # Use facet_wrap just for y-axis variations
    ylab("Value") +
    xlab("Disease") +
    labs(fill = "Model") +
    labs(fill=NULL)+
    xlab(NULL)+
    ylab(NULL)+
    theme(legend.position = "top")+
    coord_flip()+
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6)[3]) +
    scale_y_continuous(breaks = c(0.00, 0.2, 0.4, 0.6, 0.8, 1)) + # Set the desired breaks for the y-axis
    theme(legend.position = "none")+  # remove legend for now (only one model, if more models compared remove this line)
    my_base_theme()
  return(performance_barplot)
}

# performance_barplot <- my_performance_plot(mydata = performance.summary.table.sens09)
# performance_barplot

# ggsave(filename = epoxy("./output/plots/performance_sum_plot/BACC_external/performance_summary_barplot.svg"), 
#        plot = performance_barplot, 
#        width = 6, height = 3, 
#        units = "in"
# )
