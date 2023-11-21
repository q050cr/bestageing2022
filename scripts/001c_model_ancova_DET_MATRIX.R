
# run script 001c_model_de_an... first

convert_mir_name <- function(name) {
  # Replace 'mir' with 'miR'
  name <- gsub("mir", "miR", name)
  
  # Replace underscores with hyphens
  name <- gsub("_", "-", name)
  
  return(name)
}

convert_mir_name_V <- Vectorize(convert_mir_name)

# PREPARE DATA TO RUN ANOVA ----------------------------------------------------
clean_all_meta_2023_table01 <- readRDS(file = glue("{data_path_bestageing2022}/data/disease_identifier_table01/clean_all_meta_2023_table01.rds"))

anova_analysis_dat <- all_mirnas %>%
  inner_join(clean_all_meta_2023_table01 %>% select(patID, disease, sex, age), by = c("pat_id" = "patID")) %>% 
  relocate(disease, sex, age, .after = pat_id)

## FILTERING DET MATRIX ----------------------------------------------------
if (filterDetMatrix == TRUE){
  filter_dat <- anova_analysis_dat %>% select(-c(age,sex))
  
  det_mat_all_mirnas_tmp <- det_mat_all_mirnas %>% 
    filter(pat_id %in% filter_dat$pat_id) %>% 
    left_join(filter_dat %>% select(pat_id, disease), by=c("pat_id")) %>% 
    relocate(disease, .after = pat_id)
  
  print(dim(filter_dat) == dim(det_mat_all_mirnas_tmp))
  
  # EDA - Plotting a density of average log2 expressions 
  data_long <- filter_dat %>% 
    pivot_longer(cols = -c(disease, pat_id), names_to = "miRNA", values_to = "expression")
  
  # Stratified density plot
  density_plot_filter_eda_stratified <- ggplot(data_long, aes(x=expression, fill=disease)) + 
    geom_density(alpha=0.5) + 
    labs(x="Log2 Expression") +  # title="Density plot of log2 expressions", 
    facet_wrap(~ disease, nrow=2) +
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6), name=NULL) +
    theme(axis.text.y = element_blank())
  #print(density_plot_filter_eda_stratified)
  
  density_ridges_plot_filter_eda_stratified <- ggplot(data_long, aes(x=expression, y=disease, fill=disease)) + 
    geom_density_ridges(alpha=0.5) +  # , jittered_points = TRUE, position = position_points_jitter(width = 0.05)) + # too much data involved if points plotted
    labs(x="Log2 Expression", y="Density") +  # title="Density Ridge plot of log2 expressions", 
    theme_minimal(base_size = 16, base_family = 'Arial')+
    scale_fill_manual(values = thematic::okabe_ito(6), name=NULL) +
    theme(axis.text.y = element_blank())
  #print(density_ridges_plot_filter_eda_stratified)
  
  if (SAVE.files ==TRUE) {
    density_ridges_plot_filter_eda_stratified_path <- glue("{data_path_bestageing2022}/output/plots/preprocessing/filtering_miRNA/001c_ANOVA_density_ridges.svg")
    ggsave(filename = density_ridges_plot_filter_eda_stratified_path, plot = density_ridges_plot_filter_eda_stratified, 
           width = 8, height = 5, 
           units = "in"  # default
    )
  }
  
  threshold <- 0.8  # conservative threshold
  
  # Filter for control samples
  lowExprAll <- colMeans(det_mat_all_mirnas_tmp[, -(1:2)]) <= threshold
  sum(!lowExprAll)
  # Filter for disease samples
  #lowExprDisease <- colMeans(det_mat_all_mirnas_tmp[det_mat_all_mirnas_tmp$disease == disease, -(1:2)]) <= threshold
  
  # Identify miRNAs that are lowly expressed in both conditions
  #toFilter <- lowExprControl & lowExprDisease
  
  # Filter out those miRNAs from the dataset
  filteredData <- filter_dat[, c(TRUE, TRUE, !lowExprAll)] ## TRUE TRUE corresponds to keep pat_id and disease
  
  data01 <- anova_analysis_dat %>% 
    select(1:4) %>% 
    dplyr::bind_cols(filteredData[, -(1:2)])
}

## LIMMA BATCH EFFECT ------------------------------------------------------

if(runLimmaRemoveBatchEffect == TRUE) {
  batch_center <- remove_digits(data01$pat_id)
  table(batch_center, data01$disease)
  exprs_data <- as.matrix(data01[,-c(1:4)])
  
  # add metadata
  exprs_metadat <- data01 %>% 
    # batch 1
    mutate(batch_center = batch_center) %>% 
    # batch 2
    left_join(hbdx_metadat %>% select(customer_id, slide_id, sample_qc, array_qc),
              by = c("pat_id"="customer_id"))
  batch_slide_id <- exprs_metadat %>% pull(slide_id)
  
  # pheno data
  pheno_data <- data01 %>% 
    select(disease, age, sex) %>% 
    as.data.frame()
  rownames(pheno_data) <- data01$pat_id
  pheno_data$age[is.na(pheno_data$age)] <- mean(pheno_data$age, na.rm = TRUE)
  pheno_data$sex[is.na(pheno_data$sex)] <- mode_function(pheno_data$sex)
  
  
  ### PCA after batch removal -------------------------------------------------
  runPCA_batch <- TRUE
  # PCA
  if(runPCA_batch == TRUE) {
    # pca
    pca_res <- prcomp(exprs_data, center = TRUE, scale. = TRUE)
    
    # Get the first two principal components
    pca_df <- as_tibble(as.data.frame(pca_res$x[,1:5]))
    
    exprs_metadat <- exprs_metadat %>% 
      bind_cols(pca_df)
    
    # PC-plot
    # # Custom labels for the legend
    # custom_labels <- c(
    #   "acs" = "ACS",
    #   "cad" = "CAD",
    #   "control" = "Control",
    #   "dcm" = "DCM",
    #   "ref" = "HFrEF"
    # )
    
    palette <- brewer.pal(12, "Set3")
    # takes some time
    ggplot(exprs_metadat, aes(PC1, PC2, color = batch_center)) +
      geom_point(alpha=0.9, aes(shape=disease)) +
      #coord_fixed(ratio=1)+
      labs(title = NULL,  # "PCA Plot Colored by Disease",
           x = "PC1",
           y = "PC2") +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette, name=NULL) + #, labels = custom_labels)+
      scale_shape_discrete(name=NULL) +
      my_base_theme() -> plot_PC1_PC2
    # plot_PC1_PC2
    
    ggplot(exprs_metadat, aes(PC1, PC3, color = batch_center)) +
      geom_point(alpha=0.9, aes(shape=disease)) +
      labs(title = NULL,  # "PCA Plot Colored by Disease",
           x = "PC1",
           y = "PC3") +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette, name=NULL) + #, labels = custom_labels)+
      scale_shape_discrete(name=NULL,) +
      my_base_theme() -> plot_PC1_PC3
    #plot_PC1_PC3
    
    if (SAVE.files ==TRUE) {
      pca12_beforebatch <- glue("{data_path_bestageing2022}/output/plots/preprocessing/pca_batch_corr/001c_ANOVA_pca12_before_batch.svg")
      ggsave(filename = pca12_beforebatch, plot = plot_PC1_PC2, 
             width = 8, height = 6, 
             units = "in"  # default
      )
      pca13_beforebatch <- glue("{data_path_bestageing2022}/output/plots/preprocessing/pca_batch_corr/001c_ANOVA_pca13_before_batch.svg")
      ggsave(filename = pca13_beforebatch, plot = plot_PC1_PC3, 
             width = 6, height = 6, 
             units = "in"  # default
      )
    }
  }
  
  # run limma batch remove https://evayiwenwang.github.io/Managing_batch_effects/adjust.html#accounting-for-batch-effects
  # modlimma <- model.matrix(~1, data=pheno_data) 
  modlimma <- model.matrix( ~ disease, data=pheno_data)
  adjustedMatrix <- removeBatchEffect(t(exprs_data), batch=batch_slide_id, design = modlimma)  # collinearity issues with center batch:  one-to-one mapping with any condition
  limma_edata <-  adjustedMatrix %>% t() %>%  as.data.frame() %>% as_tibble() %>% 
    mutate(pat_id = data01$pat_id) %>% 
    select(pat_id, everything())
  
  # inspect random boxplots 
  set.seed(123)
  random_index <- sample(x = 1:nrow(exprs_data), size = 50, replace = FALSE) %>% sort()  # samples, not miRNAs
  
  # fast boxplot inspect
  #par(mfrow=c(1,2))
  #boxplot(as.data.frame(t(exprs_data))[, random_index],main="Original")
  #boxplot(as.data.frame(adjustedMatrix)[, random_index],main="Batch corrected")
  
  # ggboxplot 
  df_original <- as.data.frame(t(exprs_data)) %>% 
    select(all_of(random_index)) %>% 
    tibble::rownames_to_column(var = "miRNA") %>% 
    pivot_longer(cols = -miRNA, names_to = "Sample", values_to = "Expression")
  
  df_corrected <- as.data.frame(adjustedMatrix) %>% 
    select(all_of(random_index)) %>% 
    tibble::rownames_to_column(var = "miRNA") %>% 
    pivot_longer(cols = -miRNA, names_to = "Sample", values_to = "Expression")
  
  # Find the global y limits to make them same for both plots
  y_limits <- c(min(c(df_original$Expression, df_corrected$Expression)), 
                max(c(df_original$Expression, df_corrected$Expression)))
  
  # Plot original data
  p1 <- ggplot(df_original, aes(x = Sample, y = Expression)) +
    geom_boxplot(aes(group = Sample), alpha=0.5) +
    theme_minimal(base_size = 16, base_family = 'Arial')+
    labs(title = NULL, x=NULL) +
    ylim(y_limits) +
    scale_color_manual(values = palette, name=NULL) +
    theme(axis.text.x = element_blank()) +
    my_base_theme() 
  # Plot batch-corrected data
  p2 <- ggplot(df_corrected, aes(x = Sample, y = Expression)) +
    geom_boxplot(aes(group = Sample), alpha=0.5) +
    theme_minimal(base_size = 16, base_family = 'Arial')+
    labs(title = NULL, x=NULL) +
    ylim(y_limits) +
    scale_color_manual(values = palette, name=NULL) +
    theme(axis.text.x = element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank()  # remove y labels from 2nd plot
    ) +
    my_base_theme() 
  
  # Combine plots
  combined_plot <- p1 | p2
  if (SAVE.files ==TRUE) {
    boxplot_batch_correct <- glue("{data_path_bestageing2022}/output/plots/preprocessing/boxplot_batch_corr/001c_ANOVA_boxplot_before_after_batch.svg")
    ggsave(filename = boxplot_batch_correct, plot = combined_plot, 
           width = 8, height = 6, 
           units = "in"  # default
    )
  }
  # pValuesLimma = f.pvalue(adjustedMatrix,design,null_model)
  # qValuesLimma = p.adjust(pValuesLimma,method="BH")
  # sum(qValuesLimma<0.05)
  
  runPCA_batch2 <- TRUE
  # PCA
  if(runPCA_batch2 == TRUE) {
    # pca
    pca_res <- prcomp(limma_edata[-1], center = TRUE, scale. = TRUE)
    
    # Get the first two principal components
    pca_df <- as_tibble(as.data.frame(pca_res$x[,1:5]))
    
    # add metadata
    exprs_metadat <- limma_edata %>% 
      bind_cols(pca_df) %>% 
      mutate(batch_center = batch_center, 
             disease = data01$disease)
    
    palette <- brewer.pal(12, "Set3")
    # takes some time
    ggplot(exprs_metadat, aes(PC1, PC2, color = batch_center)) +
      geom_point(alpha=0.9, aes(shape=disease)) +
      labs(title = NULL,  # "PCA Plot Colored by Disease",
           x = "PC1",
           y = "PC2") +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette, name=NULL) + #, labels = custom_labels)+
      scale_shape_discrete(name=NULL) +
      my_base_theme() -> plot_PC1_PC2_corrected
    # plot_PC1_PC2_corrected
    
    ggplot(exprs_metadat, aes(PC1, PC3, color = batch_center)) +
      geom_point(alpha=0.9, aes(shape=disease)) +
      labs(title = NULL,  # "PCA Plot Colored by Disease",
           x = "PC1",
           y = "PC3") +
      theme_minimal(base_size = 16, base_family = 'Arial')+
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette, name=NULL) + #, labels = custom_labels)+
      scale_shape_discrete(name=NULL) +
      my_base_theme() -> plot_PC1_PC3_corrected
    # plot_PC1_PC3_corrected
    
    if (SAVE.files ==TRUE) {
      pca12_afterbatch <- glue("{data_path_bestageing2022}/output/plots/preprocessing/pca_batch_corr/001c_ANOVA_pca12_after_batch.svg")
      ggsave(filename = pca12_afterbatch, plot = plot_PC1_PC2_corrected, 
             width = 8, height = 6, 
             units = "in"  # default
      )
      pca13_afterbatch <- glue("{data_path_bestageing2022}/output/plots/preprocessing/pca_batch_corr/001c_ANOVA_pca13_after_batch.svg")
      ggsave(filename = pca13_afterbatch, plot = plot_PC1_PC3_corrected, 
             width = 6, height = 6, 
             units = "in"  # default
      )
    }
  }
}

if(runLimmaRemoveBatchEffect == TRUE) {
  data01 <- data01 %>% select(1:4) %>% 
    left_join(limma_edata, by = join_by(pat_id))
}

## SVA ---------------------------------------------------------------------
if (runSVA == TRUE){
  exprs_data <- as.matrix(data01[,-c(1:4)])
  rownames(exprs_data) <- data01$pat_id
  pheno_data <- data01 %>% select(disease, age, sex) %>% as.data.frame()
  rownames(pheno_data) <- data01$pat_id
  
  # mean/ mode imputation necessary
  sum(is.na(pheno_data))
  pheno_data$age[is.na(pheno_data$age)] <- mean(pheno_data$age, na.rm = TRUE)
  
  #pheno_data$disease[is.na(pheno_data$disease)] <- mode_function(pheno_data$disease)
  pheno_data$sex[is.na(pheno_data$sex)] <- mode_function(pheno_data$sex)
  
  design <- model.matrix(~ disease + sex + age, data=pheno_data)
  null_model <- model.matrix(~ sex + age, data=pheno_data)
  
  #dim(exprs_data)
  #dim(pheno_data)
  #dim(design)
  #dim(null_model)
  
  n.sv <- num.sv(t(exprs_data), design, method="leek")  # ANOVA n.sv=0
  if (n.sv != 0 & exists(x = "n.sv")) {
    sva_obj <- sva(t(exprs_data), design, null_model, n.sv=n.sv)
    
    sva_var_df <- sva_obj$sv
    colnames(sva_var_df) <- paste0("sva_", seq(n.sv))
    
    # update design matrix
    sva_design <- cbind(design, sva_var_df)
    sva_null_model <- cbind(null_model, sva_var_df)
    
    # append sva to adjust later
    data01 <- cbind(data01, sva_var_df)  # new cols 
  }
}

# if (SAVE.files ==TRUE) {
#   filename.data01 <- glue("{data_path_bestageing2022}/output/de_results/{disease}/data01.rds")
#   saveRDS(object = data01, file = filename.data01)
# }

## SUBSET of literature miRNAs
# selected_mirna_dat <- paste0("data01_", length(researchMiRNAAccession$miRNAName_v21), "mirnas")  
# assign(x = selected_mirna_dat, 
#        value = data01 %>% 
#          select(pat_id, disease, age, sex, researchMiRNAAccession$miRNAName_v21))

## STORE PROCESSED DATA -------------------------------------------------------------------


data01$age[is.na(data01$age)] <- mean(data01$age, na.rm = TRUE)
data01$sex[is.na(data01$sex)] <- mode_function(data01$sex)

if (SAVE.files ==TRUE) {
  path2dataprocessed <- glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_ANOVA_data01.rds")
  saveRDS(object = data01, file = path2dataprocessed)
}

# RUN ANCOVA  ----------------------------------------------------
# create space to store results for ALL miRNAs
all_filtered_mirnas <- data01 %>% select(-c(disease, age, sex))  # bring to same str as "all_filtered_mirnas" to work with existing code

pval.anova.disease <-rep(NA,ncol(all_filtered_mirnas)-1)
name.mir <- rep(NA,ncol(all_filtered_mirnas)-1)

median.cont <- rep(NA,ncol(all_filtered_mirnas)-1)
median.acs <- rep(NA,ncol(all_filtered_mirnas)-1)
median.cad <- rep(NA,ncol(all_filtered_mirnas)-1)
median.dcm <- rep(NA,ncol(all_filtered_mirnas)-1)
median.ref <- rep(NA,ncol(all_filtered_mirnas)-1)


mean.cont <- rep(NA, ncol(all_filtered_mirnas) - 1)
mean.acs <- rep(NA, ncol(all_filtered_mirnas) - 1)
mean.cad <- rep(NA, ncol(all_filtered_mirnas) - 1)
mean.dcm <- rep(NA, ncol(all_filtered_mirnas) - 1)
mean.ref <- rep(NA, ncol(all_filtered_mirnas) - 1)

sd.cont <- rep(NA, ncol(all_filtered_mirnas) - 1)
sd.acs <- rep(NA, ncol(all_filtered_mirnas) - 1)
sd.cad <- rep(NA, ncol(all_filtered_mirnas) - 1)
sd.dcm <- rep(NA, ncol(all_filtered_mirnas) - 1)
sd.ref <- rep(NA, ncol(all_filtered_mirnas) - 1)

# post hoc analysis grouped table
emmeans_1_pairs <- tibble(
  miRNA= NA, 
  contrast = NA_character_,
  estimate = NA,
  SE = NA,
  df = NA,
  t.ratio=NA, 
  p.value = NA)
                          
# order factors for post-hoc analyis (control-acs, control -ref, etc---)
data01 <- data01 %>% 
  mutate(disease = factor(disease, levels = c("control", "acs", "cad", "dcm", "ref")))

library(emmeans)  # for post hoc analysis of ANCOVA: estimated marginal means (EMMs), also known as least-squares means. These are means adjusted for the covariates in your model.

if (runTests == TRUE){  #takes time
  total <- ncol(all_filtered_mirnas)-1
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  for(miRNA in 1:(ncol(all_filtered_mirnas)-1)) {
    # update index since first colnames are [1] "pat_id"  "disease"  "age"   "sex"  "hsa_let_7a_3p" 
    miRNA_col <- miRNA+4
    # cont <- data01[cont.index, miRNA_col] %>% as_vector()
    # case <- data01[case.index, miRNA_col] %>% as_vector()
    
    # calc medians for each disease group
    # e Table 3 Fehlmann 2020 included medians, means, and sd
    name.mir[miRNA] <- names(data01[miRNA_col])
    
    median.cont[miRNA] <- median(data01 %>% filter(disease == "control") %>% pull(miRNA_col))
    median.acs[miRNA] <- median(data01 %>% filter(disease == "acs") %>% pull(miRNA_col))
    median.cad[miRNA] <- median(data01 %>% filter(disease == "cad") %>% pull(miRNA_col))
    median.dcm[miRNA] <- median(data01 %>% filter(disease == "dcm") %>% pull(miRNA_col))
    median.ref[miRNA] <- median(data01 %>% filter(disease == "ref") %>% pull(miRNA_col))
    
    mean.cont[miRNA] <- mean(data01 %>% filter(disease == "control") %>% pull(miRNA_col))
    mean.acs[miRNA] <- mean(data01 %>% filter(disease == "acs") %>% pull(miRNA_col))
    mean.cad[miRNA] <- mean(data01 %>% filter(disease == "cad") %>% pull(miRNA_col))
    mean.dcm[miRNA] <- mean(data01 %>% filter(disease == "dcm") %>% pull(miRNA_col))
    mean.ref[miRNA] <- mean(data01 %>% filter(disease == "ref") %>% pull(miRNA_col))
    
    sd.cont[miRNA] <- sd(data01 %>% filter(disease == "control") %>% pull(miRNA_col))
    sd.acs[miRNA] <- sd(data01 %>% filter(disease == "acs") %>% pull(miRNA_col))
    sd.cad[miRNA] <- sd(data01 %>% filter(disease == "cad") %>% pull(miRNA_col))
    sd.dcm[miRNA] <- sd(data01 %>% filter(disease == "dcm") %>% pull(miRNA_col))
    sd.ref[miRNA] <- sd(data01 %>% filter(disease == "ref") %>% pull(miRNA_col))
    
    # run ANOVA
    f <- as.formula(paste(name.mir[miRNA], " ~ disease + sex + age"))
    aov_result_1 <- aov(f, data = data01)
    aov_result_1_summary <- summary(aov_result_1)
    
    pval.anova.disease[miRNA] <- aov_result_1_summary[[1]]["disease", "Pr(>F)"]
    
    # perform post hoc test if anova shows significant differences, usually use Tukey HSD (Honestly Significant Difference) test, but here we used ANCOVA
    emmeans_1 <- emmeans(aov_result_1, ~ disease)
    emmeans_1_pairs_tmp <- pairs(emmeans_1)  # P value adjustment: tukey method for comparing a family of 5 estimates 
    emmeans_1_pairs_tmp <- emmeans_1_pairs_tmp %>% as_tibble() %>% 
      mutate(miRNA = name.mir[miRNA])
    
    emmeans_1_pairs <- emmeans_1_pairs %>% 
      bind_rows(emmeans_1_pairs_tmp)
    # status bar
    setTxtProgressBar(pb, miRNA)
  }
  
  # GATHER Results -----------------------------------------------------------
  # we conducted 2549 t-tests and 2549 glm-models for each gene
  anova.results <- tibble(miRNA = colnames(all_filtered_mirnas)[-1],  # all miRNAs without patID
                       # average.difference = average.difference,
                       pval.anova.disease = pval.anova.disease,
                       median.cont = median.cont,
                       median.acs = median.acs,
                       median.cad = median.cad,
                       median.dcm = median.dcm,
                       median.ref = median.ref,
                       
                       mean.cont = mean.cont,
                       mean.acs = mean.acs,
                       mean.cad = mean.cad,
                       mean.dcm = mean.dcm,
                       mean.ref = mean.ref,
                       
                       sd.cont = sd.cont,
                       sd.acs = sd.acs,
                       sd.cad = sd.cad,
                       sd.dcm = sd.dcm,
                       sd.ref = sd.ref)
  
  # for figure 1 vogel 2013 
  adjusting_method <- "holm" # "holm" # "BH" = Benjamini Hochberg
  anova.results <- anova.results %>% 
    mutate(padj.anova = p.adjust(pval.anova.disease, method = adjusting_method, n = length(pval.anova.disease)) )  %>% 
    mutate(sign_indicator = ifelse(padj.anova < 0.05, "p.adj≤0.05", "n.s.")) %>% 
    relocate(padj.anova, .after = pval.anova.disease) %>%
    arrange(pval.anova.disease)

  
  emmeans_1_pairs <- emmeans_1_pairs[-1, ]
  
  if (SAVE.files ==TRUE) {
    filename.anova.results <- glue("{data_path_bestageing2022}/output/de_results/anova/001c_anova_results_batch_corrected.rds")
    saveRDS(object = anova.results, file = filename.anova.results)
    filename.emmeans.results <- glue("{data_path_bestageing2022}/output/de_results/anova/001c_emmeans_results_batch_corrected.rds")
    saveRDS(object = emmeans_1_pairs, file = filename.emmeans.results)
  }
}




anova.results.xlsx <- anova.results %>% 
  select(-sign_indicator) %>% 
  rename(anova_rawp = pval.anova.disease, anova_adjp = padj.anova) %>%
  mutate(miRNA = convert_mir_name_V(miRNA))

filename.anova.xlsx <- glue("{data_path_bestageing2022}/output/de_results/anova/001c_anova_results_batch_corrected.xlsx")
write.xlsx(anova.results.xlsx, file = filename.anova.xlsx)


# read?
anova.results <- readRDS(glue("{data_path_bestageing2022}/output/de_results/anova/001c_anova_results_batch_corrected.rds"))
emmeans_1_pairs <- readRDS(glue("{data_path_bestageing2022}/output/de_results/anova/001c_emmeans_results_batch_corrected.rds"))

# top 5 dysregulated miRNAs
anova.results_top5 <- anova.results %>% 
  slice(1:5)
top5_anova_mirnas <- anova.results_top5 %>% pull(miRNA)

emmeans_1_pairs_anova_top5 <- emmeans_1_pairs %>% 
  filter(miRNA %in% top5_anova_mirnas) %>% 
  filter(p.value < 0.05) %>% 
  arrange(match(miRNA, top5_anova_mirnas))  %>% #  arranges the rows in the order specified by the top5_anova_mirnas vector
  mutate(miRNA = convert_mir_name_V(miRNA))

# Word export ----------------------------------------------------------------
# https://cran.r-project.org/web/packages/gtsummary/vignettes/rmarkdown.html#:~:text=%7Bflextable%7D%20is%20the%20default%20print,table%20printed%20with%20%7Bgt%7D.
library(flextable) # {flextable} is the default print engine for Word output, as {gt} does not support Word. If {flextable} is not installed, kable is used.
library(officer)

anova.results.xlsx$anova_rawp <- sprintf("%.2E", anova.results.xlsx$anova_rawp)  # converts to string!
anova.results.xlsx$anova_adjp <- sprintf("%.2E", anova.results.xlsx$anova_adjp)

# for word export
sect_properties <- prop_section(
  page_size = page_size(
    orient = "landscape",
    width = 8.3, height = 14
  ),
  type = "continuous",
  page_margins = page_mar()
)

# table 1
file_path <- glue("{data_path_bestageing2022}/output/de_results/anova/001c_anova_results_batch_corrected.docx")
anova_flextable <- flextable(anova.results.xlsx) %>% 
  colformat_double(
    big.mark = ",", digits = 2, na_str = "N/A"
  ) %>% 
  set_table_properties(layout = "autofit", align= "left") %>%
  fontsize(size = 8, part = "all") %>% 
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  autofit() %>% 
  save_as_docx(path=file_path, pr_section = sect_properties)


# table 2 emmeans top
emmeans_1_pairs_anova_top5.word <- emmeans_1_pairs_anova_top5
emmeans_1_pairs_anova_top5.word$p.value <- sprintf("%.2E", emmeans_1_pairs_anova_top5.word$p.value) 

file_path_emmeans <- glue("{data_path_bestageing2022}/output/de_results/anova/001c_anova_emmeansTOP_results_batch_corrected.docx")

emmeans_1_pairs_anova_top5.word_dcm_ref <- emmeans_1_pairs_anova_top5.word %>% 
  group_by(contrast) %>% 
  filter( contrast == "control - dcm" | contrast == "control - ref") %>% 
  ungroup()

anova_flextable <- flextable(emmeans_1_pairs_anova_top5.word) %>% 
  colformat_double(
    big.mark = ",", digits = 3, na_str = "N/A"
  ) %>% 
  set_table_properties(layout = "autofit") %>%
  fontsize(size = 11, part = "all") %>% 
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  autofit() %>% 
  save_as_docx(path=file_path_emmeans, pr_section = sect_properties)

# text 
# mirnas
glue("The top 5 dysregulated miRNAs were {top5_anova_mirnas[1]}, {top5_anova_mirnas[2]}, {top5_anova_mirnas[3]}, {top5_anova_mirnas[4]}, and {top5_anova_mirnas[5]}.")
# mean [sd]
glue("Mean expression for {top5_anova_mirnas[1]} [standard deviaton] of control {round(anova.results.xlsx$mean.cont[1],3)} [{round(anova.results.xlsx$sd.cont[1],3)}] vs. DCM {round(anova.results.xlsx$mean.dcm[1],3)} [{round(anova.results.xlsx$sd.dcm[1],3)}] (adj.p-val={emmeans_1_pairs_anova_top5.word_dcm_ref$p.value[1]}) vs. HFrEF {round(anova.results.xlsx$mean.ref[1],3)} [{round(anova.results.xlsx$sd.ref[1],3)}] (adj. p-val={emmeans_1_pairs_anova_top5.word_dcm_ref$p.value[2]})")
glue("Mean expression [standard deviaton] of control {round(anova.results.xlsx$mean.cont[2],3)} [{round(anova.results.xlsx$sd.cont[2],3)}] vs. DCM {round(anova.results.xlsx$mean.dcm[2],3)} [{round(anova.results.xlsx$sd.dcm[2],3)}] (adj. p-val={emmeans_1_pairs_anova_top5.word_dcm_ref$p.value[3]}) vs. HFrEF {round(anova.results.xlsx$mean.ref[2],3)} [{round(anova.results.xlsx$sd.ref[2],3)}] (adj. p-val={emmeans_1_pairs_anova_top5.word_dcm_ref$p.value[4]})")

for (i in 1:5) {
  expression_string <- glue("Mean expression for {convert_mir_name_V(top5_anova_mirnas[i])} [standard deviation] of control {round(anova.results.xlsx$mean.cont[i],3)} [{round(anova.results.xlsx$sd.cont[i],3)}] vs. DCM {round(anova.results.xlsx$mean.dcm[i],3)} [{round(anova.results.xlsx$sd.dcm[i],3)}] (adj.p-val={emmeans_1_pairs_anova_top5.word_dcm_ref$p.value[2*i-1]}) vs. HFrEF {round(anova.results.xlsx$mean.ref[i],3)} [{round(anova.results.xlsx$sd.ref[i],3)}] (adj. p-val={emmeans_1_pairs_anova_top5.word_dcm_ref$p.value[2*i]})")
  print(expression_string)
}

# p-value
glue("p-value: {anova.results.xlsx$anova_rawp[1]}")
emmeans_1_pairs_anova_top5.word_dcm_ref$p.value[1]
