


# Create hierarchical clustering plot and heatmap

# complex heatmaps
library(ComplexHeatmap)  # https://github.com/jokergoo/ComplexHeatmap
# library(pheatmap)  # https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
library(glue)
library(dplyr)
library(circlize)

data_path_bestageing2022 <- "/mnt/users/reich/rockerprojects/bestageing2022"
data_path_BestAgeing <- "/mnt/users/reich/BestAgeing"


# load data ---------------------
path2dataprocessed_acs <- glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_acs_data01.rds")
path2dataprocessed_cad <- glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_cad_data01.rds")
path2dataprocessed_dcm <- glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_dcm_data01.rds")
path2dataprocessed_hfref <- glue("{data_path_bestageing2022}/data/Rdata/processed_disease_data/001c_hfref_data01.rds")

data01_acs <- readRDS(file = path2dataprocessed_acs) %>% as_tibble()
data01_cad <- readRDS(file = path2dataprocessed_cad) %>% as_tibble()
data01_dcm <- readRDS(file = path2dataprocessed_dcm) %>% as_tibble()
data01_hfref <- readRDS(file = path2dataprocessed_hfref) %>% as_tibble()

reshape_dat <- data01_acs %>% 
  dplyr::bind_rows(data01_cad) %>% 
  dplyr::bind_rows(data01_dcm) %>% 
  dplyr::bind_rows(data01_hfref) %>% 
  dplyr::distinct(pat_id, .keep_all = TRUE) 


# UPDATE 20240125  miRetrieve miRs
table_mirna_top50bm_score_alldiseases <- readRDS(file = glue("{data_path_bestageing2022}/data-literature/miRetrieve/20240125top50mirnas_all_diseases_pmids_gpt.rds"))

# disease specific
# table_mirna_top50bm_score_alldiseases <- table_mirna_top50bm_score_alldiseases %>% 
#   filter(toupper(Topic) == toupper(all_combis$diseases[i]))

# select cardiovascular miRNAs
set.seed(1234)
reshape_dat_slice <- reshape_dat %>% 
  select(pat_id, disease, age, sex, any_of(table_mirna_top50bm_score_alldiseases$TargetName)) %>% 
  group_by(disease) %>% 
  # Sampling 100 patients from each disease group
  # If a disease group has less than 100 patients, sample that many without replacement
  # group_modify(~sample_n(., min(100, n()), replace = FALSE))  %>%  # better than slice_sample(n = 100, replace = FALSE) which wont throw an error if disease group has less than n patients
  slice_sample(n = 200, replace = FALSE) %>% 
  ungroup()

reshape_dat_mirna <- reshape_dat_slice %>% 
  select(-c(disease, age, sex))

sum(is.na(reshape_dat_mirna))  # missing due to filtering ;)

transposed_df <- t(reshape_dat_mirna[,-1])
rownames(transposed_df) <- colnames(reshape_dat_mirna)[-1]
rownames(transposed_df) <- gsub("hsa_mir_", "miR-", rownames(transposed_df))  # Replace 'hsa_mir_' with 'miR-'
rownames(transposed_df) <- gsub("_", "-", rownames(transposed_df))

colnames(transposed_df) <- reshape_dat_mirna$pat_id

#ComplexHeatmap::pheatmap(transposed_df)

# normalize rowwise (scale outputs columns -> need to transpose again)
# normalized_matrix <- t(apply(transposed_df, 1, scale))  # done by pheatmap :)

annotation_col = data.frame(
  Disease = factor(reshape_dat_slice$disease,
                   levels = c("acs", "cad", "dcm", "hfref", "control"),
                   labels = c("ACS", "CAD", "DCM", "ICM", "Control"))
)

# thematic::okabe_ito(6)
ann_colors = list(
  Disease = c(ACS = "#CC79A7", CAD = "#999999", DCM = "#0072B2", ICM = "#E69F00", Control = "#009E73")
)

pheatmap1 <- ComplexHeatmap::pheatmap(
  transposed_df, 
  name="Row Z-score",
  col = circlize::colorRamp2(c(-5, 0, 5), c("#0072B2", "white", "#E69F00")),
  scale="row",
  clustering_distance_rows="euclidean",
  clustering_method="complete", border_color=FALSE,
  annotation_col = annotation_col, # annotation_row = annotation_row, 
  annotation_colors = ann_colors,
  column_split = annotation_col$Disease,
  show_colnames = FALSE,
  heatmap_legend_param = list(
    legend_direction = "horizontal", 
    legend_width = unit(6, "cm"))
)
# pheatmap1
# https://stackoverflow.com/questions/68483646/complexheatmap-how-to-place-heatmap-legend-and-annotation-legend-differently
draw(
  pheatmap1, heatmap_legend_side="bottom", annotation_legend_side="right",
  legend_grouping = "original"
)

pheatmap2 <- ComplexHeatmap::pheatmap(
  transposed_df, 
  name="Row Z-score",
  col = circlize::colorRamp2(c(-5, 0, 5), c("#0072B2", "white", "#E69F00")),
  scale="row",
  clustering_distance_rows="euclidean",
  clustering_method="complete", border_color=FALSE,
  annotation_col = annotation_col, # annotation_row = annotation_row, 
  annotation_colors = ann_colors,
  #column_split = annotation_col$Disease,
  show_colnames = FALSE,
  heatmap_legend_param = list(
    legend_direction = "horizontal", 
    legend_width = unit(6, "cm"))
)
pheatmap2

# https://stackoverflow.com/questions/68483646/complexheatmap-how-to-place-heatmap-legend-and-annotation-legend-differently
draw(
  pheatmap2, heatmap_legend_side="bottom", annotation_legend_side="right",
     legend_grouping = "original"
  )



# pairwise ----------------

# Define disease groups for pairwise comparison
disease_groups <- unique(reshape_dat$disease)

# Loop through each disease group and compare with control
set.seed(123)
for(disease in disease_groups) {
  if(disease != "control") {
    # Subset for disease and control
    subset_dat <- reshape_dat %>%
      filter(disease %in% c(!!disease, "control")) %>%
      select(pat_id, disease, any_of(table_mirna_top50bm_score_alldiseases$TargetName)) %>%
      group_by(disease) %>% 
      slice_sample(n = 50, replace = FALSE) %>%
      ungroup()
    
    # Prepare data for heatmap
    transposed_df_loop <- t(subset_dat[,-c(1:2)])
    rownames(transposed_df_loop) <- colnames(subset_dat)[-c(1:2)]
    rownames(transposed_df_loop) <- gsub("hsa_mir_", "miR-", rownames(transposed_df_loop))
    rownames(transposed_df_loop) <- gsub("_", "-", rownames(transposed_df_loop))
    colnames(transposed_df_loop) <- subset_dat$pat_id
    
    
    # Remove rows with NA values, otherwise error in calculating distance
    transposed_df_loop <- transposed_df_loop[!rowSums(is.na(transposed_df_loop)), ]
    
    
    # Prepare annotations
    annotation_col <- data.frame(
      Disease = factor(subset_dat$disease, levels = c(disease, "control"), labels = c(toupper(disease), "Control"))
    )
    
    
    ann_colors <- list(
      Disease = setNames(c("#E69F00", "#009E73"), c(toupper(disease), "Control"))
    )
    
    # Generate heatmap
    heatmap_legend_param <- list(legend_direction = "horizontal", legend_width = unit(6, "cm"))
    
    pheatmap1_loop <- ComplexHeatmap::pheatmap(
      transposed_df_loop,
      name = "Row Z-score",
      col = circlize::colorRamp2(c(-2, 0, 2), c("#0072B2", "white", "#E69F00")),
      scale = "row",
      clustering_distance_rows = "euclidean",
      clustering_method = "complete",
      border_color = FALSE,
      annotation_col = annotation_col,
      annotation_colors = ann_colors,
      column_split = annotation_col$Disease,
      show_colnames = FALSE,
      heatmap_legend_param = heatmap_legend_param
    )
    
    # Draw heatmap with legend at the bottom
    draw(pheatmap1_loop, heatmap_legend_side = "bottom", annotation_legend_side = "right", legend_grouping = "original")
    
    
    pheatmap2_loop <- ComplexHeatmap::pheatmap(
      transposed_df_loop,
      name = "Row Z-score",
      col = circlize::colorRamp2(c(-2, 0, 2), c("#0072B2", "white", "#E69F00")),
      scale = "row",
      clustering_distance_rows = "euclidean",
      clustering_method = "complete",
      border_color = FALSE,
      annotation_col = annotation_col,
      annotation_colors = ann_colors,
      show_colnames = FALSE,
      heatmap_legend_param = heatmap_legend_param
    )
    
    # Draw heatmap with legend at the bottom
    draw(pheatmap2_loop, heatmap_legend_side = "bottom", annotation_legend_side = "right", legend_grouping = "original")
    
  }
}
