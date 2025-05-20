# Custom ggplot theme for BestAgeing2022 Project
# Author: Christoph Reich
# Date: 2025-05-20

# Load required packages
if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require(thematic)) {
  install.packages("thematic")
  library(thematic)
}

# Define custom base theme for consistent visualization
my_base_theme <- function(base_size = 12, base_family = "Arial") {
  # Create a clean, publication-ready theme
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Text elements
      plot.title = element_text(
        size = rel(1.3),
        hjust = 0.5,
        face = "bold",
        margin = margin(t = 5, b = 15)
      ),
      plot.subtitle = element_text(
        size = rel(1.1),
        hjust = 0.5,
        margin = margin(t = 0, b = 10)
      ),
      plot.caption = element_text(
        size = rel(0.8),
        hjust = 1,
        margin = margin(t = 10)
      ),

      # Axis elements
      axis.title = element_text(face = "bold", size = rel(1.1)),
      axis.text = element_text(size = rel(0.9), color = "black"),
      axis.ticks = element_line(color = "black"),

      # Panel elements
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),

      # Legend elements
      legend.title = element_text(face = "bold", size = rel(1)),
      legend.text = element_text(size = rel(0.9)),
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = NA),

      # Facet elements
      strip.background = element_rect(fill = "gray95", color = "black"),
      strip.text = element_text(face = "bold", size = rel(0.9)),

      # Remove clipping to allow drawing outside of plot area
      plot.margin = margin(15, 15, 10, 10)
    )
}

# Define color palettes for consistent use across plots
bestageing_colors <- list(
  # Disease colors
  disease = c(
    "control" = "#999999", # Gray
    "dcm" = "#E69F00", # Orange
    "acs" = "#56B4E9", # Blue
    "cad" = "#009E73", # Green
    "hfref" = "#D55E00" # Red
  ),

  # Regulation colors for volcano plots
  regulation = c(
    "Upregulated" = "#D55E00", # Red
    "Downregulated" = "#56B4E9", # Blue
    "Not significant" = "#999999" # Gray
  ),

  # Continuous color scales
  heatmap_diverging = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"),
  heatmap_sequential = c("#f7fcfd", "#e0ecf4", "#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#6e016b")
)

# Function to apply brand colors to a ggplot
apply_bestageing_colors <- function(plot, color_type = "disease") {
  if (color_type == "disease") {
    plot + scale_color_manual(values = bestageing_colors$disease) +
      scale_fill_manual(values = bestageing_colors$disease)
  } else if (color_type == "regulation") {
    plot + scale_color_manual(values = bestageing_colors$regulation) +
      scale_fill_manual(values = bestageing_colors$regulation)
  } else {
    warning("Unknown color_type. Options are 'disease' or 'regulation'")
    plot
  }
}

# Function to save plots with consistent settings
save_bestageing_plot <- function(plot, filename, width = 8, height = 6, dpi = 300,
                                 units = "in", bg = "white", device = NULL) {
  if (is.null(device)) {
    # Determine device from file extension
    ext <- tolower(tools::file_ext(filename))
    device <- switch(ext,
      "pdf" = "pdf",
      "png" = "png",
      "jpg" = "jpeg",
      "jpeg" = "jpeg",
      "svg" = "svg",
      "tiff" = "tiff",
      "pdf"
    ) # Default to PDF if extension not recognized
  }

  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    units = units,
    bg = bg,
    device = device
  )
}
