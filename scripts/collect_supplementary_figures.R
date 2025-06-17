# Created 17/06/2025 by Jennifer Astley <jennifer.astley@kennedy.ox.ac.uk>
# Last modified 17/06/2025 by Jennifer Astley

library(grid)
library(gridExtra)
library(png)
library(readr)
library(stringr)

# Set paths
figure_dir <- "/well/combat/users/bsg751/projects/eqtl/runtensorqtl/COMBAT_eQTL/figures/supplementary/figures"
caption_file <- "/well/combat/users/bsg751/projects/eqtl/runtensorqtl/COMBAT_eQTL/figures/supplementary/figure_data/figure_captions.csv"
output_pdf <- "/well/combat/users/bsg751/projects/eqtl/runtensorqtl/COMBAT_eQTL/figures/supplementary/supplementary_figures.pdf"

# Load captions
caption_data <- read_csv(caption_file)

# List and sort figure files
fig_files <- list.files(figure_dir, pattern = "^figure_S[0-9]+\\.png$", full.names = TRUE)
fig_files_sorted <- fig_files[order(as.numeric(str_extract(basename(fig_files), "\\d+")))]

# Extract captions for sorted files
captions <- sapply(basename(fig_files_sorted), function(f) {
  match <- caption_data$caption[caption_data$filename == f]
  if (length(match) == 0) {
    paste("Caption missing for", f)
  } else {
    match
  }
})

# Custom width scale factors (1 = full width)
custom_widths <- c(
  #"figure_S3.png" = 0.7,
  #"figure_S7.png" = 0.6,
  "figure_S8.png" = 0.9,
  "figure_S9.png" = 0.8,
  "figure_S10.png" = 0.8
)

# Open PDF device with Times font
pdf(output_pdf, width = 8, height = 10, family = "Times")

for (i in seq_along(fig_files_sorted)) {
  img <- readPNG(fig_files_sorted[i])
  fig_name <- basename(fig_files_sorted[i])
  
  # Determine width scale factor
  scale_width <- ifelse(fig_name %in% names(custom_widths),
                        custom_widths[fig_name], 1)
  
  grid.newpage()
  
  # Layout: 2 rows, image and caption (caption directly below image)
  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(0.8, 0.2), "npc"))))
  
  # Draw image scaled by width, maintain aspect ratio, centered
  grid.raster(img,
              vp = viewport(layout.pos.row = 1, layout.pos.col = 1),
              width = unit(scale_width, "npc"),
              height = unit(scale_width, "npc") * (nrow(img) / ncol(img)))
  
  # Draw caption just below image, with text wrapping
  grid.text(str_wrap(captions[i], width = 100),
            vp = viewport(layout.pos.row = 2, layout.pos.col = 1),
            gp = gpar(fontsize = 12, fontfamily = "Times"),
            just = "left", x = unit(0.05, "npc"))
}

dev.off()
