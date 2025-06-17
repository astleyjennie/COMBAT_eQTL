# Created 17/06/2025 by Jennifer Astley <jennifer.astley@kennedy.ox.ac.uk>
# Last modified 17/06/2025 by Jennifer Astley

library(grid)
library(gridExtra)
library(png)
library(readr)
library(stringr)

# Set paths
figure_dir <- "/well/combat/users/bsg751/projects/eqtl/runtensorqtl/COMBAT_eQTL/figures/main/figures"
caption_file <- "/well/combat/users/bsg751/projects/eqtl/runtensorqtl/COMBAT_eQTL/figures/main/figure_data/main_figure_captions.csv"
output_pdf <- "/well/combat/users/bsg751/projects/eqtl/runtensorqtl/COMBAT_eQTL/figures/main/main_figures.pdf"
library(stringr)

# Load captions
caption_data <- read_csv(caption_file)

# List and sort figure files
fig_files <- list.files(figure_dir, pattern = "^Figure[1-4]+\\.png$", full.names = TRUE)
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
  "Figure1.png" = 0.8,
  "Figure2.png" = 0.9,
  "Figure3.png" = 0.7,
  "Figure4.png" = 0.8
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
  
  # Wrap caption text to ~60 characters width for better fitting
  wrapped_caption <- str_wrap(captions[i], width = 140)
  
  # Draw caption just below image, left-justified and indented a bit
  grid.text(
    wrapped_caption,
    vp = viewport(layout.pos.row = 2, layout.pos.col = 1),
    gp = gpar(fontsize = 8, fontfamily = "Times"),
    just = "left",
    x = unit(0.05, "npc")
  )
}

dev.off()
