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
  
  # Layout: 2 rows, image and caption
  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(0.8, 0.2), "npc"))))
  
  # Draw image
  grid.raster(
    img,
    vp = viewport(layout.pos.row = 1, layout.pos.col = 1),
    width = unit(scale_width, "npc"),
    height = unit(scale_width, "npc") * (nrow(img) / ncol(img))
  )
  
  # Extract prefix (e.g., "Figure 1:") and body
  full_caption <- captions[i]
  parts <- strsplit(full_caption, ":\\s*", perl = TRUE)[[1]]
  
  if (length(parts) >= 2) {
    prefix <- paste0(parts[1], ":")
    body <- paste(parts[-1], collapse = ": ")
  } else {
    prefix <- full_caption
    body <- ""
  }
  
  # Draw bold prefix (Figure X:)
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid.text(prefix,
            x = unit(0.05, "npc"), y = unit(0.95, "npc"),
            just = c("left", "top"),
            gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Times"))
  
  # Draw wrapped regular caption below
  wrapped_body <- str_wrap(body, width = 140)
  grid.text(wrapped_body,
            x = unit(0.05, "npc"), y = unit(0.85, "npc"),
            just = c("left", "top"),
            gp = gpar(fontsize = 8, fontface = "plain", fontfamily = "Times"))
  
  popViewport()
}

dev.off()
