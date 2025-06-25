library(grid)
library(png)
library(readr)
library(stringr)
library(tidyverse)

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

# Custom width scale factors
custom_widths <- c(
  "figure_S2.png" = 0.9,
  "figure_S8.png" = 0.9,
  "figure_S9.png" = 0.8,
  "figure_S10.png" = 0.8
)

# Helper function: parse markdown and wrap text with styles, then draw on grid
draw_marked_text_wrapped <- function(text, x_start = 0.05, y_start = 0.95, 
                                     max_width = 0.9, line_height = 5.5, fontsize = 12) {
  # Split text into tokens with style info:
  # **bold**, *italic*, or plain text chunks
  pattern <- "(\\*\\*[^*]+\\*\\*|\\*[^*]+\\*|[^*]+)"
  raw_chunks <- str_extract_all(text, pattern)[[1]]
  
  # Clean chunks and detect styles
  chunks <- lapply(raw_chunks, function(chunk) {
    if (startsWith(chunk, "**") && endsWith(chunk, "**")) {
      list(text = substr(chunk, 3, nchar(chunk) - 2), style = "bold")
    } else if (startsWith(chunk, "*") && endsWith(chunk, "*")) {
      list(text = substr(chunk, 2, nchar(chunk) - 1), style = "italic")
    } else {
      list(text = chunk, style = "plain")
    }
  })
  
  # Further split chunks by space to wrap at word level
  words <- list()
  for (ch in chunks) {
    ws <- strsplit(ch$text, " ")[[1]]
    for (w in ws) {
      words[[length(words) + 1]] <- list(text = w, style = ch$style)
    }
  }
  
  # Start drawing
  cur_x <- x_start
  cur_y <- y_start
  space_width <- unit(0.005, "npc") # small space between words
  fontsize_points <- fontsize
  fontfamily <- "Times"
  
  # We need to measure text widths to wrap
  # grid::stringWidth returns 'unit' object — convert to numeric npc units
  # Define a helper function:
  measure_text <- function(txt, face) {
    convertWidth(stringWidth(txt), "npc", valueOnly = TRUE)
  }
  
  for (i in seq_along(words)) {
    word <- words[[i]]$text
    face <- switch(words[[i]]$style,
                   "bold" = "bold",
                   "italic" = "italic",
                   "plain" = "plain")
    
    word_width <- measure_text(word, face)
    space_w <- measure_text(" ", "plain")
    
    # Check if word fits on current line, else wrap
    if ((cur_x + word_width) > (x_start + max_width)) {
      # Move to next line
      cur_x <- x_start
      cur_y <- cur_y - (line_height * fontsize_points / 72) * convertHeight(unit(1,"line"), "npc", valueOnly=TRUE)
    }
    
    # Draw the word
    grid.text(word, x = unit(cur_x, "npc"), y = unit(cur_y, "npc"),
              just = "left", 
              gp = gpar(fontsize = fontsize_points, fontface = face, fontfamily = fontfamily))
    
    # Advance x by word width + space
    cur_x <- cur_x + word_width + space_w
  }
}

# Open PDF device
pdf(output_pdf, width = 8, height = 10, family = "Times")

for (i in seq_along(fig_files_sorted)) {
  img <- readPNG(fig_files_sorted[i])
  fig_name <- basename(fig_files_sorted[i])
  
  # Determine width scale factor
  scale_width <- ifelse(fig_name %in% names(custom_widths),
                        custom_widths[fig_name], 1)
  
  grid.newpage()
  
  # Layout
  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(0.8, 0.2), "npc"))))
  
  # Draw image
  grid.raster(
    img,
    vp = viewport(layout.pos.row = 1, layout.pos.col = 1),
    width = unit(scale_width, "npc"),
    height = unit(scale_width, "npc") * (nrow(img) / ncol(img))
  )
  
  # Draw caption with markdown parsing & wrapping
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  draw_marked_text_wrapped(captions[i], x_start = 0.05, y_start = 0.95, max_width = 0.9, fontsize = 12)
  popViewport()
}

dev.off()

