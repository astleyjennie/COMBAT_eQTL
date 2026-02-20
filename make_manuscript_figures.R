# Created 04/12/2025 by Jennifer Astley
# Last modified 16/12/2025 by Jennifer Astley

# Setup ####

library(pacman)
p_load(tidyverse, png, data.table, patchwork, grid, UpSetR, cowplot, scales, viridis, colorspace)

setwd('/gpfs3/well/combat/users/bsg751/projects/eqtl/runtensorqtl/october_2025/COMBAT_eQTL')
# Set chosen working directory


# Colour schemes
celltype_colors <- c(
  "CD4" = "#6657EB",
  "cMono" = "#2EC4B6",
  "CD8" = "#FF6B6B",
  "NK" = "#118AB2",
  "B" = "#EF476F",
  "ncMono" = "#FFD166",
  "Other" = "#9C96D9"
)

# ============================
#        MAIN FIGURE 1       #
# ============================

# Made in Powerpoint. Figure stored in figures/main

# ============================
#        MAIN FIGURE 2       #
# ============================


# Figure 2a ####
# Bar plot showing the number of eQTLs found per cell type

cell_count <- fread('data/main/fig2a.csv') %>%
  mutate(cell_type = factor(cell_type, levels = cell_type))

fig2a <- ggplot(cell_count, aes(x = cell_type,
                                    y = num_cis_eQTLs,
                                    fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7, colour = "white") +
  geom_text(aes(label = num_cis_eQTLs), vjust = -0.5, size = 4) +
  scale_fill_manual(values = celltype_colors) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "",
    y = "Number of Independent cis-eQTLs",
    title = "Independent cis-eQTLs per Cell Type"
  )


# Figure 2b ####
# UMAP showing cell types

umap_count <- fread("data/main/fig2b.csv")

fig2b <- ggplot(umap_count, aes(x = UMAP1,
                                y = UMAP2,
                                colour = plot_cell_count)) + 
    geom_point() + 
    scale_color_manual(
        values = setNames(umap_count$plot_color, umap_count$plot_cell_count),
        breaks = c("CD4 (212443)", "cMono (164928)", "CD8 (81638)",
                   "NK (51982)", "B (33420)", "ncMono (26831)", "Other (36138)")
    ) + 
    theme_classic(base_size = 14) + 
    theme(
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right"
    ) +
    labs(
        x = "UMAP1",
        y = "UMAP2",
        title = "UMAP Plot of Major Cell Clusters",
        colour = NULL
    ) +
    guides(colour = guide_legend(override.aes = list(size = 6)))


inset_df <- data.frame(
  cell_type = c("CD4","cMono","CD8","NK","B","ncMono","Other"),
  cells_assayed = c(212443,164928,81638,51982,33420,26831,36138),
  eqtl_count    = c(1149,894,487,388,209,173,262),
  plot_col = c("#6657EB", "#2EC4B6", "#FF6B6B", "#118AB2", "#EF476F", "#FFD166", "#9C96D9")
)

inset_plot <- ggplot(inset_df, aes(x = cells_assayed, y = eqtl_count, colour = plot_col)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_colour_identity() +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    plot.margin = margin(2,2,2,2),
    legend.position = "none"
  ) +
  labs(x = "Cells assayed", y = "eQTLs")

get_legend <- function(p) {
  g <- ggplotGrob(p)
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}

legend_b <- get_legend(fig2b)
legend_plot <- wrap_elements(full = legend_b)   # <-- FIXED

legend_column <- inset_plot / legend_plot +
  plot_layout(heights = c(1, 1.3))

fig2b <- fig2b + theme(legend.position = "none")

final_fig2b <- fig2b + legend_column +
  plot_layout(widths = c(3, 1))



# Figure 2c ####
# eQTL plots of NAPSA, ZGLP1 and KANSL1

eqtl_box <- fread('data/main/fig2c.csv')

label_df <- eqtl_box %>%
  group_by(plot_title) %>%
  summarise(
    x = 2.5, 
    y = max(Zscore, na.rm = TRUE),
    label_text = paste0(unique(label_p), "\n", unique(label_b)),
    x_axis_label = paste0(unique(lead_snp), "\nGRCh38: ", unique(variant_pos)),
    .groups = "drop"
  ) %>%
  select(-x, -y) %>%
  mutate(
    x_pos = case_when(
      plot_title == "KANSL1 in CD4 cells" ~ 0.7,
      plot_title == "NAPSA in classical monocytes" ~ 2.5,
      plot_title == "NAPSA in non-classical monocytes" ~ 2.3,
      plot_title == "ZGLP1 in classical monocytes" ~ 2.3
    ),
    y_pos = case_when(
      plot_title == "KANSL1 in CD4 cells" ~ -2,
      plot_title == "NAPSA in classical monocytes" ~ 1.5,
      plot_title == "NAPSA in non-classical monocytes" ~ 1.5,
      plot_title == "ZGLP1 in classical monocytes" ~ -1.5
    )
  )

eqtl_box$x_label <- paste0(eqtl_box$x_label1, "\nGRCh38: ", eqtl_box$x_label2)

desired_order <- c("NAPSA in non-classical monocytes", "NAPSA in classical monocytes", "ZGLP1 in classical monocytes", "KANSL1 in CD4 cells")  # replace with your actual order

plots <- eqtl_box %>%
  mutate(plot_title = factor(plot_title, levels = desired_order)) %>%
  arrange(plot_title) %>%
  group_split(plot_title)


plot_list <- lapply(plots, function(df) {
  xlab_text <- unique(df$x_label)
  title_text <- unique(df$plot_title)
  
  ggplot(df, aes(x = genotype_label, y = Zscore)) +
    geom_boxplot(aes(fill = cell_type), outlier.shape = NA,  width = 0.4, size = 0.8) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
    # Top annotation with beta/p-values
    geom_text(data = label_df %>% filter(plot_title == title_text),
              aes(x = x_pos, y = y_pos, label = label_text),
              inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 4) +
    scale_fill_manual(values = celltype_colors) +
    labs(
      x = xlab_text,        # TRUE x-axis title
      y = "Normalised Gene Expression",
      title = title_text    # Plot title
    ) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5))
})

# Combine plots in a grid (2x2)
fig2c <- wrap_plots(plot_list, ncol = 2)



# Figure 2d ####
# UMAP showing NAPSA expression

umap_napsa <- fread('data/main/fig2d.csv')
umap_napsa <- umap_napsa[order(umap_napsa$NAPSA), ]

fig2d <- ggplot(umap_napsa, aes(x = UMAP1,
                                y = UMAP2,
                                colour = NAPSA)) + 
    geom_point() + 
    scale_color_gradient2(low = "#EBFAF9", mid = "#bbf0ebff", high = "#0e786dff", midpoint=0.5) +
    theme_classic(base_size = 14) + 
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(), 
          legend.position = "none"
    ) +
    labs(x = "UMAP1", y = "UMAP2", title = "UMAP Plot of NAPSA Expression", colour = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 6)))


# Figure 2 Combined

fig2 <- (fig2a | wrap_elements(final_fig2b)) /
        (fig2c | fig2d) +
  plot_annotation(tag_levels = "a")

ggsave("figures/main/Fig2.png", fig2, width = 18, height = 14, dpi = 300)




# ============================
#        MAIN FIGURE 3       #
# ============================


# Figure 3a ####
# Heatmap showing COMBAT eQTL and GWAS variant intersections

combat_gwas <- fread('data/main/fig3a.csv')
trait_order <- read_lines('data/main/fig3_trait_order.csv')

rename_traits <- c(
  "Chronic inflammatory diseases (ankylosing spondylitis, Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)" = "Chronic inflammatory diseases",
  "Rheumatoid arthritis (rheumatoid factor and/or anti-cyclic citrullinated peptide seropositive)" = "Rheumatoid arthritis"
)


# COMBAT GWAS all traits (filtered on top 10 traits with highest number of COMBAT-only hits)
fig3a <- combat_gwas %>%
  mutate(GWAS_TRAIT = factor(GWAS_TRAIT, levels = trait_order)) %>%
  ggplot(aes(EQTL_GENE_NAME, GWAS_TRAIT)) +
  geom_tile(fill = "#60BC40", color = "white") +
  geom_text(
    data = subset(combat_gwas, EQTL_CELL_count > 1),
    aes(label = EQTL_CELL_count),
    color = "black", size = 3
  ) +
  labs(
    title = "Heatmap of eQTL Cell Count - COMBAT All Traits",
    x = "eQTL Gene", y = "GWAS Trait"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_y_discrete(
    labels = ~ str_wrap(recode(., !!!rename_traits), width = 30)
  )




# Figure 3b ####
# COMBAT ONEK1K GWAS hits barplot

top_combat_hits <- fread('data/main/combat_onek1k_gwas_barplot.csv')

fig3b <- top_combat_hits %>%
  filter(type != 'Neither') %>%
  
  mutate(type = factor(type,
                       levels = c("OneK1K only",
                                  "COMBAT and OneK1K",
                                  "COMBAT only"))) %>%
  
  # Define ordering metric per trait
  group_by(trait) %>%
  mutate(order_value = max(hit_combat)) %>%
  ungroup() %>%
  mutate(trait = factor(trait,
                        levels = unique(trait[order(-order_value)]))) %>%
  
  ggplot(aes(x = raw_value, y = trait, fill = type)) +
  geom_col(position = "stack") +
  geom_text(aes(label = ifelse(raw_value == 0, "", raw_value)),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white") +
  labs(title = "GWAS Variants Found in eQTL Datasets by Trait",
       y = "Trait",
       x = "GWAS Variants in eQTL Datasets") +
  theme_classic() +
  theme(axis.text.y = element_text(hjust = 1),
      legend.position = c(0.85, 0.85),
      legend.text = element_text(size = 10),
      legend.key.size = unit(1.2, "lines"),
      legend.background = element_blank())    + 
  scale_y_discrete(
    labels = function(x) stringr::str_wrap(recode(x, !!!rename_traits), width = 30)
  ) +
  
  scale_fill_manual(values = c(
    "COMBAT only" = "#60BC40",
    "OneK1K only" = "#2683C6",
    "COMBAT and OneK1K" = "#DEC64E"
  ),
  name = NULL)


# Figure 3c ####
# eQTL plots for COMBAT hit genes


eqtl_box_fig3 <- fread('data/main/fig3c.csv')

label_df <- eqtl_box_fig3 %>%
  group_by(plot_title) %>%
  summarise(
    x = 2.5, 
    y = max(Zscore, na.rm = TRUE),
    label_text = paste0(unique(label_p), "\n", unique(label_b)),
    x_axis_label = paste0(unique(lead_snp), "\nGRCh38: ", unique(variant_pos)),
    .groups = "drop"
  ) %>%
  select(-x, -y) %>%
  mutate(
    x_pos = case_when(
      plot_title == "REL in CD4 cells" ~ 1,
      plot_title == "IRF5 in non-classical monocytes" ~ 2.4,
      plot_title == "TRAF1 in CD8 cells" ~ 2.4
    ),
    y_pos = case_when(
      plot_title == "REL in CD4 cells" ~ -2,
      plot_title == "IRF5 in non-classical monocytes" ~ -2,
      plot_title == "TRAF1 in CD8 cells" ~ -2.5
    )
  )

eqtl_box_fig3$x_label <- paste0(eqtl_box_fig3$x_label1, "\nGRCh38: ", eqtl_box_fig3$x_label2)

desired_order <- c("REL in CD4 cells", "IRF5 in non-classical monocytes", "TRAF1 in CD8 cells")

plots3 <- eqtl_box_fig3 %>%
  mutate(plot_title = factor(plot_title, levels = desired_order)) %>%
  arrange(plot_title) %>%
  group_split(plot_title)


plot_list3 <- lapply(plots3, function(df) {
  xlab_text <- unique(df$x_label)
  title_text <- unique(df$plot_title)
  
  ggplot(df, aes(x = genotype_label, y = Zscore)) +
    geom_boxplot(aes(fill = cell_type), outlier.shape = NA, width = 0.4, size = 0.8) +  # <- narrow boxes
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
    geom_text(data = label_df %>% filter(plot_title == title_text),
              aes(x = x_pos, y = y_pos, label = label_text),
              inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 4) +
    scale_fill_manual(values = celltype_colors) +
    labs(
      x = xlab_text,
      y = "Normalised Gene Expression",
      title = title_text
    ) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5))
})


fig3c <- wrap_plots(plot_list3, ncol = 3)


# Figure 3d coloc plots for REL, IRF5 and TRAF1

gene_data_plot <- fread('data/main/fig3d_plot.csv')
gene_data_plot$facet_title <- paste0(gene_data_plot$cell, ' (', gene_data_plot$dataset, ')')

genes <- c('REL', 'IRF5', 'TRAF1')

plot_list <- lapply(genes, function(gene_i) {

    df_plot  <- gene_data_plot  %>% filter(gene == gene_i)

    ymax <- max(-log10(df_plot$P_VALUE), na.rm = TRUE) * 1.05

    p <- ggplot(df_plot, aes(x = POS, y = -log10(P_VALUE), color = dataset)) +
        geom_point(size=1) +
        facet_wrap(~facet_title, ncol = 1) +
        scale_color_manual(values = c("COMBAT"="#60BC40", "ONEK1K"="#2683C6")) +
        labs(
            x = "Position",
            y = "-log10(p-value)",
            title = gene_i
        ) +
        guides(color = guide_legend(override.aes = list(size = 5))) +
        ylim(0, ymax) +
        labs(color = NULL) +
        theme_classic() +
        theme(
            legend.position = "none",
            strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5)
        )

    p
})

fig3d <- wrap_plots(plot_list, ncol = 3)



# Figure 3e UMAP of REL, IRF5 and TRAF1

umap_3gene <- fread('data/main/fig3e.csv') %>% select(-V1)

umap_long3 <- umap_3gene %>%
  pivot_longer(
    cols = c(REL, IRF5, TRAF1,),  # all gene columns
    names_to = "GENE",
    values_to = "EXPRESSION"
  )

umap_long3$GENE <- factor(umap_long3$GENE, levels = c("REL", "IRF5", "TRAF1"))
umap_long3 <- umap_long3[order(umap_long3$EXPRESSION), ]

genes <- unique(umap_long3$GENE)

# Define three colour scales
col_scales <- list(
  scale_colour_gradient(low = "#EBFAF9", high = "#8b7823ff"),
  scale_colour_gradient(low = "#fde0dd", high = "#c51b8a"),
  scale_colour_gradient(low = "#e5f5e0", high = "#210e78ff")
)

plot_list <- lapply(seq_along(genes), function(i) {
  
  ggplot(umap_long3 %>% dplyr::filter(GENE == genes[i]),
         aes(x = UMAP1, y = UMAP2, colour = EXPRESSION)) +
    geom_point() +
    col_scales[[i]] +          # <-- different colour scale per gene
    theme_classic(base_size = 18) +
    labs(title = genes[i]) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "plain", size = 18)
    )
})

fig3e <- plot_list[[3]] | plot_list[[1]] | plot_list[[2]]


# Figure 3 combined

fig3 <- (fig3a) /
        (fig3b) /
        (fig3c) /
        (fig3d) /
        (fig3e) +
  plot_annotation(tag_levels = "a")

ggsave("figures/main/Fig3.png", fig3, width = 14, height = 20, dpi = 300)



# ============================
#        MAIN FIGURE 4       #
# ============================


# -------------------------------
# Figure 4a - Interaction eQTL boxplots: cMono
# -------------------------------
rps26_sev <- fread('data/main/fig4c_sev.csv')
hosp_colors <- c("Non-hospitalised" = "#60BC40", "Hospitalised" = "#DC320A")

label_df <- tibble(
  cell_type = c("cMono", "NK", "CD8"),
  label = c("FDR = 0.00356\nβint = 0.198", 
            "FDR = 0.0375\nβint = 0.238",
            "FDR = 0.0375\nβint = 0.889"),
  x = c(1, 1, 1),
  y = c(-2, -2, -2)
)

plot_df <- rps26_sev %>% filter(cell_type == "cMono")
plot_df$hosp <- factor(plot_df$hosp)
label_plot <- label_df %>% filter(cell_type == "cMono")

p4a <- ggplot(plot_df, aes(x = genotype_label, y = Zscore, fill = hosp)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(colour = "black", size = 1.5, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8)) +
  geom_text(data = label_plot, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 0, size = 4.5) +
  scale_fill_manual(values = hosp_colors) +
  labs(x = "rs1131017\nGRCh38: chr12:56,042,145", y = "Normalised Gene Expression",
       fill = "Infection Severity", title = "RPS26 in cMono (Infection Severity)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, vjust = 1.5),
        legend.position = "bottom")

# -------------------------------
# Figure 4b - Interaction eQTL boxplots: NK
# -------------------------------
plot_df <- rps26_sev %>% filter(cell_type == "NK")
plot_df$hosp <- factor(plot_df$hosp)
label_plot <- label_df %>% filter(cell_type == "NK")

p4b <- ggplot(plot_df, aes(x = genotype_label, y = Zscore, fill = hosp)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(colour = "black", size = 1.5, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8)) +
  geom_text(data = label_plot, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 0, size = 4.5) +
  scale_fill_manual(values = hosp_colors) +
  labs(x = "rs1131017\nGRCh38: chr12:56,042,145", y = "Normalised Gene Expression",
       fill = "Infection Severity", title = "RPS26 in NK (Infection Severity)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, vjust = 1.5),
        legend.position = "bottom")

# -------------------------------
# Figure 4c - Beta plot
# -------------------------------
beta_comp <- fread('data/main/fig4a.csv')

beta_plot <- beta_comp %>%
  mutate(sig = case_when(
    p_value < 0.0005 ~ "***",
    p_value < 0.005  ~ "**",
    p_value < 0.05   ~ "*",
    TRUE             ~ ""
  ))

slopes <- tibble(
  cell_type = c("cMono", "NK"),
  slope = c(0.198, 0.238),
  intercept = c(-1.8, -2.2)
)

facet_labels <- slopes %>%
  mutate(
    label = paste0(cell_type, " interaction\nbeta: ", slope),
    x = max(as.numeric(factor(beta_plot$inf_sev))),
    y = intercept - 0.05
  )

p4c <- ggplot(beta_plot, aes(x = factor(inf_sev), y = beta, color = cell_type)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = beta - slope_se, ymax = beta + slope_se),
                width = 0.2, size = 1, position = position_dodge(width = 0.5)) +
  geom_text(aes(label = sig, y = beta + slope_se + 0.1),
            position = position_dodge(width = 0.5), size = 5) +
  geom_abline(data = slopes, aes(slope = slope, intercept = intercept),
              color = "red", linetype = "dashed", size = 0.8) +
  geom_text(data = facet_labels,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            hjust = 1, vjust = 1,
            size = 4, color = "red") +
  facet_wrap(~cell_type) +
  scale_color_manual(values = celltype_colors) +
  labs(x = "Infection Severity", y = "Slope (beta)", color = "Cell Type",
       title = "RPS26 eQTL Effect Size by Infection Severity") +
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none")

# -------------------------------
# Figure 4d - Interaction eQTL boxplots: CD8
# -------------------------------
rps26_source <- fread('data/main/fig4c_source.csv')
plot_df <- rps26_source %>% filter(cell_type == "CD8")
plot_df$inf_source <- factor(plot_df$inf_source, levels = c("0","1"))
label_plot <- label_df %>% filter(cell_type == "CD8")
source_colors <- c("0" = "#2683C6", "1" = "#C62683")
source_labels <- c("0" = "Sepsis", "1" = "COVID-19")

p4d <- ggplot(plot_df, aes(x = genotype_label, y = Zscore, fill = inf_source)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(colour = "black", size = 1.5, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8)) +
  geom_text(data = label_plot, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 0, size = 4.5) +
  scale_fill_manual(values = source_colors, labels = source_labels) +
  labs(x = "rs1131017\nGRCh38: chr12:56,042,145", y = "Normalised Gene Expression",
       fill = "Infection Source", title = "RPS26 in CD8 (Infection Source)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, vjust = 1.5),
        legend.position = "bottom")

# -------------------------------
# Figure 4e - Scatter plot: RPS26 / CD8
# -------------------------------
scatter <- fread('data/main/fig4b.csv')
scatter_mean <- fread('data/main/fig4b_mean.csv')
label_rps26 <- tibble(
  cell_type = "CD8",
  gene = "RPS26",
  x_title = "rs1131017\nGRCh38: chr12:56,042,145",
  label = "FDR = 0.0375\n βint = 0.889",
  x = 1.5,
  y = 2.6
)
scatter_rps26 <- scatter %>% filter(gene == "RPS26", cell_type == "CD8")
mean_rps26 <- scatter_mean %>% filter(gene == "RPS26", cell_type == "CD8")
geno_levels <- unique(scatter_rps26$genotype_label)
geno_colors <- c("#046C9A", "#C641ED", "#5BBCD6")[seq_along(geno_levels)]

p4e <- ggplot(scatter_rps26, aes(x = factor(inf_source), y = Zscore, colour = genotype_label)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 4) +
  geom_line(data = mean_rps26, 
            aes(x = factor(x), y = mean_Z, group = genotype_label, colour = genotype_label),
            inherit.aes = FALSE, linetype = "solid", linewidth = 1.5) +
  geom_text(data = label_rps26, aes(x = x-1, y = y, label = label),
            inherit.aes = FALSE, hjust = -0.2, vjust = 0.5, size = 5) +
  scale_colour_manual(values = geno_colors) +
  scale_x_discrete(labels = c("0" = "Sepsis", "1" = "COVID-19")) +
  labs(x = label_rps26$x_title, y = "Normalised Gene Expression",
       colour = "Genotype", title = "RPS26 in CD8") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

# -------------------------------
# Figure 4f - Scatter plot: ADAM10 / CD4
# -------------------------------
label_adam10 <- tibble(
  cell_type = "CD4",
  gene = "ADAM10",
  x_title = "rs7161799\nGRCh38: chr15:58,478,324",
  label = "FDR = 0.0457\n βint = 1.05",
  x = 1.5,
  y = 2.6
)
scatter_adam10 <- scatter %>% filter(gene == "ADAM10", cell_type == "CD4")
mean_adam10 <- scatter_mean %>% filter(gene == "ADAM10", cell_type == "CD4")
geno_levels <- unique(scatter_adam10$genotype_label)
geno_colors <- c("#046C9A", "#C641ED", "#5BBCD6")[seq_along(geno_levels)]

p4f <- ggplot(scatter_adam10, aes(x = factor(inf_source), y = Zscore, colour = genotype_label)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 4) +
  geom_line(data = mean_adam10, 
            aes(x = factor(x), y = mean_Z, group = genotype_label, colour = genotype_label),
            inherit.aes = FALSE, linetype = "solid", linewidth = 1.5) +
  geom_text(data = label_adam10, aes(x = x-1, y = y, label = label),
            inherit.aes = FALSE, hjust = -0.2, vjust = 0.5, size = 5) +
  scale_colour_manual(values = geno_colors) +
  scale_x_discrete(labels = c("0" = "Sepsis", "1" = "COVID-19")) +
  labs(x = label_adam10$x_title, y = "Normalised Gene Expression",
       colour = "Genotype", title = "ADAM10 in CD4") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

# -------------------------------
# Combine all plots in final Figure 4 layout
# -------------------------------
fig4 <- (p4a | p4b | p4c) /
        (p4d | p4e | p4f) +
  plot_annotation(tag_levels = "a")

# Save final figure
ggsave("figures/main/Fig4.png", fig4, width = 18, height = 18, dpi = 300)



# Supplementary


# ============================
#   SUPPLEMENTARY FIGURE 1   #
# ============================
# eGene jaccard index heatmap

egene_overlap <- fread('data/supplementary/eqtl_overlap_heatmap_data.csv')

cell_vector <- c('CD4' ,'CD8', 'NK', 'cMono', 'ncMono', 'B')
egene_overlap$Cell1 <- factor(egene_overlap$Cell1, levels=cell_vector)
egene_overlap$Cell2 <- factor(egene_overlap$Cell2, levels=cell_vector)


overlap <- egene_overlap %>% 
  filter(Jaccard != 1) %>% 
  ggplot(aes(Cell1, Cell2, fill = Jaccard)) +
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, hjust=1)) +
  geom_text(aes(label = shared_egenes), color = "white", size = 4) +
  theme(
    strip.text.x = element_text(
      size = 10, color = "black", face = "italic"
    )
  ) +
  theme(plot.title = element_text(face = "italic", size=15)) +
  theme_classic(base_size = 16) +
  labs(fill = "Jaccard \nIndex (%)") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 12),
        text = element_text(family = "sans")) +
  scale_fill_gradientn(colors = c("#78B7C5", "#273046")) +
  labs(x = "", y = "", title = "")

ggsave(
  filename = 'figures/supplementary/Supplementary_Fig1.png',
  plot = overlap,
  width = 6,
  height = 5,
  dpi = 300
)

# ============================
#   SUPPLEMENTARY FIGURE 2   #
# ============================
# eGene upset

# Read the data
upset_input <- fread('data/supplementary/egene_celltype_upset_data.csv')

# Convert TRUE/FALSE to 1/0 for UpSetR
upset_matrix <- as.data.frame(lapply(upset_input[, -1, with = FALSE], function(x) as.integer(x)))

# Ensure column order
sets <- colnames(upset_matrix)

# Save the plot
png("figures/supplementary/Supplementary_Fig2.png",
    width = 2400, height = 1600, res = 200)

# Make the UpSet plot
UpSetR::upset(
  upset_matrix,
  sets = c("CD4", "cMono", "CD8", "NK", "B", "ncMono"),           # explicitly provide all sets
  keep.order = TRUE,     # preserve the column order
  nintersects = 90,
  text.scale = 1.5,
  set_size.show = TRUE,
  sets.bar.color = "#046C9A",
  main.bar.color = "#046C9A",
  matrix.color = "#046C9A"
)

grid::grid.text(
  "eGene Sharing Across Cell Types",
  x = unit(0.7, "npc"),
  y = unit(0.98, "npc"),
  just = c("right", "top"),
  gp = grid::gpar(fontsize = 20)
)

dev.off()




# ============================
#   SUPPLEMENTARY FIGURE 3   #
# ============================
# UMAP showing ZGLP1, KANSL1, ADAM10, RPS26 expression

umap_4gene <- fread('data/supplementary/4gene_umap.csv') %>% select(-V1)

umap_long <- umap_4gene %>%
  pivot_longer(
    cols = c(ZGLP1, KANSL1, ADAM10, RPS26),
    names_to = "GENE",
    values_to = "EXPRESSION"
  )

# Set facet order
umap_long$GENE <- factor(
  umap_long$GENE,
  levels = c("ZGLP1", "KANSL1", "ADAM10", "RPS26")
)

umap_long <- umap_long[order(umap_long$EXPRESSION), ]

facet_labels <- umap_long %>%
  group_by(GENE) %>%
  summarise(
    x = min(UMAP1),
    y = max(UMAP2)
  ) %>%
  mutate(
    label = letters[seq_len(n())]
  )

supp_fig_umap_4gene <- ggplot(
  umap_long,
  aes(
    x = UMAP1,
    y = UMAP2,
    colour = EXPRESSION
  )
) + 
  geom_point() + 
  scale_color_gradient2(
    low = "#EBFAF9",
    mid = "#bbf0ebff",
    high = "#0e786dff",
    midpoint = 0.5
  ) +
  theme_classic(base_size = 18) + 
  facet_wrap(~GENE, scales = "free") +
  geom_text(
    data = facet_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.4,
    vjust = 1.2,
    size = 6
  ) +
  theme(
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(), 
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "plain", size = 18)
  )

ggsave(
  "figures/supplementary/Supplementary_Fig3.png",
  supp_fig_umap_4gene,
  width = 16, height = 12, dpi = 300)


# ============================
#   SUPPLEMENTARY FIGURE 4   #
# ============================
# REPLICATION ANALYSIS

# Bar plot showing cohort composition across three studies

# Data on cohort
sample_counts <- tribble(
  ~dataset, ~phenotype, ~count,
  "COMBAT", "healthy", 8,
  "COMBAT", "asymptomatic COVID-19", 11,
  "COMBAT", "mild COVID-19", 17,
  "COMBAT", "critical COVID-19", 13,
  "COMBAT", "severe COVID-19", 32,
  "COMBAT", "sepsis", 19,
  
  "WANG", "most severe COVID-19", 501,
  "WANG", "severe COVID-19", 494,
  "WANG", "mild COVID-19", 332,
  "WANG", "asymptomatic COVID-19", 78,
  
  "ONEK1K", "healthy", 982
)

# --- Define phenotype order (mild → severe) ---
sample_counts <- sample_counts %>%
  mutate(
    phenotype = factor(
      phenotype,
      levels = c(
        "healthy",
        "asymptomatic COVID-19",
        "mild COVID-19",
        "severe COVID-19",
        "most severe COVID-19",
        "critical COVID-19",
        "sepsis"
      ),
      ordered = TRUE
    )
  )

# --- Define color palette ---
phenotype_colors <- c(
  "healthy" = "#a6cee3",
  "asymptomatic COVID-19" = "#b2df8a",
  "mild COVID-19" = "#33a02c",
  "severe COVID-19" = "#e31a1c",
  "most severe COVID-19" = "#6a3d9a",
  "critical COVID-19" = "#fb9a99",
  "sepsis" = "#ff7f00"
)

# --- Plot ---
p_cohort <- ggplot(sample_counts, aes(x = dataset, y = count, fill = phenotype)) +
  geom_bar(stat = "identity", position = "stack",  linewidth = 0.2) +
  scale_fill_manual(values = phenotype_colors, drop = FALSE) +
  labs(
    title = "Sample composition across studies",
    x = "Study",
    y = "Number of samples",
    fill = "Clinical phenotype"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = c(0.2, 0.6),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# --- Save outputs ---
write_csv(sample_counts, "data/supplementary/cohort_comparison_barchart_data.csv")


# Upset plot showing replicated eGenes between datasets

upset_data <- fread('data/supplementary/egene_cohort_replication_upset_data.csv')

# Define custom colors for each dataset and intersection

set_order <- c("ONEK1K", "WANG", "COMBAT")

color_set <- c(
  COMBAT = "#60BC40",
  ONEK1K = "#2683C6",
  WANG = "#C62683"
)[set_order]

color_intersects <- c(
  "ONEK1K" = "#2683C6",
  "WANG" = "#C62683",
  "ONEK1K&COMBAT" = "#DEC64E",
  "ONEK1K&WANG" = "#8C4EDC",   
  "COMBAT" = "#60BC40",
  "ONEK1K&WANG&COMBAT" = "#555555",
  "WANG&COMBAT" = "#C66926"

)

p_upset_genes <- UpSetR::upset(
  upset_data,
  sets = set_order,
  nsets = 3,
  nintersects = NA,
  order.by = "freq",
  sets.bar.color = color_set,     
  main.bar.color = unname(color_intersects), 
  text.scale = 1.3,
  mainbar.y.label = "Shared eGenes",
  sets.x.label = "Total eGenes per dataset"
)

# Save
png("figures/supplementary/Supplementary_Fig4b.png",
    width  = 8 * 300,
    height = 5 * 300,
    res    = 300)
print(p_upset_genes)
dev.off()


# Comparison scatter plot

scatter_compare <- fread('data/supplementary/comparison_scatter_top_6_data.csv')
scatter_compare_wang_onek <- fread('data/supplementary/comparison_scatter_wang_onekek.csv')


# COMBAT-ONEK1K slope comparison plot
scatter_compare_onek1k <- ggplot(scatter_compare %>% filter(other_dataset=='ONEK1K'), 
                  aes(x = COMBAT_slope, y = other_slope,
                      color = COMBAT_cell)) +
  geom_point(alpha=0.7) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", color = "grey40") +
  scale_color_manual(values = celltype_colors) +
  guides(color = "none",
    shape = guide_legend(title = "Dataset")) +
  facet_wrap(~ COMBAT_cell, scales = "free", ncol = 3) +
  labs(
    x = "COMBAT slope",
    y = "ONEK1K slope",
    color = "COMBAT cell type",
    shape = "ONEK1K cell type",
    title = "Slope comparison between COMBAT and ONEK1K"
  ) +
  theme_classic() +
  theme(
    #axis.text.x = element_blank(),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 16),
    strip.background = element_blank(),            # remove facet background
    strip.text = element_text(size = 14),
    #axis.ticks.x = element_blank()
  )

# COMBAT-WANG slope comparison plot
scatter_compare_wang <- ggplot(scatter_compare %>% filter(other_dataset=='ONEK1K'),
                 aes(x = COMBAT_slope, y = other_slope,
                     color = COMBAT_cell)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", color = "grey40") +
  scale_color_manual(values = celltype_colors) +
  guides(color = "none",
    shape = guide_legend(title = "Dataset")) +
  facet_wrap(~ COMBAT_cell, scales = "free", ncol = 3) +
  labs(
    x = "COMBAT slope",
    y = "WANG slope",
    color = "COMBAT cell type",
    #shape = "WANG cell type",
    title = "Slope comparison between COMBAT and WANG"
  ) +
  theme_classic() +
  theme(
    #axis.text.x = element_blank(),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 16),
    strip.background = element_blank(),            # remove facet background
    strip.text = element_text(size = 14),
    #axis.ticks.x = element_blank()
  )

# WANG-ONEK1K slope comparison plot
# Note: It is for the above gene/snp combinations
scatter_compare_wang_onek1k <- ggplot(scatter_compare_wang_onek,
                 aes(x = onek1k_slope, y = wang_slope,
                     color = COMBAT_cell)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", color = "grey40") +
  scale_color_manual(values = celltype_colors) +
  guides(color = "none",
    shape = guide_legend(title = "Dataset")) +
  facet_wrap(~ COMBAT_cell, scales = "free", ncol = 3) +
  labs(
    x = "ONEK1K slope",
    y = "WANG slope",
    color = "COMBAT cell type",
    #shape = "WANG cell type",
    title = "Slope comparison between ONEK1K and WANG"
  ) +
  theme_classic() +
  theme(
    #axis.text.x = element_blank(),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 16),
    strip.background = element_blank(),            # remove facet background
    strip.text = element_text(size = 14),
    #axis.ticks.x = element_blank()
  )


combined_plot <- scatter_compare_onek1k / scatter_compare_wang /  scatter_compare_wang_onek1k + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")


# Replicated egenes between cohorts (stacked barplot)

egene_rep <- fread('data/supplementary/egene_replication_barplot.csv')

egene_rep <- egene_rep %>%
  mutate(COMBAT_cell = factor(COMBAT_cell, 
                              levels = c("CD4", "cMono", "CD8", "NK", "B", "ncMono")))

rep_plot <- ggplot(egene_rep, aes(x = COMBAT_cell, y = prop, fill = status)) +
  geom_col() +
  geom_text(aes(label = N), 
            position = position_stack(vjust = 0.5), size = 3, color = "black") +
  facet_wrap(~ dataset, scales = "free_x") +
  scale_fill_manual(values = c(
    "COMBAT only"        = "#60BC40",
    "COMBAT and OneK1K"  = "#DEC64E",
    "COMBAT and Wang"    = "#C66926"
  )) +
  labs(
    x = "COMBAT cell type",
    y = "Proportion of eGenes",
    fill = "Replication Status  ",
    title = "Replication of COMBAT eGenes in OneK1K and WANG"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )


# Combine into supplementary figure

upset_png <- png::readPNG("figures/supplementary/Supplementary_Fig4b.png")
upset_grob <- grid::rasterGrob(upset_png, width = unit(1, "npc"), height = unit(1, "npc"))
upset_plot <- ggdraw() + draw_grob(upset_grob)

left_col <- p_cohort / upset_plot / rep_plot
right_col <- scatter_compare_onek1k / scatter_compare_wang / scatter_compare_wang_onek1k

final_layout <- (
  (p_cohort / upset_plot / rep_plot) |             # left column
  (scatter_compare_onek1k / scatter_compare_wang / scatter_compare_wang_onek1k) # right column
) + 
  plot_layout(
    widths = c(1, 1),
    heights = c(1, 1.5, 1) # middle-left taller
  ) +
  plot_annotation(tag_levels = 'a')

# Save figure
ggsave(
  filename = "figures/supplementary/Supplementary_Fig4.png",
  plot = final_layout,
  width = 16,
  height = 15,
  dpi = 300
)



# ============================
#   SUPPLEMENTARY FIGURE 5   #
# ============================

# --- Load processed data ---
plot_data <- readRDS("data/supplementary/coloc_plot_data.rds")
coloc_ok_filtered <- plot_data$coloc_ok_filtered
highlight_rows <- as.data.table(plot_data$highlight_rows)
coloc_wang_filtered <- plot_data$coloc_wang_filtered

cell_order <- c("CD4", "CD8", "NK", "cMono", "ncMono", "B")
coloc_ok_filtered[, combat_cell := factor(combat_cell, levels = cell_order)]
coloc_wang_filtered[, combat_cell := factor(combat_cell, levels = cell_order)]
highlight_rows[, combat_cell := factor(combat_cell, levels = cell_order)]

# Define mapping from combat_cell to oneK1K subtypes
combat_to_onek1k <- list(
  CD4   = c("CD4 Effector memory/TEMRA", "CD4 Naive/Central memory T cell", "CD4 SOX4 T cell"),
  CD8   = c("CD8 Effector memory", "CD8 Naive/Central memory T cell", "CD8 S100B T cell"),
  cMono = c("Classic Monocyte"),
  ncMono = c("Non-classic Monocyte"),
  NK = c("Natural Killer Cell", "Natural Killer Recruiting Cell"),
  B = c("Naïve/Immature B Cell", "Memory B Cell", "Plasma Cell")
)



# Generate oneK1K colors
onek1k_colors <- list()

for (combat in names(combat_to_onek1k)) {
  subtypes <- combat_to_onek1k[[combat]]
  n <- length(subtypes)
  main_col <- celltype_colors[combat]
  
  if (n == 1) {
    # Single subtype: just use the main color
    onek1k_colors[[combat]] <- setNames(main_col, subtypes)
  } else {
    # Generate visually distinct shades around the main color
    # Convert main color to HCL to get its hue
    hcl_main <- as(hex2RGB(main_col), "polarLUV")@coords
    hue_main <- hcl_main[1, "H"]  # hue in degrees
    
    # Generate colors evenly spaced in lightness, keeping same hue and chroma
    lightness_vals <- seq(40, 80, length.out = n)  # adjust min/max lightness
    chroma_val <- 50  # adjust chroma to taste
    cols <- hcl(h = hue_main, c = chroma_val, l = lightness_vals)
    
    onek1k_colors[[combat]] <- setNames(cols, subtypes)
  }
}

# Combine into a single vector for ggplot
onek1k_colors_vec <- unlist(onek1k_colors)
names(onek1k_colors_vec) <- sub("^[^.]+\\.", "", names(onek1k_colors_vec))


coloc_ok_filtered[, onek1k_cell := as.character(onek1k_cell)]

# --- Create individual plots per facet for coloc_ok_filtered ---
p_list_ok <- lapply(unique(coloc_ok_filtered$combat_cell), function(facet){
  df <- subset(coloc_ok_filtered, combat_cell == facet)
  highlight_df <- subset(highlight_rows, combat_cell == facet)
  
  ggplot(df, aes(x = PP.H4, color = onek1k_cell, fill = onek1k_cell)) +
    geom_density(alpha = 0.4, size = 0.8) +
    scale_color_manual(values = onek1k_colors_vec) +
    scale_fill_manual(values = onek1k_colors_vec) +
    labs(x = "PP.H4",
         y = "Density",
         title = paste("", facet),
         color = "OneK1K Cell Type",  # <-- new legend title
         fill = "OneK1K Cell Type" ) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = c(0.5, 0.8),
          plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(fill = alpha("white", 0.9), color = NA)) +
    geom_vline(
      data = highlight_df,
      aes(xintercept = PP.H4),
      color = "black",
      linetype = "dashed",
      size = 0.8
    ) +
    geom_text(
      data = highlight_df,
      aes(x = PP.H4, y = y_pos, label = label_text, hjust = hjust_val),
      inherit.aes = FALSE,
      angle = 0,
      size = 4,
      color = "black"
    )
})

# Combine all per-facet plots into one patchwork for coloc_ok_filtered
p_coloc_ok_patch <- wrap_plots(p_list_ok, ncol = 3)

# --- Create individual plots per facet for Wang dataset ---
p_list_wang <- lapply(unique(coloc_wang_filtered$combat_cell), function(facet){
  df <- subset(coloc_wang_filtered, combat_cell == facet)
  
  ggplot(df, aes(x = PP.H4, color = combat_cell, fill = combat_cell)) +
    geom_density(alpha = 0.6, size = 0.8) +
    scale_color_manual(values = celltype_colors) +
    scale_fill_manual(values = celltype_colors) +
    labs(x = "PP.H4", y = "Density", title = paste("", facet)) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5))
})

# Combine all per-facet plots into one patchwork for Wang dataset
p_coloc_wang_patch <- wrap_plots(p_list_wang, ncol = 3)


# --- Combine patchwork plots ---
final_plot <- p_coloc_ok_patch / p_coloc_wang_patch + 
  plot_layout(heights = c(1,1)) +
  plot_annotation(tag_levels = "a")

ggsave("figures/supplementary/Supplementary_Fig5.png",
       final_plot,
       width = 20, height = 20, dpi = 300)

# ============================
#   SUPPLEMENTARY FIGURE 6   #
# ============================

# Combat GWAS heatmap COVID-19 traits
combat_gwas_covid <- fread('data/supplementary/combat_gwas_heatmap_covid.csv')

# COMBAT GWAS covid traits
combat_covid_heatmap <- combat_gwas_covid %>%
  ggplot(aes(x = EQTL_GENE_NAME, y = GWAS_TRAIT, fill = EQTL_CELL_count)) +
  geom_tile(color = "white") + 
  scale_fill_gradient(
  low = "#8AE06C", high = "#164008", na.value = "white",
  breaks = c(1, 2),
  labels = c("1", "2")
) +
  labs(title = "Heatmap of eQTL Cell Count - COMBAT COVID Traits",
       x = "eQTL Gene",
       y = "GWAS Trait",
       fill = "eQTL Cell Count") +
  theme_classic() +
  #scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(
    labels = function(x) stringr::str_wrap(recode(x, !!!rename_traits), width = 30)
  ) 

ggsave(
  filename = 'figures/supplementary/Supplementary_Fig6.png',
  plot = combat_covid_heatmap,
  width = 12,
  height = 6,
  dpi = 300
)


# ============================
#   SUPPLEMENTARY FIGURE 7   #
# ============================


# Onek1k GWAS all traits
onek1k_gwas_all <- fread('data/supplementary/onek1k_gwas_heatmap_all.csv')

# ONEK1K GWAS all traits (filtered on above studies)
onek1k_gwas_heatmap <- onek1k_gwas_all %>%
  arrange(EQTL_GENE_NAME) %>%
  #mutate(GWAS_TRAIT = factor(GWAS_TRAIT, levels = trait_order)) %>%
  mutate(gene_group = if_else(
    dense_rank(EQTL_GENE_NAME) <= max(dense_rank(EQTL_GENE_NAME)) / 2,
    "Group 1",
    "Group 2"
  )) %>%
  ggplot(aes(x = EQTL_GENE_NAME, y = GWAS_TRAIT, fill = EQTL_CELL_count)) +
  geom_tile(color = "white") + 
  scale_fill_gradient(low = "#72BCF2", high = "#083454", na.value = "white", labels = scales::number_format(accuracy = 1)) +# Use white borders for tiles
  labs(title = "Heatmap of eQTL Cell Count - OneK1K All Traits",
       x = "eQTL Gene",
       y = "GWAS Trait",
       fill = "eQTL Cell Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_blank()) +
  #scale_y_discrete(labels = function(x) str_wrap(x, width = 30))# Rotate x-axis labels for better visibility
  scale_y_discrete(
    labels = function(x) stringr::str_wrap(recode(x, !!!rename_traits), width = 30)
  ) +
  facet_wrap(~gene_group, scales = "free_x", ncol=1)

ggsave(
  filename = 'figures/supplementary/Supplementary_Fig7.png',
  plot = onek1k_gwas_heatmap,
  width = 16,
  height = 8,
  dpi = 300
)





# ============================
#   SUPPLEMENTARY FIGURE 8   #
# ============================

# Onek1k GWAS covid traits

onek1k_gwas_covid <- fread('data/supplementary/onek1k_gwas_heatmap_covid.csv')

onek1k_covid_heatmap <- onek1k_gwas_covid %>%
  arrange(EQTL_GENE_NAME) %>%
  #mutate(GWAS_TRAIT = factor(GWAS_TRAIT, levels = trait_order)) %>%
  mutate(gene_group = if_else(
    dense_rank(EQTL_GENE_NAME) <= max(dense_rank(EQTL_GENE_NAME)) / 2,
    "Group 1",
    "Group 2"
  )) %>%
  ggplot(aes(x = EQTL_GENE_NAME, y = GWAS_TRAIT, fill = EQTL_CELL_count)) +
  geom_tile(color = "white") + 
  scale_fill_gradient(low = "#72BCF2", high = "#083454", na.value = "white", labels = scales::number_format(accuracy = 1)) +# Use white borders for tiles
  labs(title = "Heatmap of eQTL Cell Count - ONEK1K COVID Traits",
       x = "eQTL Gene",
       y = "GWAS Trait",
       fill = "eQTL Cell Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_blank()) +
  #scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
  scale_y_discrete(
    labels = function(x) stringr::str_wrap(recode(x, !!!rename_traits), width = 30)
  ) +
  facet_wrap(~gene_group, scales = "free_x", ncol=1)


ggsave(
  filename = 'figures/supplementary/Supplementary_Fig8.png',
  plot = onek1k_covid_heatmap,
  width = 19,
  height = 15,
  dpi = 300
)

# ============================
#   SUPPLEMENTARY FIGURE 9   #
# ============================

# COMBAT WANG GWAS BAR PLOT

combat_wang_gwas <- fread('data/supplementary/combat_wang_gwas_barplot.csv')

combat_wang_barplot <- combat_wang_gwas %>%
  filter(type != 'Neither') %>%
  
  mutate(type = factor(type,
                       levels = c("Wang only",
                                  "COMBAT and Wang",
                                  "COMBAT only"))) %>%
 
  # Define ordering metric per trait
  group_by(trait) %>%
  mutate(order_value = max(hit_combat)) %>%   # or sum(hit_combat)
  ungroup() %>%
  
  # Reorder trait factor levels BEFORE plotting
  mutate(trait = factor(trait,
                        levels = unique(trait[order(-order_value)]))) %>%
  filter(!(hit_combat == 0 & hit_wang == 0 & hit_both == 0)) %>%
  mutate(trait_group = if_else(
    dense_rank(trait) <= max(dense_rank(trait)) / 2,
    "Group 1",
    "Group 2"
  )) %>%
  
  ggplot(aes(x = raw_value, y = trait, fill = type)) +
  geom_col(position = "stack") +
  geom_text(aes(label = ifelse(raw_value == 0, "", raw_value)),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white") +
  labs(title = "GWAS Variants Found in eQTL Datasets by Trait",
       y = "Trait",
       x = "GWAS Variants in eQTL Datasets") +
  theme_classic() +
  theme(axis.text.y = element_text(hjust = 1),
      legend.position = c(0.85, 0.85),  # inside top-right
      legend.text = element_text(size = 14),   # increase text size
      legend.key.size = unit(1.6, "lines"),    # increase colour box size
      legend.background = element_blank())   +
  theme(strip.text = element_blank()) +
  # Apply readable pretty trait names here (wrapped)
  scale_y_discrete(
    labels = function(x) stringr::str_wrap(recode(x, !!!rename_traits), width = 30)
  ) +
  
  scale_fill_manual(values = c(
    "COMBAT only" = "#60BC40",
    "Wang only" = "#C62683",
    "COMBAT and Wang" = "#C66926"
  ),
  name = NULL) +
  facet_wrap(~trait_group, scales = "free_y", nrow=1)


ggsave(
  filename = 'figures/supplementary/Supplementary_Fig9.png',
  plot = combat_wang_barplot,
  width = 20,
  height = 20,
  dpi = 300
)

# ============================
#   SUPPLEMENTARY FIGURE 10  #
# ============================

# Wang 

trait_order_wang <- read_lines('data/supplementary/wang_gwas_heatmap_trait_order.csv')

rename_traits_wang <- c(
  "Chronic inflammatory diseases (ankylosing spondylitis, Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)" = "Chronic inflammatory diseases",
  "Rheumatoid arthritis (rheumatoid factor and/or anti-cyclic citrullinated peptide seropositive)" = "Rheumatoid arthritis",
  "Systemic seropositive rheumatic diseases (Systemic sclerosis or systemic lupus erythematosus or rheumatoid arthritis or idiopathic inflammatory myopathies)" = "Systemic seropositive rheumatic diseases",
  "Juvenile idiopathic arthritis (oligoarticular or rheumatoid factor-negative polyarticular)" = "Juvenile idiopathic arthritis",
  "Rheumatoid arthritis (rheumatoid factor and anti-cyclic citrullinated peptide seronegative)" = "Rheumatoid arthritis"
)

# WANG GWAS all traits (filtered on above studies)

wang_gwas_filter <- fread('data/supplementary/wang_gwas_heatmap_all.csv')

wang_heatmap_all <- wang_gwas_filter %>%
  arrange(EQTL_GENE_NAME) %>%
  mutate(trait_group = if_else(
    dense_rank(GWAS_TRAIT) <= max(dense_rank(GWAS_TRAIT)) / 2,
    "Group 1",
    "Group 2"
  )) %>%
  ggplot(aes(x = EQTL_GENE_NAME, y = GWAS_TRAIT)) +
  geom_tile(color = "white", fill = "#C62683") + 
  labs(title = "Heatmap of eQTL Cell Count - Wang All Traits",
       x = "eQTL Gene",
       y = "GWAS Trait") +
       #fill = "eQTL Cell Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_blank()) +
  #scale_y_discrete(labels = function(x) str_wrap(x, width = 30))# Rotate x-axis labels for better visibility
  scale_y_discrete(
    labels = function(x) stringr::str_wrap(recode(x, !!!rename_traits_wang), width = 30)
  ) +
 facet_wrap(~trait_group, scales = "free_y", nrow=1)

ggsave(
  filename = 'figures/supplementary/Supplementary_Fig10.png',
  plot = wang_heatmap_all,
  width = 25,
  height = 16,
  dpi = 300
)


# ============================
#   SUPPLEMENTARY FIGURE 11  #
# ============================

# WANG GWAS COVID HEATMAP

wang_gwas_covid <- fread('data/supplementary/wang_gwas_heatmap_covid.csv')

wang_heatmap_covid <- wang_gwas_covid %>%
  arrange(EQTL_GENE_NAME) %>%
  ggplot(aes(x = EQTL_GENE_NAME, y = GWAS_TRAIT)) +
  geom_tile(color = "white", fill = "#C62683") + 
  labs(title = "Heatmap of eQTL Cell Count - Wang COVID Traits",
       x = "eQTL Gene",
       y = "GWAS Trait",
       fill = "eQTL Cell Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_blank()) +
  #scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
  scale_y_discrete(
    labels = function(x) stringr::str_wrap(recode(x, !!!rename_traits), width = 30)
  ) 


ggsave(
  filename = 'figures/supplementary/Supplementary_Fig11.png',
  plot = wang_heatmap_covid,
  width = 7,
  height = 4,
  dpi = 300
)


# ============================
#   SUPPLEMENTARY FIGURE 12  #
# ============================
# FGSEA

fgsea_data <- fread('data/supplementary/fgsea_heatmap_data.csv')

fgsea_heatmap <- ggplot(fgsea_data,
                        aes(x = cell, y = pathway, color = NES, size = -log10(pval))) +  # size mapped to p-value
    geom_point() +
    #facet_grid(context ~ NES_sign, scales = "free") +  # NES_sign rows, context columns
    facet_wrap(~description) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # NES still color
    theme_classic() +
    labs(
      x = "COMBAT cell type",
      y = "Pathway",
      fill = "NES",  # still keeping label as NES for legend
      title = "FGSEA Significant Pathways (Signed Analysis)"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "bottom",
      strip.background = element_blank()
    )

ggsave(
  filename = 'figures/supplementary/Supplementary_Fig12.png',
  plot = fgsea_heatmap,
  width = 10,
  height = 5,
  dpi = 300
)


# ============================
#   SUPPLEMENTARY FIGURE 13  #
# ============================

# Made by others. Figures stored in figures/supplementary

# ============================
#   SUPPLEMENTARY FIGURE 14  #
# ============================

# Made by others. Figures stored in figures/supplementary

# ============================
#   SUPPLEMENTARY FIGURE 15  #
# ============================
# CELL TYPE MARKERS UMAP

marker_umap <- fread('data/supplementary/umap_marker_genes.csv')

desired_order <- c(
  "CD3D (T cell marker)",
  "CD4 (CD4+ T cell marker)",
  "CD8A (CD8+ T cell marker)",
  "LST1 (Monocyte marker)",
  "CD14 (Classical Monocyte marker)",
  "FCGR3A (Non-classical Monocyte marker)",
  "NCAM1 (NK cell marker)",
  "NKG7 (NK cell marker)",
  "CD6 (NK cell marker)",
  "CD19 (B cell marker)",
  "MS4A1 (B cell marker)",
  "CD79A (B cell marker)"
)

marker_umap$title <- factor(marker_umap$title, levels = desired_order)

marker_umap <- marker_umap[order(marker_umap$expression), ]

vmin <- quantile(marker_umap$expression, 0.01)
vmax <- quantile(marker_umap$expression, 0.99)

supp_fig_umap_12gene <- ggplot(marker_umap, aes(x = UMAP1,
                                             y = UMAP2,
                                             colour = expression)) + 
    geom_point() + 
    scale_color_viridis(
      option = "plasma",
      limits = c(vmin, vmax),
      oob = squish
    ) +
    theme_classic(base_size = 18) + 
    facet_wrap(~title, scales = "free", ncol=3) +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(), 
          legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(face = "plain", size = 18))


ggsave("figures/supplementary/Supplementary_Fig15.png",
       supp_fig_umap_12gene,
       width = 18, height = 20, dpi = 300)



# ============================
#   SUPPLEMENTARY FIGURE 16  #
# ============================
# Marker gene dot plot

marker_dot <- fread('data/supplementary/dotplot_marker_genes.csv')

marker_panel <- list(
    "CD4" = c("CD3D", "CD4", "IL7R"),
    "CD8" = c("CD3D", "CD8A", "CD8B"),
    "cMono" = c("LST1", "CD14", "S100A8"),
    "ncMono" = c("LST1", "FCGR3A", "MS4A7"),
    "NK" = c("NKG7", "NCAM1", "CD6"),
    "B" = c("CD19", "MS4A1", "CD79A")
)

# Flatten to create the desired gene_label order
desired_order <- unlist(lapply(names(marker_panel), function(ct) {
    paste0(marker_panel[[ct]], " (", ct, ")")
}))

marker_dot$gene_label <- factor(marker_dot$gene_label, levels = desired_order)
marker_dot$cell_type <- factor(marker_dot$cell_type, levels = names(marker_panel))

marker_dotplot <- ggplot(marker_dot, aes(x=cell_type, y=gene_label)) +
    geom_point(aes(size=pct_expr, color=avg_expr)) +
    scale_color_viridis(option="plasma") +
    scale_size(range = c(1, 8)) +
    theme_classic(base_size = 16) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 16)
    ) +
    labs(
        x="Cell type",
        y="Marker gene",
        color="Average\nexpression",
        size="Fraction\nexpressing"
        )

ggsave("figures/supplementary/Supplementary_Fig16.png",
       marker_dotplot,
       width = 7, height = 8, dpi = 300)


# ============================
#   SUPPLEMENTARY FIGURE 17  #
# ============================

# Made in Powerpoint. Figure stored in figures/supplementary


# ============================
#   SUPPLEMENTARY FIGURE 18  #
# ============================

# CELL TYPE COUNTS
assay_count <- fread('data/supplementary/cell_counts.csv')

assay_count[, fill_cell_type :=
  ifelse(cell_type %in% names(celltype_colors),
         cell_type,
         "Other")
]

assay_count[, cell_type := reorder(cell_type, -cell_count)]


fig18a <- ggplot(assay_count, aes(x = cell_type,
                                  y = cell_count,
                                  fill = fill_cell_type)) +
  geom_bar(stat = "identity", width = 0.7, colour = "white") +
  geom_text(aes(label = cell_count), vjust = -0.5, size = 4) +
  scale_fill_manual(values = celltype_colors) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(size = 14)
  ) +
  labs(
    x = "",
    y = "Number of Cells Assayed",
    title = "Number of Cells Assayed per Cell Type"
  )


# CELL TYPE PROPORTIONS

cell_prop <- fread('data/supplementary/sample_cell_prop.csv') %>%
  mutate(COMBAT_ID = as.character(COMBAT_ID)) %>%
  group_by(Diagnosis) %>%
  mutate(
    COMBAT_ID = factor(COMBAT_ID,
                       levels = filter(pick(everything()), cell_type_plot == "CD4") %>%
                                arrange(desc(prop)) %>%
                                pull(COMBAT_ID))
  ) %>%
  ungroup() %>%
  mutate(cell_type_plot = factor(cell_type_plot,
                                 levels = c('Other', 'CD4', 'cMono', 'CD8', 'NK', 'B', 'ncMono')))

cell_prop_bar <- ggplot(cell_prop, aes(x = COMBAT_ID, y = prop, fill = cell_type_plot)) +
  geom_bar(stat = "identity", width = 0.8) +  # set bar width
  scale_fill_manual(values = celltype_colors,
                    breaks = c('Other', 'CD4', 'cMono', 'CD8', 'NK', 'B', 'ncMono')) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.background = element_blank(),            # remove facet background
    strip.text = element_text(size = 14),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 14)
  ) +
  labs(
    x = "Sample",
    y = "Proportion of cells",
    title = "Proportion of Cells per Sample",
    fill = ""
  ) +
  facet_grid(. ~ Diagnosis, scales = "free_x", space = "free")


# Supplementary Figure 18 combined

fig_s18 <- (fig18a) / (cell_prop_bar) +
  plot_annotation(tag_levels = "a")

ggsave("figures/supplementary/Supplementary_Fig18.png", fig_s18, width = 13, height = 12, dpi = 300)


# ============================
#   SUPPLEMENTARY FIGURE 19  #
# ============================

# Made in powerpoint. Figure stored in figures/supplementary

# ============================
#   SUPPLEMENTARY FIGURE 20  #
# ============================

# Finemapping sensitivity analysis

sa <- fread('data/supplementary/sensitivity_analysis_results.txt')
sa$title <- paste0(sa$Gene_name, ', ', sa$Cell, ', ', sa$Chr, ',\n', sa$lead_SNP)

line_df <- sa %>%
  distinct(title, snps_finemapped, multiplier)

finemap_sa <- ggplot(sa, aes(x = multiplier, y = cred_set_size)) +
  geom_point(colour = "#8CD160") +
  geom_hline(data = line_df, aes(yintercept = snps_finemapped), linetype = "dashed", colour = "#608BD1", linewidth = 1) +
  geom_vline(data = line_df, aes(xintercept = 0.2), linetype = "dashed", colour = "#608BD1", linewidth = 1) +
  facet_wrap(~title, scales = "free") +
  labs(
    x = "Value of Prior Variance Multiplier",
    y = "Size of Credible Set"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, vjust = 1.5),
    legend.position = "bottom"
  )

ggsave(filename = "figures/supplementary/Supplementary_Fig20.png", plot = finemap_sa,
       width = 10, height = 8, dpi = 300)

