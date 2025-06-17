# Make manuscript plots
# Created 18/03/2025 by Jennifer Astley <jennifer.astley@balliol.ox.ac.uk
# Last modified 06/05/2025 by Jennifer Astley

# Setup ####

library(pacman)
p_load(tidyverse, png, data.table)

setwd('/gpfs3/well/combat/users/bsg751/projects/eqtl/runtensorqtl/COMBAT_eQTL/figures/main')
# Set chosen working directory
# Directory contains 1. scripts 2. figures
# figures directory contains 1. main and 2. supplementary
# Manuscript directory contains two directories: main and supplementary
# Each contains two directories: figure_data and figures

# ----- Figure 2 ----- ####

# Figure 2a ####

cell_count <- read.table('figure_data/fig2a.csv', header = T)
# Cell | cell_numbers | Count | Mode

scale_factor <- max(cell_count$cell_numbers) / max(cell_count$Count)

png('figures/figure2a.png', width = 10, height = 6, units = "in", res = 300)
p2a <- cell_count %>% 
  mutate(Cell = fct_reorder(Cell, cell_numbers, .desc = TRUE)) %>% 
  ggplot(aes(x = Cell)) +
  geom_bar(aes(y = Count, fill = "eQTL Signal Count"), stat = "identity", position = "dodge") +
  geom_text(aes(y = Count, label = Count), vjust = -0.5) +
  geom_line(aes(y = cell_numbers / scale_factor, color = "Number of Cells Collected"), group = 1, linewidth = 1, linetype = "solid") +
  geom_point(aes(y = cell_numbers / scale_factor, color = "Number of Cells Collected"), size = 3) +
  scale_y_continuous(
    name = "eQTL Signal Count",
    sec.axis = sec_axis(~ . * scale_factor, name = "")
  ) +
  scale_fill_manual(name = "", values = "#3B9AB2", labels = c("eQTL Signal Count")) +
  scale_color_manual(name = "", values = c("Number of Cells Collected" = "red"), labels = c("Number of Cells Collected")) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "italic", size = 15, hjust = 0.5),
    strip.text.x = element_text(size = 15, color = "black", face = "italic"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 70, hjust = 1, size = 12),
    text = element_text(family = "sans")
  ) +
  labs(title = "eQTLs Found Per Cell Type", x = "Cell Type", y = "eQTL Signal Count")

print(p2a)
dev.off()


# Figure 2b ####

# UMAP with number of cells labelled for each cell type

p2b <- readPNG("figures/legend_annotated_clusters.png")
p2b <- rasterGrob(p2b, interpolate=TRUE)

# Figure 2c ####

boxplot_data <- fread('figure_data/eqtl_boxplot_data.csv')
# Genotype | Z_Score | Cell | Gene_name | boxplot

df1 <- boxplot_data %>% filter(boxplot %in% c('boxplot1'))
df2 <- boxplot_data %>% filter(boxplot %in% c('boxplot2'))
df3 <- boxplot_data %>% filter(boxplot %in% c('boxplot3'))
df4 <- boxplot_data %>% filter(boxplot %in% c('boxplot4'))

# P-values, betas, SNPs and positions

# NAPSA in ncMono
pval1 <- 9.64e-22
beta1 <- -1.19
snp1 <- "rs60182980"
snp_pos1 <- "GRCh38: 19: 50,342,447"

# NAPSA in cMono
pval2 <- 8.35e-19
beta2 <- -1.20
snp2 <- "rs60182980"
snp_pos2 <- "GRCh38: 19: 50,342,447"

# ZGLP1 in cMono
pval3 <- 1.56e-9
beta3 <- 1.44
snp3 <- "rs73510898"
snp_pos3 <- "GRCh38: 19: 10,305,768"

# KANSL1 in CD4
pval4 <- 1.47e-14
beta4 <- -1.06
snp4 <- "rs17660167"
snp_pos4 <- "GRCh38: 17: 46,088,945"

df1$gene_cell <- paste0("NAPSA in \n non-classical monocytes")
df2$gene_cell <- paste0("NAPSA in \n classical monocytes")
df3$gene_cell <- paste0("ZGLP1 in \n classical monocytes")
df4$gene_cell <- paste0("KANSL1 in \n CD4 cells")

# eQTL boxplot 1 - NAPSA ncMono

genotype_levels1 <- unique(df1$Genotype)

p2c1 <- df1 %>% 
  mutate(Genotype = factor(Genotype, levels = genotype_levels1)) %>%
  ggplot(aes(x = Genotype, y = Z_Score, fill = Cell)) +
  geom_boxplot(lwd = 1, width = 0.7, outlier.shape = NA) +
  facet_wrap(~ gene_cell, scales = "free_x", ncol = 2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05)) +
  theme_classic() +
  ylim(c(-2.3, 3)) +
  theme(
    plot.title = element_text(face = "bold.italic", size = 15),
    plot.subtitle = element_text(face = "italic", size = 14),
    strip.text.x = element_text(size = 15, color = "black", face = "italic"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.position = "none",
    strip.background = element_rect(fill = "#EBCC2A"),
    strip.placement = "outside"
  ) +
  scale_fill_manual(values = c("#EBCC2A")) +
  geom_text(aes(label = paste0('p = ', pval1, '\n \U03B2 = ', beta1)), 
            x = 3, y = 0.5, hjust = 0.5, vjust = -0.8, size = 4, fontface = 'italic') +
  labs(title = "", x = paste0(snp1, '\n', snp_pos1), y = "Normalised Gene Expression") +
  guides(fill = guide_legend(keywidth = unit(1.5, "cm"), keyheight = unit(1.5, "cm")))

# eQTL boxplot 2

genotype_levels2 <- unique(df2$Genotype)

p2c2 <- df2 %>% 
  mutate(Genotype = factor(Genotype, levels = genotype_levels2)) %>%
  ggplot(aes(x = Genotype, y = Z_Score, fill = Cell)) +
  geom_boxplot(lwd = 1, width = 0.7, outlier.shape = NA) +
  facet_wrap(~ gene_cell, scales = "free_x", ncol = 2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05)) +
  theme_classic() +
  ylim(c(-2.3, 3)) +
  theme(
    plot.title = element_text(face = "bold.italic", size = 15),
    plot.subtitle = element_text(face = "italic", size = 14),
    strip.text.x = element_text(size = 15, color = "black", face = "italic"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.position = "none",
    strip.background = element_rect(fill = "#60BC40"),
    strip.placement = "outside"
  ) +
  scale_fill_manual(values = c("#60BC40")) +
  geom_text(aes(label = paste0('p = ', pval2, '\n \U03B2 = ', beta2)), 
            x = 3, y = 0.5, hjust = 0.5, vjust = -0.8, size = 4, fontface = 'italic') +
  labs(title = "", x = paste0(snp2, '\n', snp_pos2), y = "") +
  guides(fill = guide_legend(keywidth = unit(1.5, "cm"), keyheight = unit(1.5, "cm")))

# eQTL boxplot 3

genotype_levels3 <- unique(df3$Genotype)

p2c3 <- df3 %>% 
  mutate(Genotype = factor(Genotype, levels = genotype_levels3)) %>%
  ggplot(aes(x = Genotype, y = Z_Score, fill = Cell)) +
  geom_boxplot(lwd = 1, width = 0.7, outlier.shape = NA) +
  facet_wrap(~ gene_cell, scales = "free_x", ncol = 2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05)) +
  theme_classic() +
  ylim(c(-2.3, 3)) +
  theme(
    plot.title = element_text(face = "bold.italic", size = 15),
    plot.subtitle = element_text(face = "italic", size = 14),
    strip.text.x = element_text(size = 15, color = "black", face = "italic"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.position = "none",
    strip.background = element_rect(fill = "#60BC40"),
    strip.placement = "outside"
  ) +
  scale_fill_manual(values = c("#60BC40")) +
  geom_text(aes(label = paste0('p = ', pval3, '\n \U03B2 = ', beta3)), 
            x = 3, y = -1.5, hjust = 0.5, vjust = -0.8, size = 4, fontface = 'italic') +
  labs(title = "", x = paste0(snp3, '\n', snp_pos3), y = "") +
  guides(fill = guide_legend(keywidth = unit(1.5, "cm"), keyheight = unit(1.5, "cm")))

# eQTL boxplot 4

genotype_levels4 <- unique(df4$Genotype)

p2c4 <- df4 %>% 
  mutate(Genotype = factor(Genotype, levels = genotype_levels4)) %>%
  ggplot(aes(x = factor(Genotype), y = Z_Score, fill = Cell)) +
  geom_boxplot(lwd = 1, width = 0.7, outlier.shape = NA) +
  facet_wrap(~ gene_cell, scales = "free_x", ncol = 2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05)) +
  theme_classic() +
  ylim(c(-2.3, 3)) +
  geom_text(aes(label = paste0('p = ', pval4, '\n \U03B2 = ', beta4)), 
            x = 1, y = -2.5, hjust = 0.5, vjust = -0.8, size = 4, fontface = 'italic') +
  theme(
    plot.title = element_text(face = "bold.italic", size = 15),
    plot.subtitle = element_text(face = "italic", size = 14),
    strip.text.x = element_text(size = 15, color = "black", face = "italic"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.position = "none",
    strip.background = element_rect(fill = "#5670A6"),
    strip.placement = "outside"
  ) +
  scale_fill_manual(values = c("#5670A6")) +
  labs(title = "", x = paste0(snp4, '\n', snp_pos4), y = "") +
  guides(fill = guide_legend(keywidth = unit(1.5, "cm"), keyheight = unit(1.5, "cm")))

png('figures/figure2c.png', width = 1500, height = 1300, res = 200)
p2c <- grid.arrange(p2c1, p2c2, p2c3, p2c4, nrow = 2)

dev.off()

p2c <- readPNG("figures/figure2c.png")
p2c <- rasterGrob(p2c, interpolate=TRUE)

# Figure 2d ####

# Differential expression of NAPSA

p2d <- readPNG("figures/umap_napsa.png")
p2d <- rasterGrob(p2d, interpolate=TRUE)


# Figure 2 combined ####

png('figures/Figure2.png', height = 2300, width = 2600, res = 200)
figure2 <- grid.arrange(p2a, p2b, p2c, p2d
)
grid.text("(a)", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 16, fontfamily = "Times"))
grid.text("(b)", x = unit(0.52, "npc"), y = unit(0.98, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 16, fontfamily = "Times"))
grid.text("(c)", x = unit(0.02, "npc"), y = unit(0.50, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 16, fontfamily = "Times"))
grid.text("(d)", x = unit(0.52, "npc"), y = unit(0.50, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 16, fontfamily = "Times"))
dev.off()


# ----- Figure 3 ----- ####

combat_gwas_filter <- fread('figure_data/combat_gwas_all_traits_heatmap.csv')
# EQTL_SNP |	GWAS_SNP |	GWAS_TRAIT |	EQTL_CELL_TYPE |	EQTL_GENE |	EQTL_CELL_count
combat_gwas_covid <- fread('figure_data/combat_gwas_covid_traits_heatmap.csv')
# EQTL_SNP |	GWAS_SNP |	GWAS_TRAIT |	STUDY |	EQTL_CELL_TYPE |	EQTL_GENE |	EQTL_CELL_count
onek1k_gwas_filter <- fread('figure_data/onek1k_gwas_all_traits_heatmap.csv')
# EQTL_SNP |	GWAS_SNP |	GWAS_TRAIT |	EQTL_CELL_TYPE |	EQTL_GENE |	EQTL_CELL_count
onek1k_gwas_covid <- fread('figure_data/onek1k_gwas_covid_traits_heatmap.csv')
# EQTL_SNP |	GWAS_SNP |	GWAS_TRAIT |	STUDY |	EQTL_CELL_TYPE |	EQTL_GENE |	EQTL_CELL_count


# Note: these plots are not filtered manually for one study per trait. In the COMBAT case, each COVID-related trait has only 1 study.
# In the ONEK1K case, there are 18 covid-related traits reported across 14 studies.
# All of these are in the heatmaps but this can be modified.

trait_order <- c("Rheumatoid arthritis (rheumatoid factor and/or anti-cyclic citrullinated peptide seropositive)", "Chronic inflammatory diseases (ankylosing spondylitis, Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)", "COVID-19 (critical illness vs population or mild symptoms)", "Hay fever and/or eczema", "Ulcerative colitis", "Alzheimer's disease or family history of Alzheimer's disease", "Amyotrophic lateral sclerosis", "Autoimmune neurological syndromes with anti-GAD65 autoantibodies", "Autoimmune traits", "COVID-19 (hospitalized covid vs population)")

# Figure 3a ####

# COMBAT GWAS all traits (filtered on above studies)
png('figures/figure3a.png', width = 1500, height = 1000)
p3a <- combat_gwas_filter %>%
  mutate(GWAS_TRAIT = factor(GWAS_TRAIT, levels = trait_order)) %>%
  ggplot(aes(x = EQTL_GENE, y = GWAS_TRAIT, fill = EQTL_CELL_count)) +
  geom_tile(color = "white") + 
  scale_fill_gradient(low = "#9D26C6", high = "#C64D26", na.value = "white") +
  labs(title = "Heatmap of eQTL Cell Count - COMBAT All Traits",
       x = "eQTL Gene",
       y = "GWAS Trait",
       fill = "eQTL Cell Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
dev.off()




# Figure 3b ####

top_traits_long <- fread('figure_data/barchart_data_plotting_oct24.csv', sep = '\t')
# trait	hit_combat |	hit_onek1k |	hit_both |	hit_neither |	type |	proportion |	raw_value

# Assuming top_traits_long is already grouped and arranged by hit_combat
png('figures/figure3b.png', width = 1000)
p3b <- top_traits_long %>%
          filter(!type == 'Neither') %>%
          mutate(trait = str_wrap(trait, width = 50),
                 # Set the order of the 'type' variable for stacking
                 type = factor(type, levels = c("OneK1K only", "COMBAT and OneK1K", "COMBAT only")),
                 # Set the order of traits based on the value of hit_combat
                 trait = factor(trait, levels = unique(trait[order(-hit_combat)]))) %>%  # Use hit_combat to determine order
          ggplot(aes(x = raw_value, y = trait, fill = type)) + 
          geom_bar(stat = "identity", position = "stack") +
          labs(title = "GWAS Variants Found in eQTL Datasets by Trait",
               y = "Trait",
               x = "GWAS Variants in eQTL Datasets") +
          geom_text(aes(label = ifelse(raw_value == 0, "", raw_value)), 
                    position = position_stack(vjust = 0.5), 
                    size = 3.5, color = "white") +
          theme_minimal() +
          theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
          scale_fill_manual(values = c(
            "COMBAT only" = "#60BC40",
            "OneK1K only" = "#2683C6",
            "COMBAT and OneK1K" = "#C62683"
          ))
dev.off()


# Figure 3c ####
# eQTL boxplots for GWAS intersection genes

whole_result <- read.table('figure_data/gwas_gene_eqtl_boxplot_data.csv', header = T)
# Genotype |	Z_Score |	Cell |	Gene_name |	boxplot

df1 <- whole_result %>% filter(boxplot %in% c('boxplot1'))
df2 <- whole_result %>% filter(boxplot %in% c('boxplot2'))
df3 <- whole_result %>% filter(boxplot %in% c('boxplot3'))

df1$gene_cell <- "REL in \n CD4 cells"
df2$gene_cell <- "IRF5 in \n non-classical monocytes"
df3$gene_cell <- "TRAF1 in \n CD8 cells"

snp1 <- "rs13017599"
snp_pos1 <- "60,937,196"
chr1 <- 2

snp2 <- "rs3778754"
snp_pos2 <- "128,935,498"
chr2 <- 7

snp3 <- "rs2159778"
snp_pos3 <- "120,915,851"
chr3 <- 9

# eQTL boxplot 1

genotype_levels1 <- unique(df1$Genotype)
p3c1 <- df1 %>% 
  mutate(Genotype = factor(Genotype, levels = genotype_levels1)) %>%
  ggplot(aes(x = Genotype, y = Z_Score, fill = Cell)) +
  geom_boxplot(lwd = 1, width = 0.7, outlier.shape = NA) +
  facet_wrap(~ gene_cell, scales = "free_x", ncol = 2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05)) +
  theme_minimal() +
  ylim(c(-2.3, 3)) +
  theme(
    plot.title = element_text(face = "bold.italic", size = 11),
    plot.subtitle = element_text(face = "italic", size = 11),
    strip.text.x = element_text(size = 11, color = "black", face = "italic"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11),
    legend.position = "none",
    strip.background = element_rect(fill = '#df65b0'),
    strip.placement = "outside"
  ) +
  scale_fill_manual(values = c('#df65b0')) + 
  labs(title = "", x = paste0(snp1, '\n GRCh38: ', chr1, ':', snp_pos1), y = "Z Score") +
  guides(fill = guide_legend(keywidth = unit(1.5, "cm"), keyheight = unit(1.5, "cm")))

png('figures/figure3c1.png', width = 300, height = 250)
p3c1
dev.off()


# eQTL boxplot 2

genotype_levels2 <- unique(df2$Genotype)

p3c2 <- df2 %>% 
  mutate(Genotype = factor(Genotype, levels = genotype_levels2)) %>%
  ggplot(aes(x = Genotype, y = Z_Score, fill = Cell)) +
  geom_boxplot(lwd = 1, width = 0.7, outlier.shape = NA) +
  facet_wrap(~ gene_cell, scales = "free_x", ncol = 2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05)) +
  theme_minimal() +
  ylim(c(-2.3, 3)) +
  theme(
    plot.title = element_text(face = "bold.italic", size = 11),
    plot.subtitle = element_text(face = "italic", size = 11),
    strip.text.x = element_text(size = 11, color = "black", face = "italic"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11),
    legend.position = "none",
    strip.background = element_rect(fill = '#7fcdbb'),
    strip.placement = "outside"
  ) +
  scale_fill_manual(values = c('#7fcdbb')) +
  labs(title = "", x = paste0(snp2, '\n GRCh38: ', chr2, ':', snp_pos2), y = "") +
  guides(fill = guide_legend(keywidth = unit(1.5, "cm"), keyheight = unit(1.5, "cm")))

png('figures/figure3c2.png', width = 300, height = 250)
p3c2
dev.off()

# eQTL boxplot 3

genotype_levels3 <- unique(df3$Genotype)

p3c3 <- df3 %>% 
  mutate(Genotype = factor(Genotype, levels = genotype_levels3)) %>%
  ggplot(aes(x = Genotype, y = Z_Score, fill = Cell)) +
  geom_boxplot(lwd = 1, width = 0.7, outlier.shape = NA) +
  facet_wrap(~ gene_cell, scales = "free_x", ncol = 2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05)) +
  theme_minimal() +
  ylim(c(-2.3, 3)) +
  theme(
    plot.title = element_text(face = "bold.italic", size = 11),
    plot.subtitle = element_text(face = "italic", size = 11),
    strip.text.x = element_text(size = 11, color = "black", face = "italic"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11),
    legend.position = "none",
    #legend.text = element_text(size = 12),
    strip.background = element_rect(fill = '#fe9929'),
    strip.placement = "outside"
  ) +
  scale_fill_manual(values = c('#fe9929')) +
  labs(title = "", x = paste0(snp3, '\n GRCh38: ', chr3, ':', snp_pos3), y = "") +
  guides(fill = guide_legend(keywidth = unit(1.5, "cm"), keyheight = unit(1.5, "cm")))

png('figures/figure3c3.png', width = 300, height = 250)
p3c3
dev.off()


png('figures/figure3c.png', width = 800, height = 250)
p3c <- grid.arrange(p3c1, p3c2, p3c3, nrow = 1)
dev.off()

# Figure 3d ####
# Pseudobulk expression violin plots by infection severity

final_data <- read.table('figure_data/gwas_gene_violin_plot_data.csv', header = T)
# Z_score | COMBAT_ID | Source | ventilation_assistance | COVID_DIAG | DIAG | INFECTION_SEVERITY | Cell_Type | Gene | infection

rel_colors <- c("Healthy" = "#f1eef6", "Asymptomatic" = "#d4b9da", "Mild" = "#df65b0", 
                "Severe" = "#e31a1c", "Critical" = "#67001f") 

irf5_colors <- c("Healthy" = "#f7fcf0", "Asymptomatic" = "#c7e9b4", "Mild" = "#7fcdbb", 
                 "Severe" = "#2c7fb8", "Critical" = "#253494")

traf1_colors <- c("Healthy" = "#fff7bc", "Asymptomatic" = "#fec44f", "Mild" = "#fe9929", 
                  "Severe" = "#d95f0e", "Critical" = "#993404")


# Assign colors dynamically within final_data
final_data <- final_data %>%
  mutate(fill_color = case_when(
    Gene == "REL"  ~ rel_colors[infection],
    Gene == "IRF5" ~ irf5_colors[infection],
    Gene == "TRAF1" ~ traf1_colors[infection]
  ))

# Ensure the facets appear in the desired order
final_data <- final_data %>%
  mutate(Gene = factor(Gene, levels = c("REL", "IRF5", "TRAF1")),
         Cell_Type = factor(Cell_Type, levels = c("CD4", "ncMono", "CD8")),
         infection = factor(infection, levels = c("Healthy", "Asymptomatic", "Mild", "Severe", "Critical")))

# Create the violin plot with facet_wrap and correctly applied colors
png('figures/figure3d.png', width = 800, height = 250)
p3d <- ggplot(final_data, aes(x = infection, y = Z_score, fill = fill_color)) +
  geom_violin() +
  labs(
    x = "Infection Severity",
    y = "Z Score",
    fill = ""
  ) +
  facet_wrap(~ Cell_Type, scales = "free_y", nrow = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    strip.text = element_blank()
  ) +
  scale_fill_identity()
dev.off()


# Figure 3e ####
# UMAP expression of three GWAS intersection genes

p3e1 <- readPNG("figures/umap_REL.png")
p3e1 <- rasterGrob(p3e1, interpolate=TRUE)

p3e2 <- readPNG("figures/umap_IRF5.png")
p3e2 <- rasterGrob(p3e2, interpolate=TRUE)

p3e3 <- readPNG("figures/umap_TRAF1.png")
p3e3 <- rasterGrob(p3e3, interpolate=TRUE)

p3e <- grid.arrange(p3e1, p3e2, p3e3, nrow = 1)

png("figures/figure_3e.png", width = 3600, height = 1200, res = 300)
p3e
dev.off()


# Figure 3 combined ####

png("figures/Figure3.png", width = 4200, height = 6000, res = 300)
p3 <- grid.arrange(p3a, p3b, p3c, p3d, p3e, nrow = 5, heights = c(1.5, 1.5, 1, 1, 1))
grid.text("(a)", x = unit(0.02, "npc"), y = unit(0.99, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 16, fontfamily = "Times"))
grid.text("(b)", x = unit(0.02, "npc"), y = unit(0.75, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 16, fontfamily = "Times"))
grid.text("(c)", x = unit(0.02, "npc"), y = unit(0.50, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 16, fontfamily = "Times"))
grid.text("(d)", x = unit(0.02, "npc"), y = unit(0.34, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 16, fontfamily = "Times"))
grid.text("(e)", x = unit(0.02, "npc"), y = unit(0.15, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 16, fontfamily = "Times"))
dev.off()

# ----- Figure 4 ----- ####

# Figure 4a ####

grid_score <- read.table('figure_data/grid_score_2_cells_jun_24.csv', header = T)
# Cell |	INFECTION_SEVERITY |	pval |	beta |	pval_BH
boxplt_result <- read.table('figure_data/context_specific_boxplot_data_rps26_adam10.csv', header = T)
# Genotype |	Z_Score	DIAG |	INFECTION_SEVERITY |	HOSP |	Cell |	Gene

# RPS26
gene <- 'ENSG00000197728' 
gene_name <- "RPS26"
snp <- "rs10876864"
chr <- 12



order_levels <- c('Healthy\ncontrols', 'Community\ncases', 'Mild', 'Severe', 'Critical\n(Intubated)')

grid_score <- grid_score %>%
  mutate(Infection_Label = case_when(
    INFECTION_SEVERITY == 0 ~ 'Healthy\ncontrols',
    INFECTION_SEVERITY == 1 ~ 'Community\ncases',
    INFECTION_SEVERITY == 2 ~ 'Mild',
    INFECTION_SEVERITY == 3 ~ 'Severe',
    INFECTION_SEVERITY == 4 ~ 'Critical\n(Intubated)',
    TRUE ~ 'Unknown'
  )) %>% mutate(Infection_Label = factor(Infection_Label, levels = order_levels))

png('figures/figure4a.png', width = 450, height = 650)
grid_score %>%  
  ggplot(aes(Cell, Infection_Label, fill = beta)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  theme_classic(base_size = 20) +
  labs(fill = "eQTL\neffect size") +
  scale_fill_gradientn(colors = c("#00A08A", "#FFFC46", "#F21A00"), limits = c(-2.3, -0.5)) +
  labs(x = "", y = "") +  
  theme(
    plot.title = element_text(face = "bold.italic", hjust = 0.5, vjust = 1.5, size = 20),
    legend.position = "bottom",
    legend.box="horizontal",
    legend.key.size = unit(1.2, "cm")
  ) +
  ggtitle("RPS26") +
  annotate("text", x = 2.75, y = 1.5, label = "Non-hospitalised", size = 7, color = "#60BC40", angle = 90) +
  annotate("text", x = 2.75, y = 4, label = "Hospitalised", size = 7, color = "#F21A00", angle = 90) +
  annotate("text", x = 3.1, y = 4, label = "dummy line", size = 6, color = "white", angle = 90) +
  geom_segment(aes(x = 2.6, y = 0.5, xend = 2.6, yend = 2.5), color = "#60BC40", linewidth = 3, linetype = "solid") +
  geom_segment(aes(x = 2.6, y = 2.6, xend = 2.6, yend = 5.5), color = "#F21A00", linewidth = 3, linetype = "solid") +
  guides(col = guide_legend(title.position = "top",title.hjust =0.5))
dev.off()

p4a <- readPNG("figures/figure4a.png")
p4a <- rasterGrob(p4a, interpolate=TRUE)


# Figure 4b ####


# First plot (RPS26, 2 cell types)
p4b1 <- boxplt_result %>% 
  filter(! DIAG %in% c('HEALTHY')) %>% 
  filter(Cell %in% c('CD8', 'DC')) %>% 
  mutate(Diag = factor(DIAG)) %>%
  ggplot(aes(x = Diag, y = Z_Score, colour = Genotype)) +
  geom_point(position = position_jitter(width = 0.1), size = 4, shape = 19) + 
  facet_wrap(~Cell, scales='free_y', nrow=1)+
  stat_summary(fun = mean, geom = "point", shape = 17, size = 4, position = position_dodge(width = 0.2)) +
  geom_line(stat = "summary", fun = mean, aes(group = Genotype), position = position_dodge(width = 0.2), linewidth=1.5) +
  theme_classic()+
  theme(
    plot.title = element_text(face = "bold.italic", size = 15, hjust = 0.5),
    plot.subtitle = element_text(face = "italic", size = 14),
    plot.caption = element_text(face = "italic", size = 10, hjust = 0.5),
    strip.text.x = element_text(size = 15, color = "black", face = "italic"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14)
  ) +
  scale_colour_manual(values = c("#5BBCD6", "#C641ED", "#046C9A")) +
  labs(title = "RPS26", fill = "", x = "rs10876864\nGRCh38: 12:56,007,301", y = expression(paste("Normalised ", italic("RPS26"), " Expression"))) +
  guides(color = guide_legend(override.aes = list(size = 4)))

# Second plot (ADAM10, 1 cell type)
p4b2 <- boxplt_result %>% 
  filter(! DIAG %in% c('HEALTHY')) %>% 
  filter(Cell %in% c('CD4')) %>% 
  mutate(Diag = factor(DIAG)) %>%
  ggplot(aes(x = Diag, y = Z_Score, colour = Genotype)) +
  geom_point(position = position_jitter(width = 0.1), size = 4, shape = 19) + 
  facet_wrap(~Cell, scales='free_y', nrow=1)+
  stat_summary(fun = mean, geom = "point", shape = 19, size = 4, position = position_dodge(width = 0.2)) +
  geom_line(stat = "summary", fun = mean, aes(group = Genotype), position = position_dodge(width = 0.2), linewidth=1.5) +
  theme_classic()+
  theme(
    plot.title = element_text(face = "bold.italic", size = 15, hjust = 0.5),
    plot.subtitle = element_text(face = "italic", size = 14),
    plot.caption = element_text(face = "italic", size = 10, hjust = 0.5),
    strip.text.x = element_text(size = 15, color = "black", face = "italic"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14)
  ) +
  scale_colour_manual(values = c("#5BBCD6", "#C641ED", "#046C9A")) +
  labs(title = "ADAM10", fill = "", x = "rs7161799\nGRCh38: 15:58,478,324", y = expression(paste("Normalised ", italic("ADAM10"), " Expression"))) +
  guides(color = guide_legend(override.aes = list(size = 4)))

# Figure 4c ####

leg_p4c1 <- boxplt_result %>%
  filter(Cell %in% c('NK', 'cMono')) %>% 
  mutate(HOSP = factor(HOSP, levels = c('Non-hospitalised', 'Hospitalised'))) %>%
  ggplot(aes(x = Genotype, y = Z_Score, fill = HOSP)) +
  geom_boxplot(lwd = 1, width = 0.7, outlier.shape = NA) + 
  facet_wrap(~Cell, scales = 'free_y', ncol = 3) +
  theme_classic() +
  theme(
    strip.text.x = element_text(size = 15, color = "black", face = "italic"),
    legend.position = "right", 
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22), 
    legend.key.size = unit(2.5, "lines")
  ) +
  scale_fill_manual(values = c("#60BC40", "#F21A00")) +
  labs(fill = "")

legend_p4c1 <- get_legend(leg_p4c1)

p4c1_plot <- boxplt_result %>%
  filter(Cell %in% c('NK', 'cMono')) %>% 
  mutate(HOSP = factor(HOSP, levels = c('Non-hospitalised', 'Hospitalised'))) %>%
  ggplot(aes(x = Genotype, y = Z_Score, fill = HOSP)) +
  geom_boxplot(lwd = 1, width = 0.7, outlier.shape = NA) +
  facet_wrap(~Cell, scales = 'free_y', ncol = 1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05)) +
  geom_smooth(aes(group = HOSP, linetype = HOSP, color = HOSP), method = "lm", formula = y ~ x, 
              se = F, linewidth = 1.2, linetype = "dashed") +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold.italic", size = 15, hjust = 0.5),
    plot.subtitle = element_text(face = "italic", size = 14),
    plot.caption = element_text(face = "italic", size = 10, hjust = 0.5),
    strip.text.x = element_text(size = 15, color = "black", face = "italic"),
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14)
  ) +
  scale_fill_manual(values = c("#60BC40", "#F21A00")) +
  scale_color_manual(values = c("#60BC40", "#F21A00")) +
  labs(title = "", fill = "", x = snp, y = expression(paste("Normalised ", italic("RPS26"), " Expression")))

leg_p4c2 <- boxplt_result %>%
  filter(! DIAG %in% c('HEALTHY')) %>% 
  filter(Cell %in% c('CD8', 'DC')) %>% 
  mutate(DIAG = factor(DIAG)) %>%
  ggplot(aes(x = Genotype, y = Z_Score, fill = DIAG)) +
  geom_boxplot(lwd = 1, width = 0.7, outlier.shape = NA) +
  facet_wrap(~Cell, scales = 'free_y', ncol = 3) +
  theme_classic() +
  theme(
    strip.text.x = element_text(size = 15, color = "black", face = "italic"),
    legend.position = "right",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.key.size = unit(2.5, "lines")
  ) +
  scale_fill_manual(values = c("#C62683", "#2683C6")) +
  labs(fill = "")

legend_p4c2 <- get_legend(leg_p4c2)

p4c2_plot <- boxplt_result %>%
  filter(! DIAG %in% c('HEALTHY')) %>% 
  filter(Cell %in% c('CD8', 'DC')) %>% 
  mutate(DIAG = factor(DIAG)) %>%
  ggplot(aes(x = Genotype, y = Z_Score, fill = DIAG)) +
  geom_boxplot(lwd = 1, width = 0.7, outlier.shape = NA) +
  facet_wrap(~Cell, scales = 'free_y', ncol = 1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05)) +
  geom_smooth(aes(group = DIAG, linetype = DIAG, color = DIAG), method = "lm", formula = y ~ x, 
              se = F, linewidth = 1.2, linetype = "dashed") +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold.italic", size = 15, hjust = 0.5),
    plot.subtitle = element_text(face = "italic", size = 14),
    plot.caption = element_text(face = "italic", size = 10, hjust = 0.5),
    strip.text.x = element_text(size = 15, color = "black", face = "italic"),
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14)
  ) +
  scale_fill_manual(values = c("#C62683", "#2683C6")) +
  scale_color_manual(values = c("#C62683", "#2683C6")) +
  labs(title = "", fill = "", x = snp, y = expression(paste("Normalised ", italic("RPS26"), " Expression")))

# Figure 4 combined ####

fig4_matrix_array <- rbind(c(1,1,1,1,2,2,2,2,3,3),
                           c(1,1,1,1,2,2,2,2,3,3),
                           c(1,1,1,1,2,2,2,2,3,3),
                           c(1,1,1,1,4,4,4,5,5,5),
                           c(6,6,6,6,6,7,7,7,7,7),
                           c(6,6,6,6,6,7,7,7,7,7),
                           c(6,6,6,6,6,7,7,7,7,7)
)

png('figures/Figure4.png', height = 1400, width = 1200)
figure4 <- grid.arrange(p4a, p4b1, p4b2, legend_p4c1, legend_p4c2, p4c1_plot, p4c2_plot,
                        nrow = 7,
                        heights = c(1, 1, 1, 0.5, 1.2, 1.2, 1.2),
                        layout_matrix = fig4_matrix_array
)
grid.text("(a)", x = unit(0.02, "npc"), y = unit(0.99, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 16, fontfamily = "Times"))
grid.text("(b)", x = unit(0.4, "npc"), y = unit(0.99, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 16, fontfamily = "Times"))
grid.text("(c)", x = unit(0.02, "npc"), y = unit(0.52, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 16, fontfamily = "Times"))
print(figure4)
dev.off()




# Supplementary Figures ####
setwd("/gpfs3/well/combat/users/bsg751/projects/eqtl/runtensorqtl/COMBAT_eQTL/figures/supplementary/")

# Figure S1 ####
# Intersection of eGenes per cell type (heatmap)
# Heatmap showing overlapping eGenes between cell types (major resolution)

egene_overlap <- fread('figure_data/egene_overlap.csv')
# Mode | Cell1 | Cell2 | overlap_percent | overlap_raw

cell_vector <- c("CD4", "cMono", "CD8", "NK", "ncMono", "B", "DC", "DP", "GDT", "HSC", "DN", "MAIT", "PB", "PLT", "iNKT")
egene_overlap$Cell1 <- factor(egene_overlap$Cell1, levels=cell_vector)
egene_overlap$Cell2 <- factor(egene_overlap$Cell2, levels=cell_vector)


png("figures/figure_S1.png", width = 11, height = 9, units = "in", res = 300)
fs1 <- egene_overlap %>% 
  filter(overlap_percent != 100) %>% 
  ggplot(aes(Cell1, Cell2, fill = overlap_percent)) +
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, hjust=1)) +
  geom_text(aes(label = overlap_raw), color = "white", size = 4) +
  theme(
    strip.text.x = element_text(
      size = 10, color = "black", face = "italic"
    )
  ) +
  theme(plot.title = element_text(face = "italic", size=15)) +
  theme_classic(base_size = 16) +
  labs(fill = "Overlapping \neGenes (%)") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 12),
        text = element_text(family = "sans")) +
  scale_fill_gradientn(colors = c("#78B7C5", "#273046")) +
  labs(x = "", y = "", title = "")
print(fs1)
dev.off()

# Figure S2 ####

# Intersection of eGenes per cell type (upset plot)

# Upset plot showing eGene cell type intersection/overlap
# Filtering intersections containing only one Gene-SNP pair for image clarity

cis_results <- fread('figure_data/all_results_cis_indep.csv')
# Mode |	Cell |	Chr |	Gene |	SNP |	SNP_TSS |	maf |	pq_val |	slope |	slope_se

filtered_upset_data <- new_dataframe %>%
  distinct(Gene, SNP, .keep_all = TRUE) %>%
  unnest(cols = Cell_List) %>%
  mutate(CellMember = 1) %>%
  pivot_wider(
    names_from = Cell_List,
    values_from = CellMember,
    values_fill = list(CellMember = 0)
  )

filtered_upset_data <- filtered_upset_data %>%
  mutate(intersection_key = apply(select(., CD4, cMono, CD8, NK, ncMono, B, DC, DP, GDT, HSC, DN, MAIT, PB, PLT, iNKT), 1, paste, collapse = "")) %>%
  group_by(intersection_key) %>%
  filter(n() > 1) %>%   # Keep only combinations that appear more than once
  ungroup()

# Plot
png("figures/figure_S2.png", width = 22, height = 9, units = "in", res = 300)
UpSetR::upset(
  as.data.frame(filtered_upset_data),
  sets = c("CD4", "cMono", "CD8", "NK", "ncMono", "B", "DC", "DP", "GDT", "HSC", "DN", "MAIT", "PB", "PLT", "iNKT"),
  keep.order = FALSE,
  nintersects = 90,
  text.scale = 1.5,
  set_size.show = TRUE,
  sets.bar.color = "#046C9A",
  main.bar.color = "#046C9A",
  matrix.color = "#046C9A"
)
dev.off()




# Figure S3 ####

# Differential expression of ZGLP1, KANSL1, ADAM10, RPS26
# Made in python

# Figure S4 ####

# ONEK1K GWAS all traits (filtered on above studies)
png("figures/figure_S4.png", width = 24, height = 9, units = "in", res = 300)
onek1k_gwas_filter %>%
  mutate(GWAS_TRAIT = factor(GWAS_TRAIT, levels = trait_order)) %>%
  ggplot(aes(x = EQTL_GENE, y = GWAS_TRAIT, fill = EQTL_CELL_count)) +
  geom_tile(color = "white") + 
  scale_fill_gradient(low = "#269FC6", high = "#4FC626", na.value = "white") +
  labs(title = "Heatmap of eQTL Cell Count - OneK1K All Traits",
       x = "eQTL Gene",
       y = "GWAS Trait",
       fill = "eQTL Cell Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
dev.off()

# Figure S5 ####

# COMBAT GWAS covid traits
png("figures/figure_S5.png", width = 10, height = 5, units = "in", res = 300)
combat_gwas_covid %>%
  ggplot(aes(x = EQTL_GENE, y = GWAS_TRAIT, fill = EQTL_CELL_count)) +
  geom_tile(color = "white") + 
  scale_fill_gradient(
    low = "#9D26C6",
    high = "#C64D26",
    na.value = "white",
    breaks = c(1, 2, 3),
    limits = c(1, 3)
  ) +
  labs(title = "Heatmap of eQTL Cell Count - COMBAT COVID Traits",
       x = "eQTL Gene",
       y = "GWAS Trait",
       fill = "eQTL Cell Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
dev.off()

# Figure S6 ####

# ONEK1K GWAS covid traits
png("figures/figure_S6.png", width = 22, height = 9, units = "in", res = 300)
onek1k_gwas_covid %>%
  ggplot(aes(x = EQTL_GENE, y = GWAS_TRAIT, fill = EQTL_CELL_count)) +
  geom_tile(color = "white") + 
  scale_fill_gradient(
    low = "#269FC6",
    high = "#4FC626",
    na.value = "white",
    breaks = seq(1, 13, by = 2),
    limits = c(1, 13)
  ) +
  labs(title = "Heatmap of eQTL Cell Count - ONEK1K COVID Traits",
       x = "eQTL Gene",
       y = "GWAS Trait",
       fill = "eQTL Cell Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
dev.off()

# Figure S7 ####

# Genetic principal component analysis

genetic_pcs <- fread('figure_data/PCA_assignments.csv')

png("figures/figure_S7.png", width = 6, height = 4, units = "in", res = 300)
ggplot(genetic_pcs, aes(x = PC1, y = PC2)) +
  geom_point(
    data = subset(data_combined, group == "1000G_EUR"),
    aes(color = Population),
    shape = 19, size = 2
  ) +
  geom_point(
    data = subset(data_combined, group == "1000G_nonEUR"),
    aes(x = PC1, y = PC2),
    color = "lightgrey",
    shape = 19, size = 2
  ) +
  geom_point(
    data = subset(data_combined, group == "COMBAT_EUR"),
    aes(x = PC1, y = PC2),
    color = "black", fill = "black",
    shape = 21, size = 2, stroke = 0.5
  ) +
  geom_point(
    data = subset(data_combined, group == "COMBAT_nonEUR"),
    aes(x = PC1, y = PC2),
    color = "black", fill = NA,
    shape = 21, size = 2, stroke = 0.5
  ) +
  theme_minimal() +
  labs(
    title = "",
    x = "PC1", y = "PC2",
    color = "Population\n(1000G EUR)"
  ) +
  theme(legend.position = "right")
dev.off()

# Figure S8 ####

# Imputation quality
# Made by others

# Figure S9 ####

# Clinical phenotype categorisation
# Made in PowerPoint

# Figure S10 ####

# Sensitivity analysis of the effect of prior variance multiplier on credible set size

cred_set_data <- fread('figure_data/finemapping_sensitivity_analysis.csv')
# Cell | Chr |  Gene |   lead_SNP | multiplier | snps_finemapped | cred_set_size

png('figures/figure_S10.png', height=1600, width=1600, res=150)
credset_plot <- cred_set_data %>% 
  ggplot(aes(x=multiplier)) + 
  facet_wrap(~Gene + Cell + Chr + lead_SNP, labeller = label_wrap_gen(multi_line=FALSE)) + 
  geom_point(aes(y=cred_set_size), size=1.5, color="#60BC40") +
  geom_line(aes(y=snps_finemapped, color="variants finemapped"), linetype='dashed', linewidth=1.) +
  theme_classic() +
  scale_color_manual(name = "", values = "#273046") +
  theme(plot.title = element_text(face = "bold.italic", size=15, hjust=0.5),
        plot.subtitle = element_text(face = "italic", size=14),
        plot.caption = element_text(face = "italic", size=10, hjust=0.5),
        strip.text.x = element_text(size = 15, color = "black", face = "italic"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 15),
        legend.position = 'bottom') + 
  labs(x = "Value of Prior Variance Multiplier", y = "Size of Credible Set")
print(credset_plot)
dev.off()

