rm(list = ls())

############################################################
# CONFIGURATION: edit these paths before running
############################################################
# Run from the repository root so that "data/" resolves.
DATA_DIR <- "data"                   # files included in this repository
EXT_DIR  <- "/YOURDIRECTORY"         # inputs not included here (see README)
OUT_DIR  <- "figures"                # figures are written here
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
MAG_SUMMARY_FILE <- file.path(DATA_DIR, "MAGs_aridity_proteome_summary.csv")
BAC120_SUMMARY_FILE <- file.path(DATA_DIR, "bac120_common69_per_genome_summary.csv")

library(data.table)
library(arrow)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(stringr)
library(tidyr)
library(dplyr)

##########################################################################################
# Analysis of the amino acid composition of bac120 marker genes only for global bacterial MAGs with aridity 
mag.dbb = read.csv(MAG_SUMMARY_FILE, header = T)
head(mag.dbb)
bac120.pepstats = read.csv(BAC120_SUMMARY_FILE, header = T)
bac120.pepstats$Proteome <- sub(".*(GC[AF]_[0-9]+).*", "\\1", bac120.pepstats$sample)
head(bac120.pepstats) ## WE CHOSE THE 69 BAC120 GENES THAT APPEAR IN AT LEAST 90% OF THE MAGs

bac120.dbb = merge(mag.dbb, bac120.pepstats, by = "Proteome")
head(bac120.dbb)

bac120.dbb <- bac120.dbb %>%
  # Categorize aridity values into groups
  filter(Aridity >= 0, Aridity < 6) %>%
  mutate(Aridity_group = case_when(
    Aridity >= 0    & Aridity < 0.2  ~ "Arid",
    Aridity >= 0.2  & Aridity < 0.65 ~ "Semiarid",
    Aridity >= 0.65 & Aridity < 6    ~ "Humid"
  )) %>%
  mutate(Aridity_group = factor(Aridity_group, levels = c("Arid", "Semiarid", "Humid"))) %>%
  filter(!is.na(Aridity_group))
head(bac120.dbb)

COLS_ARID <- c(
  "Arid"      = "#e31a1c",
  "Semiarid"  = "#f89540",
  "Humid"     = "#666666"
)

# Acidic AAs
acidic.boxplot <- ggplot(bac120.dbb, 
                         aes(x = Aridity_group, y = mole_pct_acidic)) + geom_boxplot(aes(fill = Aridity_group), 
               width = 0.55, 
               alpha = 0.6,
               color = "grey20", 
               outlier.shape = NA) +
  geom_jitter(aes(color = Aridity_group), 
              width = 0.18, 
              size = 1.0,
              alpha = 0.8, 
              show.legend = FALSE) +
  scale_fill_manual(values = COLS_ARID, guide = "none") +
  scale_color_manual(values = COLS_ARID, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  ggsignif::geom_signif(
    comparisons = list(c("Arid", "Semiarid"), 
                       c("Arid", "Humid"), 
                       c("Semiarid", "Humid")),
    test = "wilcox.test",
    map_signif_level = FALSE,   
    textsize = 3.5,
    step_increase = 0.065,       
    tip_length = 0.01,
    vjust = 0.12                 
  ) +
  labs(
    x = "Habitat origin",
    y = "Proportion of acidic amino acids (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14))
acidic.boxplot

tiff(file.path(OUT_DIR, "bac120_acidic_aridity_boxplot_mags.tiff"), units="in", width=4, height=4, res=300)
acidic.boxplot
dev.off()


# Basic AAs
basic.boxplot <- ggplot(bac120.dbb, 
                         aes(x = Aridity_group, y = mole_pct_basic)) + geom_boxplot(aes(fill = Aridity_group), 
                                                                                     width = 0.55, 
                                                                                     alpha = 0.6,
                                                                                     color = "grey20", 
                                                                                     outlier.shape = NA) +
  geom_jitter(aes(color = Aridity_group), 
              width = 0.18, 
              size = 1.0,
              alpha = 0.8, 
              show.legend = FALSE) +
  scale_fill_manual(values = COLS_ARID, guide = "none") +
  scale_color_manual(values = COLS_ARID, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  ggsignif::geom_signif(
    comparisons = list(c("Arid", "Semiarid"), 
                       c("Arid", "Humid"), 
                       c("Semiarid", "Humid")),
    test = "wilcox.test",
    map_signif_level = FALSE,   
    textsize = 3.5,
    step_increase = 0.065,       
    tip_length = 0.01,
    vjust = 0.12                 
  ) +
  labs(
    x = "Habitat origin",
    y = "Proportion of basic amino acids (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14))
basic.boxplot

tiff(file.path(OUT_DIR, "bac120_basic_aridity_boxplot_mags.tiff"), units="in", width=4, height=4, res=300)
basic.boxplot
dev.off()


# Charged AA
charged.boxplot <- ggplot(bac120.dbb, 
                         aes(x = Aridity_group, y = mole_pct_charged)) + geom_boxplot(aes(fill = Aridity_group), 
                                                                                     width = 0.55, 
                                                                                     alpha = 0.6,
                                                                                     color = "grey20", 
                                                                                     outlier.shape = NA) +
  geom_jitter(aes(color = Aridity_group), 
              width = 0.18, 
              size = 1.0,
              alpha = 0.8, 
              show.legend = FALSE) +
  scale_fill_manual(values = COLS_ARID, guide = "none") +
  scale_color_manual(values = COLS_ARID, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  ggsignif::geom_signif(
    comparisons = list(c("Arid", "Semiarid"), 
                       c("Arid", "Humid"), 
                       c("Semiarid", "Humid")),
    test = "wilcox.test",
    map_signif_level = FALSE,   
    textsize = 3.5,
    step_increase = 0.065,       
    tip_length = 0.01,
    vjust = 0.12                 
  ) +
  labs(
    x = "Habitat origin",
    y = "Proportion of charged amino acids (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14))
charged.boxplot

tiff(file.path(OUT_DIR, "bac120_charged_aridity_boxplot_mags.tiff"), units="in", width=4, height=4, res=300)
charged.boxplot
dev.off()


# Isoelectric point (pI)
ip.boxplot <- ggplot(bac120.dbb, 
                          aes(x = Aridity_group, y = Ip)) + geom_boxplot(aes(fill = Aridity_group), 
                                                                                       width = 0.55, 
                                                                                       alpha = 0.6,
                                                                                       color = "grey20", 
                                                                                       outlier.shape = NA) +
  geom_jitter(aes(color = Aridity_group), 
              width = 0.18, 
              size = 1.0,
              alpha = 0.8, 
              show.legend = FALSE) +
  scale_fill_manual(values = COLS_ARID, guide = "none") +
  scale_color_manual(values = COLS_ARID, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  ggsignif::geom_signif(
    comparisons = list(c("Arid", "Semiarid"), 
                       c("Arid", "Humid"), 
                       c("Semiarid", "Humid")),
    test = "wilcox.test",
    map_signif_level = FALSE,   
    textsize = 3.5,
    step_increase = 0.065,       
    tip_length = 0.01,
    vjust = 0.12                 
  ) +
  labs(
    x = "Habitat origin",
    y = "Isoelectric point (Ip)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14))
ip.boxplot

tiff(file.path(OUT_DIR, "bac120_ip_aridity_boxplot_mags.tiff"), units="in", width=4, height=4, res=300)
ip.boxplot
dev.off()

