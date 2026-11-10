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
MAG_TREE_FILE <- file.path(DATA_DIR, "MAGs_bac120_GTDB_phylogeny.tree")
MAG_PI_TIPS_FILE <- file.path(DATA_DIR, "MAGs_pI_phylogeny_tips.csv")

############################################################
# LIBRARIES
############################################################
library(terra)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(rnaturalearth)
library(rnaturalearthdata)
library(arrow)
library(purrr)
library(data.table)
library(ggsignif)
library(ggpubr)

df_summarized_complete = read.csv(MAG_SUMMARY_FILE,
                                  header = T)
head(df_summarized_complete)

pi.model = lm(Aridity ~ Ip, data = df_summarized_complete)
AIC(pi.model) # AIC = 290.69
acid.percent.model = lm(Aridity ~ acid.aa.percent, data = df_summarized_complete)
AIC(acid.percent.model) # AIC = 285.13 BEST MODEL FOR AA-SEQUENCE LEVEL ANALYSIS***
basic.percent.model = lm(Aridity ~ basic.aa.percent, data = df_summarized_complete)
AIC(basic.percent.model) # AIC = 331.51

# Boxplot of Ip as a function of aridity levels
df_summarized_complete = as.data.frame(df_summarized_complete)

df_summarized_complete <- df_summarized_complete %>%
  mutate(aridity_group = case_when(
    Aridity >= 0    & Aridity < 0.2  ~ "Arid",
    Aridity >= 0.2  & Aridity < 0.65 ~ "Semiarid",
    Aridity >= 0.65 & Aridity < 6    ~ "Humid"
  )) %>%
  mutate(aridity_group = factor(aridity_group, levels = c("Arid", "Semiarid", "Humid")))
head(df_summarized_complete)

COLS_ARID <- c(
  "Arid"     = "#e31a1c",
  "Semiarid" = "#f89540",
  "Humid"    = "#666666"
)

top10_phyla <- df_summarized_complete %>%
  group_by(Phylum) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

# Subset dataframe to keep only top 10 phyla
df_top10_phyla <- df_summarized_complete %>%
  filter(Phylum %in% top10_phyla) %>%
  mutate(Phylum = factor(Phylum, levels = top10_phyla)) 
head(df_top10_phyla)

df_top10_phyla <- df_top10_phyla %>%
  mutate(basic_acid_ratio = basic.aa.percent / acid.aa.percent)

signif_df <- df_top10_phyla %>%
  group_by(Phylum) %>%
  summarise(
    p.value = tryCatch({
      wilcox.test(
        acid.aa.percent[aridity_group == "Arid"],
        acid.aa.percent[aridity_group == "Humid"]
      )$p.value
    }, error = function(e) NA_real_)
  ) %>%
  mutate(
    p.signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value >= 0.05 ~ "ns",
      TRUE            ~ NA_character_
    ),
    # Place label just above the top of each boxplot group
    y_pos = max(df_top10_phyla$acid.aa.percent, na.rm = TRUE) * 1.05
  )

nrow(df_top10_phyla) #1918 MAGs belong to the top 10 phyla

act.chl.bacil.df_top10_phyla = subset(df_top10_phyla, Phylum %in% c("Actinobacteriota", "Chloroflexota", "Bacillota"))
head(act.chl.bacil.df_top10_phyla)

acid.bact.pseud.df_top10_phyla = subset(df_top10_phyla, Phylum %in% c("Acidobacteriota", "Bacteroidota", "Pseudomonadota"))
head(acid.bact.pseud.df_top10_phyla)

other.df_top10_phyla = subset(df_top10_phyla, !Phylum %in% c("Actinobacteriota", "Chloroflexota", "Bacillota"))
head(other.df_top10_phyla)

median(act.chl.bacil.df_top10_phyla$Ip)
IQR(act.chl.bacil.df_top10_phyla$Ip)

median(acid.bact.pseud.df_top10_phyla$Ip)
IQR(acid.bact.pseud.df_top10_phyla$Ip)

median(other.df_top10_phyla$Ip)
IQR(other.df_top10_phyla$Ip)

wilcox.test(act.chl.bacil.df_top10_phyla$Ip, other.df_top10_phyla$Ip)

sum(df_top10_phyla$Phylum == "Patescibacteriota", na.rm = TRUE)
sum(df_top10_phyla$Phylum == "Verrucomicrobiota", na.rm = TRUE)

pates.df_top10_phyla = subset(df_top10_phyla, Phylum == "Patescibacteriota")
verr.df_top10_phyla = subset(df_top10_phyla, Phylum == "Verrucomicrobiota")

arid.pates.df_top10_phyla = subset(pates.df_top10_phyla, pates.df_top10_phyla$aridity_group %in%
                                     c("Arid", "Semiarid"))
arid.verr.df_top10_phyla = subset(verr.df_top10_phyla, verr.df_top10_phyla$aridity_group %in%
                                     c("Arid", "Semiarid"))
humid.pates.df_top10_phyla = subset(pates.df_top10_phyla, pates.df_top10_phyla$aridity_group %in%
                                     c("Humid"))
humid.verr.df_top10_phyla = subset(verr.df_top10_phyla, verr.df_top10_phyla$aridity_group %in%
                                    c("Humid"))

wilcox.test(arid.pates.df_top10_phyla$Ip, humid.pates.df_top10_phyla$Ip)
wilcox.test(arid.verr.df_top10_phyla$Ip, humid.verr.df_top10_phyla$Ip)

kruskal.test(pates.df_top10_phyla$Ip ~ pates.df_top10_phyla$aridity_group)
kruskal.test(verr.df_top10_phyla$Ip ~ verr.df_top10_phyla$aridity_group)

acid.percent.boxplot.phyla <- ggplot(df_top10_phyla, aes(x = Phylum, y = acid.aa.percent, fill = aridity_group)) +
  geom_boxplot(width = 0.7, 
               alpha = 0.6,
               color = "grey20", 
               outlier.shape = NA,
               position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = aridity_group), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 0.6,
              alpha = 0.4, 
              show.legend = FALSE) +
  scale_fill_manual(values = COLS_ARID, name = "Aridity group") +
  scale_color_manual(values = COLS_ARID, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  labs(
    x = "Phylum",
    y = "Proportion of acidic amino acids (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 5, r = 5, b = 20, l = 30, unit = "pt")
  ) +
  geom_text(
    data = signif_df,
    aes(x = Phylum, y = y_pos, label = p.signif),
    inherit.aes = FALSE,
    size = 4,
    vjust = 0
  )
acid.percent.boxplot.phyla

tiff(file.path(OUT_DIR, "Fig4B_acidpercent_boxplot_mags_phyla.tiff"), 
     units="cm", width=12, height=10, res=300)
acid.percent.boxplot.phyla
dev.off()


basic.percent.boxplot.phyla <- ggplot(df_top10_phyla, aes(x = Phylum, y = basic.aa.percent, fill = aridity_group)) +
  geom_boxplot(width = 0.7, 
               alpha = 0.6,
               color = "grey20", 
               outlier.shape = NA,
               position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = aridity_group), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 0.6,
              alpha = 0.4, 
              show.legend = FALSE) +
  scale_fill_manual(values = COLS_ARID, name = "Aridity group") +
  scale_color_manual(values = COLS_ARID, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    x = "Phylum",
    y = "Proportion of basic amino acids (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 5, r = 5, b = 20, l = 30, unit = "pt")
  ) +
  stat_compare_means(
    aes(group = aridity_group),
    method = "kruskal.test",       # one p-value per phylum
    label = "p.signif",
    label.y.npc = 0.95,
    size = 3.2
  )
basic.percent.boxplot.phyla

tiff(file.path(OUT_DIR, "Fig4C_basicpercent_boxplot_mags_phyla.tiff"), 
     units="cm", width=12, height=10, res=300)
basic.percent.boxplot.phyla
dev.off()


acid.basic.ratio.boxplot.phyla <- ggplot(df_top10_phyla, aes(x = Phylum, y = basic.aa.percent/acid.aa.percent, fill = aridity_group)) +
  geom_boxplot(width = 0.7, 
               alpha = 0.6,
               color = "grey20", 
               outlier.shape = NA,
               position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = aridity_group), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 0.6,
              alpha = 0.4, 
              show.legend = FALSE) +
  scale_fill_manual(values = COLS_ARID, name = "Aridity group") +
  scale_color_manual(values = COLS_ARID, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    x = "Phylum",
    y = "Basic / Acidic amino acids"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 5, r = 5, b = 20, l = 30, unit = "pt")
  ) +
  stat_compare_means(
    aes(group = aridity_group),
    method = "kruskal.test",       # one p-value per phylum
    label = "p.signif",
    label.y.npc = 0.95,
    size = 3.2
  )
acid.basic.ratio.boxplot.phyla

tiff(file.path(OUT_DIR, "Fig4D_basicacidratio_boxplot_mags_phyla.tiff"), 
     units="cm", width=12, height=10, res=300)
acid.basic.ratio.boxplot.phyla
dev.off()


# pI

signif_df <- df_top10_phyla %>%
  group_by(Phylum) %>%
  summarise(
    p.value = tryCatch({
      wilcox.test(
        acid.aa.percent[aridity_group == "Arid"],
        acid.aa.percent[aridity_group == "Humid"]
      )$p.value
    }, error = function(e) NA_real_)
  ) %>%
  mutate(
    p.signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value >= 0.05 ~ "ns",
      TRUE            ~ NA_character_
    ),
    # Place label just above the top of each boxplot group
    y_pos = max(df_top10_phyla$Ip, na.rm = TRUE) * 1.05
  )

pi.boxplot.phyla <- ggplot(df_top10_phyla, aes(x = Phylum, y = Ip, fill = aridity_group)) +
  geom_boxplot(width = 0.7, 
               alpha = 0.6,
               color = "grey20", 
               outlier.shape = NA,
               position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = aridity_group), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 0.6,
              alpha = 0.4, 
              show.legend = FALSE) +
  scale_fill_manual(values = COLS_ARID, name = "Aridity group") +
  scale_color_manual(values = COLS_ARID, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  labs(
    x = "Phylum",
    y = "Isoelectric point (pI)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 5, r = 5, b = 20, l = 30, unit = "pt")
  ) +
  geom_text(
    data = signif_df,
    aes(x = Phylum, y = y_pos, label = p.signif),
    inherit.aes = FALSE,
    size = 4,
    vjust = 0
  )
pi.boxplot.phyla

tiff(file.path(OUT_DIR, "Fig4E_pI_boxplot_mags_phyla.tiff"), 
     units="cm", width=12, height=10, res=300)
pi.boxplot.phyla
dev.off()


### Phylogenetic analysis
library(ape)
library(phytools)

mag.tree = read.tree(MAG_TREE_FILE)
mag.pi = read.csv(MAG_PI_TIPS_FILE,
                  header = T)
head(mag.pi)

pI <- setNames(mag.pi$pI, mag.pi$Genome)
common <- intersect(mag.tree$tip.label, names(pI))
tree <- drop.tip(mag.tree, setdiff(mag.tree$tip.label, common))
trait <- pI[common]

# Blomberg's K
K <- phylosig(tree,
              trait,
              method = "K",
              test = TRUE)
K

# Phylogenetic signal K : 0.945807 
# P-value (based on 1000 randomizations) : 0.001

# Pagel's lambda
lambda <- phylosig(tree,
                   trait,
                   method = "lambda",
                   test = TRUE)
lambda

# Phylogenetic signal lambda : 0.986799 
# logL(lambda) : 575.391 
# LR(lambda=0) : 3840.75 
# P-value (based on LR test) : 0.001 

# Plot
plot(contMap(tree, trait))
