rm(list = ls())

############################################################
# CONFIGURATION: edit these paths before running
############################################################
# Run from the repository root so that "data/" resolves.
DATA_DIR <- "data"                   # files included in this repository
EXT_DIR  <- "/YOURDIRECTORY"         # inputs not included here (see README)
OUT_DIR  <- "figures"                # figures are written here
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
ARTHRO_DIR <- file.path(EXT_DIR, "Arthrobacter_test")

library(data.table)
library(arrow)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(stringr)
library(tidyr)
library(dplyr)

##########################################################################################
# Explorative analysis of 53 Arthrobacter isolates tested for desiccation tolerance and proteome properties 
arthrobacter.id = read.csv(file.path(ARTHRO_DIR, "Arthrobacter_IDs.csv"), header = T)
arthrobacter.pepstats = read.csv(file.path(ARTHRO_DIR, "Arthrobacter_proteomes_summary.csv"), header = T)

arthrobacter.dbb = merge(arthrobacter.id, arthrobacter.pepstats, by = "Genome")
arthrobacter.dbb$Group = gsub(" strains", "", arthrobacter.dbb$Group)
head(arthrobacter.dbb)

# Acidic AAs
wt <- wilcox.test(mole_pct_acidic ~ Group, data = arthrobacter.dbb)
p_label <- paste0("Mann-Whitney U, P = ", formatC(wt$p.value, format = "e", digits = 1))

acidic.origin <- ggplot(arthrobacter.dbb, aes(x = Group, y = mole_pct_acidic, fill = Group)) + 
  theme_classic() +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_jitter(aes(color = Group), size = 2, alpha = 0.5, width = 0.2) +
  annotate("text", x = 0.55, y = min(arthrobacter.dbb$mole_pct_acidic) * 1.125,
           label = p_label, size = 4, hjust = 0) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12), 
        legend.position = "null",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)), 
        axis.title.y = element_text(size = 14, margin = margin(r = 10)), 
        plot.title = element_text(size = 14, face = "bold", colour = "Black"), 
        plot.margin = margin(1, 1, 1, 1, "cm")) + 
  ylab("Proportion of acidic amino acids (%)") + 
  xlab("Genome origin") + 
  scale_fill_manual(values = c("purple", "orange")) + 
  scale_color_manual(values = c("purple", "orange"))
acidic.origin

tiff(file.path(OUT_DIR, "Acidic_AntarcticGroup_Boxplot_2026.tiff"), units="cm", width=10, height=14, res = 300, compression = 'lzw')
acidic.origin
dev.off()



# Charged AA
wt <- wilcox.test(mole_pct_basic ~ Group, data = arthrobacter.dbb)
p_label <- paste0("Mann-Whitney U, P = ", formatC(wt$p.value, format = "e", digits = 1))

basic.origin <- ggplot(arthrobacter.dbb, aes(x = Group, y = mole_pct_basic, fill = Group)) + 
  theme_classic() +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_jitter(aes(color = Group), size = 2, alpha = 0.5, width = 0.2) +
  annotate("text", x = 0.55, y = min(arthrobacter.dbb$mole_pct_basic) * 1.1,
           label = p_label, size = 4, hjust = 0) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12), 
        legend.position = "null",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)), 
        axis.title.y = element_text(size = 14, margin = margin(r = 10)), 
        plot.title = element_text(size = 14, face = "bold", colour = "Black"), 
        plot.margin = margin(1, 1, 1, 1, "cm")) + 
  ylab("Proportion of basic amino acids (%)") + 
  xlab("Genome origin") + 
  scale_fill_manual(values = c("purple", "orange")) + 
  scale_color_manual(values = c("purple", "orange"))
basic.origin

tiff(file.path(OUT_DIR, "Basic_AntarcticGroup_Boxplot_2026.tiff"), units="cm", width=10, height=14, res = 300, compression = 'lzw')
basic.origin
dev.off()


# Charged AA
wt <- wilcox.test(mole_pct_charged ~ Group, data = arthrobacter.dbb)
p_label <- paste0("Mann-Whitney U, P = ", formatC(wt$p.value, format = "e", digits = 1))

charged.origin <- ggplot(arthrobacter.dbb, aes(x = Group, y = mole_pct_charged, fill = Group)) + 
  theme_classic() +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_jitter(aes(color = Group), size = 2, alpha = 0.5, width = 0.2) +
  annotate("text", x = 0.55, y = min(arthrobacter.dbb$mole_pct_charged) * 1.1,
           label = p_label, size = 4, hjust = 0) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12), 
        legend.position = "null",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)), 
        axis.title.y = element_text(size = 14, margin = margin(r = 10)), 
        plot.title = element_text(size = 14, face = "bold", colour = "Black"), 
        plot.margin = margin(1, 1, 1, 1, "cm")) + 
  ylab("Proportion of charged amino acids (%)") + 
  xlab("Genome origin") + 
  scale_fill_manual(values = c("purple", "orange")) + 
  scale_color_manual(values = c("purple", "orange"))
charged.origin

tiff(file.path(OUT_DIR, "Charged_AntarcticGroup_Boxplot_2026.tiff"), units="cm", width=10, height=14, res = 300, compression = 'lzw')
charged.origin
dev.off()


# Isoelectric point (pI)
wt <- wilcox.test(mean_pI ~ Group, data = arthrobacter.dbb)
p_label <- paste0("Mann-Whitney U, P = ", formatC(wt$p.value, digits = 3))

pi.origin <- ggplot(arthrobacter.dbb, aes(x = Group, y = mean_pI, fill = Group)) + 
  theme_classic() +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_jitter(aes(color = Group), size = 2, alpha = 0.5, width = 0.2) +
  annotate("text", x = 0.55, y = min(arthrobacter.dbb$mean_pI) * 1.11,
           label = p_label, size = 4, hjust = 0) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12), 
        legend.position = "null",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)), 
        axis.title.y = element_text(size = 14, margin = margin(r = 10)), 
        plot.title = element_text(size = 14, face = "bold", colour = "Black"), 
        plot.margin = margin(1, 1, 1, 1, "cm")) + 
  ylab("Isoelectric point (pI)") + 
  xlab("Genome origin") + 
  scale_fill_manual(values = c("purple", "orange")) + 
  scale_color_manual(values = c("purple", "orange"))
pi.origin

tiff(file.path(OUT_DIR, "Ip_AntarcticGroup_Boxplot_2026.tiff"), units="cm", width=10, height=14, res = 300, compression = 'lzw')
pi.origin
dev.off()


# Test association between proteome features and desiccation tolerance
desic.dbb = read.csv(file.path(ARTHRO_DIR, "Arthrobacter_desiccation_data.csv"), header = T)
head(desic.dbb)

arthrobacter.dbb.desic = merge(arthrobacter.dbb, desic.dbb, by = "Genome")
head(arthrobacter.dbb.desic)

# Spearman correlations
cor_test <- cor.test(
  arthrobacter.dbb.desic$Predicted.50..Growth.Reduction,
  arthrobacter.dbb.desic$mole_pct_basic,
  method = "spearman"
)
cor_test

cor_test.acid <- cor.test(
  arthrobacter.dbb.desic$Predicted.50..Growth.Reduction,
  arthrobacter.dbb.desic$mole_pct_acidic,
  method = "spearman"
)
cor_test.acid

cor_test.pI <- cor.test(
  arthrobacter.dbb.desic$Predicted.50..Growth.Reduction,
  arthrobacter.dbb.desic$mean_pI,
  method = "spearman"
)
cor_test.pI


label_cor <- paste0(
  "Spearman rho = ", round(cor_test$estimate, 2),
  ", P = ", formatC(cor_test$p.value, digits = 2)
)

medians_df <- arthrobacter.dbb.desic %>%
  group_by(Predicted.50..Growth.Reduction) %>%
  summarise(median_y = median(mole_pct_basic, na.rm = TRUE))

# Plot

COLS_DESIC <- c(
  "-2.5" = "#7F0000",
  "-2"   = "#B2182B",
  "-1.5" = "#D6604D",
  "-1"   = "#F4A582",
  "-0.5" = "#D9D9D9" 
)

arthrobacter.dbb.desic$Predicted.50..Growth.Reduction = as.factor(arthrobacter.dbb.desic$Predicted.50..Growth.Reduction)
arthrobacter.dbb.desic$Predicted.50..Growth.Reduction = factor(arthrobacter.dbb.desic$Predicted.50..Growth.Reduction,
                                                               levels = c("-2.5", "-2", "-1.5", "-1", "-0.5"))


p <-
  ggplot(
    arthrobacter.dbb.desic,
    aes(
      x = factor(Predicted.50..Growth.Reduction, 
                 levels = c("-2.5", "-2", "-1.5", "-1", "-0.5")),
      y = mole_pct_basic,
      color = factor(Predicted.50..Growth.Reduction)
    )
  ) +
  geom_jitter(
    size = 3,
    alpha = 0.75,
    width = 0.15
  ) +
  geom_segment(
    data = medians_df,
    aes(
      x = as.numeric(factor(Predicted.50..Growth.Reduction,
                            levels = c("-2.5", "-2", "-1.5", "-1", "-0.5"))) - 0.2,
      xend = as.numeric(factor(Predicted.50..Growth.Reduction,
                               levels = c("-2.5", "-2", "-1.5", "-1", "-0.5"))) + 0.2,
      y = median_y,
      yend = median_y
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = 0.75,
    y = Inf,
    label = label_cor,
    hjust = 0,
    vjust = 1.5,
    size = 4
  ) +
  scale_color_manual(
    values = COLS_DESIC,
    limits = c("-2.5", "-2", "-1.5", "-1", "-0.5")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    x = "Predicted 50% growth reduction",
    y = "Proportion of basic amino acids (%)"
  ) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"
  )
p

tiff(file.path(OUT_DIR, "Arthrobacter_BasicAA_desiccation_spearman_2026.tiff"), units="cm", width=10, height=10, res = 300, compression = 'lzw')
p
dev.off()


################################################################################
###################FIGURE 1D: Density plot for protein pI#######################
################################################################################
## Read in the per protein pepstats results of the Arthrobacter isolates 

file_paths <- list.files(path = file.path(ARTHRO_DIR, "pepstats_arthrobacter"), 
                         pattern = "\\.csv$", 
                         full.names = TRUE)
file_paths

# Read all files and add filename column
combined_data <- rbindlist(lapply(file_paths, function(file) {
  dt <- fread(file)
  dt[, source_file := basename(file)]
}), fill = TRUE)
head(combined_data)

combined_data$Genome = gsub("_per_protein.csv", "", combined_data$source_file)
head(combined_data)

# Read in the desiccation tolerance data

desic.dbb = read.csv(file.path(ARTHRO_DIR, "Arthrobacter_desiccation_data.csv"), header = T)
head(desic.dbb)

arthrobacter.dbb.desic = merge(combined_data, desic.dbb, by = "Genome")
head(arthrobacter.dbb.desic)

#arthrobacter.dbb.desic = subset(arthrobacter.dbb.desic, Predicted.50..Growth.Reduction > -2.5) # Optional removal of single -2.5 value
#head(arthrobacter.dbb.desic)

arthrobacter.dbb.desic$Predicted.50..Growth.Reduction = as.factor(arthrobacter.dbb.desic$Predicted.50..Growth.Reduction)
arthrobacter.dbb.desic$Predicted.50..Growth.Reduction = factor(arthrobacter.dbb.desic$Predicted.50..Growth.Reduction,
                                                               levels = c("-2.5", "-2", "-1.5", "-1", "-0.5"))
#arthrobacter.dbb.desic$Predicted.50..Growth.Reduction = factor(arthrobacter.dbb.desic$Predicted.50..Growth.Reduction,
#                                                               levels = c("-2", "-1.5", "-1", "-0.5")) # removal -2.5
head(arthrobacter.dbb.desic)


COLS_DESIC <- c(
  "-2.5" = "#7F0000",
  "-2"   = "#B2182B",
  "-1.5" = "#D6604D",
  "-1"   = "#F4A582",
  "-0.5" = "#D9D9D9" 
)

#COLS_DESIC <- c(
#  "-2"   = "#B2182B",
#  "-1.5" = "#D6604D",
#  "-1"   = "#F4A582",
#  "-0.5" = "#D9D9D9" 
#)


dp = ggplot(arthrobacter.dbb.desic, aes(x = pI, colour = Predicted.50..Growth.Reduction)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  scale_colour_manual(values = COLS_DESIC, name = "Predicted 50% growth reduction") +
  labs(
    x = "Protein isoelectric point (pI)",
    y = "Density"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 14))
dp

tiff(file.path(OUT_DIR, "Arthrobacter_desiccation_density_plot_no-2p5.tiff"), 
     units="cm", width=18, height=10, res=300)
dp
dev.off()


################################################################################
#######################SIGNALP ARTHROBACTER ANALYSIS############################
################################################################################

sp.arth = read.csv(file.path(ARTHRO_DIR, "SignalP_arthrobacter/signalp_arthrobacter_pepstats_database.csv"),
                   header = T)
head(sp.arth)

sp.arth.desic = merge(sp.arth, desic.dbb, by = "Genome")
unique(sp.arth.desic$localization) # "cytoplasmic", "secreted", "transmembrane"
head(sp.arth.desic)

# Cytoplasmic (pI, percentage acidic, percentage basic)
cyto.sp.arth.desic = subset(sp.arth.desic, localization == "cytoplasmic")
head(cyto.sp.arth.desic)

cyto.dataset <- 
  cyto.sp.arth.desic %>%
  group_by(Genome) %>%
  summarise(
    acid.aa.percent  = mean(acidic_pct),
    basic.aa.percent = mean(basic_pct),
    mean.pI = mean(pI),
    Predicted.50..Growth.Reduction = first(Predicted.50..Growth.Reduction),
    .groups = "drop"
  )
head(cyto.dataset)

# Spearman correlation: basic rho = 0.212, P = 0.128; acid rho = -0.034, P = 0.806; pI rho = 0.119, P = 0.397;  
cor_test <- cor.test(
  cyto.dataset$Predicted.50..Growth.Reduction,
  cyto.dataset$acid.aa.percent,
  method = "spearman"
)

## No significant associations with desiccation tolerance in cytoplasmic proteins

# Transmembrane (pI, percentage acidic, percentage basic)
trans.sp.arth.desic = subset(sp.arth.desic, localization == "transmembrane")
head(trans.sp.arth.desic)

trans.dataset <- 
  trans.sp.arth.desic %>%
  group_by(Genome) %>%
  summarise(
    acid.aa.percent  = mean(acidic_pct),
    basic.aa.percent = mean(basic_pct),
    mean.pI = mean(pI),
    Predicted.50..Growth.Reduction = first(Predicted.50..Growth.Reduction),
    .groups = "drop"
  )
head(trans.dataset)

# Spearman correlation: basic rho = 0.268, P = 0.053; acid rho = 0.104, P = 0.458; pI rho = -0.031, P = 0.823;  
cor_test <- cor.test(
  trans.dataset$Predicted.50..Growth.Reduction,
  trans.dataset$basic.aa.percent,
  method = "spearman"
)

## Marginal association of basic AA percentage with desiccation tolerance in transmembrane proteins


# Secreted (pI, percentage acidic, percentage basic)
secret.sp.arth.desic = subset(sp.arth.desic, localization == "secreted")
head(secret.sp.arth.desic)

secret.dataset <- 
  secret.sp.arth.desic %>%
  group_by(Genome) %>%
  summarise(
    acid.aa.percent  = mean(acidic_pct),
    basic.aa.percent = mean(basic_pct),
    mean.pI = mean(pI),
    Predicted.50..Growth.Reduction = first(Predicted.50..Growth.Reduction),
    .groups = "drop"
  )
head(secret.dataset)

# Spearman correlation: basic rho = 0.193, P = 0.166; acid rho = 0.011, P = 0.939; pI rho = 0.192, P = 0.167;  
cor_test <- cor.test(
  secret.dataset$Predicted.50..Growth.Reduction,
  secret.dataset$mean.pI,
  method = "spearman"
)

## No significant associations with desiccation tolerance in secreted proteins


## Acid/Basic ratio calculation

file_paths <- list.files(path = file.path(ARTHRO_DIR, "pepstats_arthrobacter"), 
                         pattern = "\\.csv$", 
                         full.names = TRUE)
file_paths

# Read all files and add filename column
combined_data <- rbindlist(lapply(file_paths, function(file) {
  dt <- fread(file)
  dt[, source_file := basename(file)]
}), fill = TRUE)
head(combined_data)

combined_data$Genome = gsub("_per_protein.csv", "", combined_data$source_file)
head(combined_data)

# Read in the desiccation tolerance data

desic.dbb = read.csv(file.path(ARTHRO_DIR, "Arthrobacter_desiccation_data.csv"), header = T)
head(desic.dbb)

arthrobacter.dbb.desic = merge(combined_data, desic.dbb, by = "Genome")
head(arthrobacter.dbb.desic)

arthrobacter.dbb.desic$Predicted.50..Growth.Reduction = as.factor(arthrobacter.dbb.desic$Predicted.50..Growth.Reduction)
arthrobacter.dbb.desic$Predicted.50..Growth.Reduction = factor(arthrobacter.dbb.desic$Predicted.50..Growth.Reduction,
                                                               levels = c("-2.5", "-2", "-1.5", "-1", "-0.5"))
head(arthrobacter.dbb.desic)

# Filter to keep only those proteins with a Ip between 4.5-5.0 or 9.0-10.0 (the most distinctive peaks)
filtered_data <- arthrobacter.dbb.desic[
  (arthrobacter.dbb.desic$pI >= 4.5 & arthrobacter.dbb.desic$pI <= 5) | 
    (arthrobacter.dbb.desic$pI >= 9 & arthrobacter.dbb.desic$pI <= 10),]
head(filtered_data)

# Count number of proteins within each metaproteome which have Ip 4.5-5.0 or 9.0-10.0
range_counts <- filtered_data %>%
  group_by(Genome) %>%
  summarise(
    Acid.prots = sum(pI >= 4.5 & pI <= 5),
    Basic.prots = sum(pI >= 9 & pI <= 10),
    total_count = n()  # Optional: total entries per file
  ) %>%
  ungroup()
range_counts$Basic.acid.ratio = range_counts$Basic.prots / range_counts$Acid.prots
head(range_counts)

acid.basic.dbb.desic = merge(range_counts, desic.dbb, by = "Genome")
head(acid.basic.dbb.desic)

# Spearman correlation
cor_test <- cor.test(
  acid.basic.dbb.desic$Predicted.50..Growth.Reduction,
  acid.basic.dbb.desic$Basic.acid.ratio,
  method = "spearman"
)

## No significant association of the basic/acid ratio with desiccation tolerance in Arthrobacter


## Basic / Acid ratio

ratio.dataset = read.csv(file.path(ARTHRO_DIR, "Basic_acid_ratio_dataset_arthrobacter.csv"), header = T)
head(ratio.dataset)

#ratio.dataset = subset(ratio.dataset, Predicted.50..Growth.Reduction > -2.5) # optional removal -2.5
#head(ratio.dataset) # optional removal -2.5

# Spearman correlation
cor_test <- cor.test(
  ratio.dataset$Predicted.50..Growth.Reduction,
  ratio.dataset$Acid.prots,
  method = "spearman"
) # rho = -0.248, P = 0.074
cor_test

label_cor <- paste0(
  "Spearman rho = ", round(cor_test$estimate, 2),
  ", P = ", formatC(cor_test$p.value, digits = 2)
)

medians_df <- ratio.dataset %>%
  group_by(Predicted.50..Growth.Reduction) %>%
  summarise(median_y = median(Acid.prots, na.rm = TRUE))

# Plot

COLS_DESIC <- c(
  "-2.5" = "#7F0000",
  "-2"   = "#B2182B",
  "-1.5" = "#D6604D",
  "-1"   = "#F4A582",
  "-0.5" = "#D9D9D9" 
)

ratio.dataset$Predicted.50..Growth.Reduction = as.factor(ratio.dataset$Predicted.50..Growth.Reduction)
ratio.dataset$Predicted.50..Growth.Reduction = factor(ratio.dataset$Predicted.50..Growth.Reduction,
                                                               levels = c("-2.5", "-2", "-1.5", "-1", "-0.5"))

p <-
  ggplot(
    ratio.dataset,
    aes(
      x = factor(Predicted.50..Growth.Reduction, 
                 levels = c("-2.5", "-2", "-1.5", "-1", "-0.5")),
      y = Acid.prots,
      color = factor(Predicted.50..Growth.Reduction)
    )
  ) +
  geom_jitter(
    size = 3,
    alpha = 0.75,
    width = 0.15
  ) +
  geom_segment(
    data = medians_df,
    aes(
      x = as.numeric(factor(Predicted.50..Growth.Reduction,
                            levels = c("-2.5", "-2", "-1.5", "-1", "-0.5"))) - 0.2,
      xend = as.numeric(factor(Predicted.50..Growth.Reduction,
                               levels = c("-2.5", "-2", "-1.5", "-1", "-0.5"))) + 0.2,
      y = median_y,
      yend = median_y
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1
  ) +
  scale_y_log10(  # Log10 scale with custom labels
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = function(x) round(log10(x), 1)  # Show only the log10 values, rounded to 1 decimal
  ) +
  annotate(
    "text",
    x = 0.75,
    y = Inf,
    label = label_cor,
    hjust = 0,
    vjust = 1.5,
    size = 4
  ) +
  scale_color_manual(
    values = COLS_DESIC,
    limits = c("-2.5", "-2", "-1.5", "-1", "-0.5")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    x = "Predicted 50% growth reduction",
    y = "Number of acidic proteins"
  ) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"
  )
p

tiff(file.path(OUT_DIR, "Arthrobacter_TotalAcidProts_desiccation_spearman_2026.tiff"), units="cm", width=10, height=10, res = 300, compression = 'lzw')
p
dev.off()



# Spearman correlation
cor_test <- cor.test(
  ratio.dataset$Predicted.50..Growth.Reduction,
  ratio.dataset$Basic.prots,
  method = "spearman"
) # rho = -0.066, P = 0.637

label_cor <- paste0(
  "Spearman rho = ", round(cor_test$estimate, 2),
  ", P = ", formatC(cor_test$p.value, digits = 2)
)

medians_df <- ratio.dataset %>%
  group_by(Predicted.50..Growth.Reduction) %>%
  summarise(median_y = median(Basic.prots, na.rm = TRUE))

# Plot

COLS_DESIC <- c(
  "-2.5" = "#7F0000",
  "-2"   = "#B2182B",
  "-1.5" = "#D6604D",
  "-1"   = "#F4A582",
  "-0.5" = "#D9D9D9" 
)

ratio.dataset$Predicted.50..Growth.Reduction = as.factor(ratio.dataset$Predicted.50..Growth.Reduction)
ratio.dataset$Predicted.50..Growth.Reduction = factor(ratio.dataset$Predicted.50..Growth.Reduction,
                                                      levels = c("-2.5", "-2", "-1.5", "-1", "-0.5"))

p <-
  ggplot(
    ratio.dataset,
    aes(
      x = factor(Predicted.50..Growth.Reduction, 
                 levels = c("-2.5", "-2", "-1.5", "-1", "-0.5")),
      y = Acid.prots,
      color = factor(Predicted.50..Growth.Reduction)
    )
  ) +
  geom_jitter(
    size = 3,
    alpha = 0.75,
    width = 0.15
  ) +
  geom_segment(
    data = medians_df,
    aes(
      x = as.numeric(factor(Predicted.50..Growth.Reduction,
                            levels = c("-2.5", "-2", "-1.5", "-1", "-0.5"))) - 0.2,
      xend = as.numeric(factor(Predicted.50..Growth.Reduction,
                               levels = c("-2.5", "-2", "-1.5", "-1", "-0.5"))) + 0.2,
      y = median_y,
      yend = median_y
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1
  ) +
  scale_y_log10(  # Log10 scale with custom labels
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = function(x) round(log10(x), 1)  # Show only the log10 values, rounded to 1 decimal
  ) +
  annotate(
    "text",
    x = 0.75,
    y = Inf,
    label = label_cor,
    hjust = 0,
    vjust = 1.5,
    size = 4
  ) +
  scale_color_manual(
    values = COLS_DESIC,
    limits = c("-2.5", "-2", "-1.5", "-1", "-0.5")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    x = "Predicted 50% growth reduction",
    y = "Number of basic proteins"
  ) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"
  )
p

tiff(file.path(OUT_DIR, "Arthrobacter_TotalBasicProts_desiccation_spearman_2026.tiff"), units="cm", width=10, height=10, res = 300, compression = 'lzw')
p
dev.off()

