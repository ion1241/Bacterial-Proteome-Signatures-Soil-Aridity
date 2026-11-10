rm(list = ls())

############################################################
# CONFIGURATION: edit these paths before running
############################################################
# Run from the repository root so that "data/" resolves.
DATA_DIR <- "data"                   # files included in this repository
EXT_DIR  <- "/YOURDIRECTORY"         # inputs not included here (see README)
OUT_DIR  <- "figures"                # figures are written here
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
AI_FILE <- file.path(EXT_DIR, "Global-AI_ET0__annual_v3_1", "ai_v31_yr.tif")
MAG_METADATA_FILE <- file.path(EXT_DIR, "metadata_mags.tsv")
MAG_SUMMARY_FILE <- file.path(DATA_DIR, "MAGs_aridity_proteome_summary.csv")
MAG_FEATURES_PARQUET_DIR <- file.path(EXT_DIR, "Data_mags_parquet", "proteome_ai6500")

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

############################################################
# FILE PATHS
############################################################
# Global MAGs
POINTS_GLOBAL_FILE <- MAG_METADATA_FILE

############################################################
# COLOR PALETTE (CONSISTENT ACROSS MAPS)
############################################################
aridity_colors <- c(
  "Hyper-arid"    = "#e31a1c",
  "Arid"          = "#fd8d3c",
  "Semi-arid"     = "#fecc5c",
  "Dry sub-humid" = "#ffffb2",
  "Humid"         = "#666666"
)

############################################################
# READ ARIDITY INDEX RASTER
############################################################

ai_rast <- terra::rast(AI_FILE)

# Global land polygons (for masking oceans)
world_sf <- rnaturalearth::ne_countries(
  scale = "medium", returnclass = "sf"
)
world_vect <- terra::vect(world_sf)

# Aggregate + mask oceans
ai_global <- terra::aggregate(
  ai_rast, fact = 6, fun = "mean", na.rm = TRUE
)
ai_global <- terra::mask(ai_global, world_vect)

# Raster to dataframe (oceans removed)
ai_global_df <- as.data.frame(ai_global, xy = TRUE, na.rm = TRUE) |>
  rename(ai = 3) |>
  mutate(
    ai = ai / 10000,
    category = case_when(
      is.infinite(ai) | ai >= 0.65 ~ "Humid",
      ai >= 0.50 & ai < 0.65       ~ "Dry sub-humid",
      ai >= 0.20 & ai < 0.50       ~ "Semi-arid",
      ai >= 0.05 & ai < 0.20       ~ "Arid",
      ai <  0.05                   ~ "Hyper-arid",
      TRUE                         ~ NA_character_
    ),
    category = factor(
      category,
      levels = c("Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid", "Humid")
    )
  ) |>
  filter(!is.na(category))

# Parse MAG coordinates like "50.958 N 6.856 E"
points_global_df <- read_tsv(POINTS_GLOBAL_FILE, show_col_types = FALSE) |>
  mutate(
    lat = as.numeric(str_extract(lat_lon, "^[0-9.]+")) *
      ifelse(str_detect(lat_lon, " S"), -1, 1),
    lon = as.numeric(str_extract(lat_lon, "(?<= )[0-9.]+")) *
      ifelse(str_detect(lat_lon, " W"), -1, 1)
  ) |>
  select(-lat_lon)

# Global map
p_map_global <- ggplot() +
  geom_tile(data = ai_global_df, aes(x = x, y = y, fill = category)) +
  scale_fill_manual(values = aridity_colors, name = "Aridity class") +
  geom_sf(data = world_sf, fill = NA, colour = "grey30",
          linewidth = 0.3, inherit.aes = FALSE) +
  geom_jitter(
    data = points_global_df,
    aes(x = lon, y = lat),
    shape = 21, fill = "white", colour = "black",
    size = 1.5, width = 0.7, height = 0.7
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold")
  )

p_map_global

ggsave(
  filename = file.path(OUT_DIR, "Fig2A_map_mags.tiff"),
  plot     = p_map_global,
  width    = 14, height = 16,
  units    = "cm", dpi = 300,
  compression = "lzw"
)


################################################################################
##############FIGURE 2B: Facetted boxplot for acidic and basic AAs##############
################################################################################

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
  mutate(aridity_groups = factor(aridity_group, levels = c("Arid", "Semiarid", "Humid")))
head(df_summarized_complete)

COLS_ARID <- c(
  "Arid"     = "#e31a1c",
  "Semiarid" = "#f89540",
  "Humid"    = "#666666"
)

plot_df_long <- df_summarized_complete %>%
    pivot_longer(
    cols = c(acid.aa.percent, basic.aa.percent),
    names_to = "aa_type",
    values_to = "aa_percent"
  ) %>%
  mutate(
    aa_type = recode(
      aa_type,
      "acid.aa.percent"  = "Acidic amino acids",
      "basic.aa.percent" = "Basic amino acids"
    )
  )

plot_df_long$aridity_group = factor(plot_df_long$aridity_group, levels = c("Arid", "Semiarid", "Humid"))

# Test for normality
arid.plot_df_long = subset(plot_df_long, aridity_group == "Arid")
semiarid.plot_df_long = subset(plot_df_long, aridity_group == "Semiarid")
humid.plot_df_long = subset(plot_df_long, aridity_group == "Humid")

shapiro.test(arid.plot_df_long$aa_percent)
shapiro.test(semiarid.plot_df_long$aa_percent)
shapiro.test(humid.plot_df_long$aa_percent)

shapiro.test(arid.plot_df_long$Ip)
shapiro.test(semiarid.plot_df_long$Ip)
shapiro.test(humid.plot_df_long$Ip)

# Non-normality is confirmed across all variables

## Kruskal-Wallis tests

#Acidic AAs
acidic.plot_df_long = subset(plot_df_long, aa_type == "Acidic amino acids")
kruskal.test(aa_percent ~ aridity_group, data = acidic.plot_df_long)

#Basic AAs
basic.plot_df_long = subset(plot_df_long, aa_type == "Basic amino acids")
kruskal.test(aa_percent ~ aridity_group, data = basic.plot_df_long)

#pI
kruskal.test(mean.pI ~ aridity_group, data = plot_df_long)

## Get medians and interquartile ranges (IQR)
acidic.plot_df_long %>%
  group_by(aridity_group) %>%
  summarise(n = n(), median = median(aa_percent, na.rm = TRUE),
            Q1 = quantile(aa_percent, 0.25, na.rm = TRUE), Q3 = quantile(aa_percent, 0.75, na.rm = TRUE),
            IQR = IQR(aa_percent, na.rm = TRUE))

basic.plot_df_long %>%
  group_by(aridity_group) %>%
  summarise(n = n(), median = median(aa_percent, na.rm = TRUE),
            Q1 = quantile(aa_percent, 0.25, na.rm = TRUE), Q3 = quantile(aa_percent, 0.75, na.rm = TRUE),
            IQR = IQR(aa_percent, na.rm = TRUE))

## Histogram with distribution across aridity values
head(plot_df_long)

COLS_ARID <- c(
  "Arid"     = "#e31a1c",
  "Semiarid" = "#f89540",
  "Humid"    = "#666666"
)

breaks <- c(seq(0, 1.6, by = 0.1), 5.8)

df_binned <- plot_df_long %>%
  distinct(Proteome, Aridity, aridity_group) %>%
  mutate(
    aridity_bin = cut(
      Aridity,
      breaks = breaks,
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  count(aridity_bin, aridity_group)

p = ggplot(df_binned,
           aes(x = aridity_bin, y = n, fill = aridity_group)) +
  geom_col(color = "black", alpha = 0.6) +
  scale_fill_manual(values = COLS_ARID) +
  labs(
    x = "Aridity index (AI)",
    y = "Number of MAGs",
    fill = "Aridity category"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14))
p

tiff(file.path(OUT_DIR, "SupplFig3_distribution_MAGs.tiff"), 
     units="cm", width=16, height=10, res=300)
p
dev.off()


## Stacked barplot with distribution across phyla
COLS_ARID <- c(
  "Arid"     = "#e31a1c",
  "Semiarid" = "#f89540",
  "Humid"    = "#666666"
)

phyla_counts <- plot_df_long %>%
  distinct(Proteome, Phylum, Aridity) %>%  # remove duplicates
  mutate(
    aridity_group = case_when(
      Aridity < 0.2  ~ "Arid",
      Aridity < 0.65 ~ "Semiarid",
      TRUE           ~ "Humid"
    )
  ) %>%
  count(Phylum, aridity_group)

phyla_counts$aridity_group = factor(phyla_counts$aridity_group, levels = c("Arid", "Semiarid", "Humid"))

pd = ggplot(
  phyla_counts,
  aes(x = reorder(Phylum, n, sum), y = n, fill = aridity_group)
) +
  geom_col(alpha = 0.6) +
  coord_flip() +
  scale_fill_manual(values = COLS_ARID) +
  labs(
    x = "Phylum",
    y = "Number of MAGs",
    fill = "Aridity category"
  ) +
  theme_classic() +
  theme(
    x.axis.text = element_text(size = 10),
    y.axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

tiff(file.path(OUT_DIR, "SupplFig3B_distribution_MAGs.tiff"), 
     units="cm", width=14, height=14, res=300)
pd
dev.off()


acid.basic.percent.boxplot = ggplot(plot_df_long, aes(x = aridity_group, y = aa_percent)) +
  geom_boxplot(
    aes(fill = aridity_group),
    width = 0.55,
    alpha = 0.6,
    color = "grey20",
    outlier.shape = NA
  ) +
  
  geom_jitter(
    aes(color = aridity_group),
    width = 0.18,
    size = 1.0,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  
  facet_wrap(~ aa_type, nrow = 1, scales = "free_y") +
  
  scale_fill_manual(values = COLS_ARID, guide = "none") +
  scale_color_manual(values = COLS_ARID, guide = "none") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  
  ggsignif::geom_signif(
    comparisons = list(
      c("Arid", "Semiarid"),
      c("Arid", "Humid"),
      c("Semiarid", "Humid")
    ),
    test = "wilcox.test",
    map_signif_level = FALSE,
    textsize = 3.2,
    step_increase = 0.08,
    tip_length = 0.01,
    vjust = 0.2
  ) +
  
  labs(
    x = "Aridity preference",
    y = "Proportion of amino acids (%)"
  ) +
  
  theme_classic() +
  theme(
    strip.text = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

tiff(file.path(OUT_DIR, "Fig2B_AcidicBasicPercent_boxplot_mags.tiff"), 
     units="cm", width=16, height=10, res=300)
acid.basic.percent.boxplot
dev.off()


################################################################################
####################FIGURE 2C: Linear fit for pI with log-equation##############
################################################################################

model_log <- lm(Ip ~ log(Aridity), data = df_summarized_complete %>% filter(Aridity > 0, Aridity < 6))
intercept <- round(coef(model_log)[1], 3)
slope <- round(coef(model_log)[2], 3)
r2 <- round(summary(model_log)$r.squared, 3)
eq_text <- paste0("y = ", intercept, " + ", slope, "·ln(x)\nR² = ", r2)
summary(model_log)

plot_df <- df_summarized_complete %>%
  filter(Aridity > 0, Aridity < 6) %>%
  mutate(
    aridity_group = case_when(
      Aridity >= 0    & Aridity < 0.2  ~ "Arid",
      Aridity >= 0.2  & Aridity < 0.65 ~ "Semiarid",
      Aridity >= 0.65 & Aridity < 6    ~ "Humid"
    )
  )

plot_df$aridity_group = factor(plot_df$aridity_group, levels = c("Arid", "Semiarid", "Humid"))

pi.linear = ggplot(plot_df, aes(x = Aridity, y = Ip)) +
  geom_point(
    aes(colour = aridity_group),
    alpha = 0.4,
    size = 3
  ) +
  
  geom_smooth(
    method = "lm",
    formula = y ~ log(x),
    se = TRUE,
    colour = "black",
    fill = "black",
    alpha = 0.2
  ) +
  
  annotate(
    "text",
    x = min(plot_df$Aridity),
    y = max(plot_df$Ip) * 1.03,
    label = eq_text,
    hjust = 0,
    vjust = 1,
    size = 5,
    colour = "black"
  ) +
  
  scale_colour_manual(
    values = COLS_ARID,
    name = "Aridity class"
  ) +
  
  labs(
    x = "Aridity (AI)",
    y = "Isoelectric point (pI)"
  ) + xlim(c(0,3)) +
  
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )

tiff(file.path(OUT_DIR, "Fig2C_Ip_regression_mags.tiff"), 
     units="cm", width=12, height=10, res=300)
pi.linear
dev.off()


################################################################################
###################FIGURE 2D: Density plot for protein pI#######################
################################################################################

base_dir <- MAG_FEATURES_PARQUET_DIR

combined_df <- open_dataset(base_dir) %>%
  collect()
head(combined_df)
names(combined_df)[names(combined_df) == "proteome_id"] = "Genome"

names(df_summarized_complete)[names(df_summarized_complete) == "Proteome"] = "Genome"
head(df_summarized_complete)

dt = merge(combined_df, df_summarized_complete, by = "Genome")
head(dt)

dt <- as.data.table(dt)[, .(pI, aridity)]

dt[, aridity_cat := fifelse(
  aridity <= 0.2, "Arid",
  fifelse(aridity <= 0.65, "Semiarid", "Humid")
)]

dt$aridity_cat = factor(dt$aridity_cat, levels = c("Arid", "Semiarid", "Humid"))

#set.seed(1)
#dt_sample <- dt[, .SD[sample(.N, min(.N, 1e5))], by = aridity_cat]
#head(dt_sample)

COLS_ARID <- c(
  "Arid"         = "#e31a1c",
  "Semiarid"  = "#f89540",
  "Humid"       = "#666666"
)

dp = ggplot(dt, aes(x = pI, colour = aridity_cat)) +
  geom_density(alpha = 0.2, linewidth = 1) +
  scale_colour_manual(values = COLS_ARID, name = "Aridity") +
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

tiff(file.path(OUT_DIR, "Fig2D_density_plot_mags.tiff"), 
     units="cm", width=12, height=10, res=300)
dp
dev.off()
