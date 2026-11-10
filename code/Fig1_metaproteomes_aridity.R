rm(list = ls())

############################################################
# CONFIGURATION: edit these paths before running
############################################################
# Run from the repository root so that "data/" resolves.
DATA_DIR <- "data"                   # files included in this repository
EXT_DIR  <- "/YOURDIRECTORY"         # inputs not included here (see README)
OUT_DIR  <- "figures"                # figures are written here
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
METAPROTEOME_FEATURES_DIR <- file.path(EXT_DIR, "Metaproteome_protparams")
METAPROTEOME_PI_FILE <- file.path(EXT_DIR, "Data_metaproteomes_parquet", "pI.parquet")
METAGENOME_METADATA_FILE <- file.path(DATA_DIR, "metagenomes_metadata.tsv")
AI_FILE <- file.path(EXT_DIR, "Global-AI_ET0__annual_v3_1", "ai_v31_yr.tif")

library(data.table)
library(dplyr)
library(arrow)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(stringr)
library(tidyr)

#####################################
########### METAPROTEOMES ###########
#####################################
# Get list of TSV files
file_paths <- list.files(path = METAPROTEOME_FEATURES_DIR, 
                         pattern = "\\.tsv$", 
                         full.names = TRUE)

# Read all files and add filename column
combined_data <- rbindlist(lapply(file_paths, function(file) {
  dt <- fread(file)
  dt[, source_file := basename(file)]
}), fill = TRUE)

head(combined_data)

combined_data$metagenome_id = gsub(".protparam.tsv", "", combined_data$source_file)
head(combined_data)

# Get true percentage of acidic and basic AAs per proteome
proteome.aa <- combined_data[, .(
  acid.aa.percent  = sum(acid_aa  * length / 100) / sum(length) * 100,
  basic.aa.percent = sum(basic_aa * length / 100) / sum(length) * 100,
  mean.pI = mean(pI)), by = metagenome_id]
head(proteome.aa)

# Read the aridity data and metadata of the metaproteomes
pi.parquet = read_parquet(METAPROTEOME_PI_FILE)
head(pi.parquet)

# Merge the datasets by metagenome id
full.dataset = merge(proteome.aa, pi.parquet, by = "metagenome_id")
full.dataset$acid.basic.percent.ratio = full.dataset$acid.aa.percent / full.dataset$basic.aa.percent
full.dataset$acid.basic.aa.ratio = full.dataset$mean.acid.aa.percent / full.dataset$mean.basic.aa.percent
head(full.dataset)

pi.model = lm(aridity ~ mean.pI, data = full.dataset)
AIC(pi.model) # AIC = 290.69
acid.percent.model = lm(aridity ~ acid.aa.percent, data = full.dataset)
AIC(acid.percent.model) # AIC = 285.13 BEST MODEL FOR AA-SEQUENCE LEVEL ANALYSIS***
basic.percent.model = lm(aridity ~ basic.aa.percent, data = full.dataset)
AIC(basic.percent.model) # AIC = 331.51

################################################################################
##############FIGURE 1B: Facetted boxplot for acidic and basic AAs##############
################################################################################

COLS_ARID <- c(
  "Arid"     = "#e31a1c",
  "Semiarid" = "#f89540",
  "Humid"    = "#666666"
)

plot_df_long <- full.dataset %>%
  mutate(
    aridity_group = case_when(
      aridity >= 0    & aridity < 0.2  ~ "Arid",
      aridity >= 0.2  & aridity < 0.65 ~ "Semiarid",
      aridity >= 0.65 & aridity < 6    ~ "Humid"
    )
  ) %>%
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

shapiro.test(arid.plot_df_long$mean.pI)
shapiro.test(semiarid.plot_df_long$mean.pI)
shapiro.test(humid.plot_df_long$mean.pI)

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

breaks <- c(seq(0, 1.6, by = 0.1), 3)

df_binned <- plot_df_long %>%
  distinct(metagenome_id, aridity, aridity_group) %>%
  mutate(
    aridity_bin = cut(
      aridity,
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
    y = "Number of metagenomes",
    fill = "Aridity category"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14))
p

tiff(file.path(OUT_DIR, "SupplFig1_distribution_metaproteomes.tiff"), 
     units="cm", width=16, height=10, res=300)
p
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

tiff(file.path(OUT_DIR, "Fig1B_AcidicBasicPercent_boxplot_metaproteomes.tiff"), 
     units="cm", width=16, height=10, res=300)
acid.basic.percent.boxplot
dev.off()


################################################################################
####################FIGURE 1C: Linear fit for pI with log-equation##############
################################################################################

model_log <- lm(mean.pI ~ log(aridity), data = full.dataset %>% filter(aridity > 0, aridity < 6))
intercept <- round(coef(model_log)[1], 3)
slope <- round(coef(model_log)[2], 3)
r2 <- round(summary(model_log)$r.squared, 3)
eq_text <- paste0("y = ", intercept, " + ", slope, "·ln(x)\nR² = ", r2)

plot_df <- full.dataset %>%
  filter(aridity > 0, aridity < 6) %>%
  mutate(
    aridity_group = case_when(
      aridity >= 0    & aridity < 0.2  ~ "Arid",
      aridity >= 0.2  & aridity < 0.65 ~ "Semiarid",
      aridity >= 0.65 & aridity < 6    ~ "Humid"
    )
  )

plot_df$aridity_group = factor(plot_df$aridity_group, levels = c("Arid", "Semiarid", "Humid"))

pi.linear = ggplot(plot_df, aes(x = aridity, y = mean.pI)) +
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
    x = min(plot_df$aridity) + 0.01,
    y = max(plot_df$mean.pI) * 1.03,
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
  ) +
  
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )

tiff(file.path(OUT_DIR, "Fig1C_Ip_regression_metaproteomes.tiff"), 
     units="cm", width=12, height=10, res=300)
pi.linear
dev.off()


################################################################################
###################FIGURE 1D: Density plot for protein pI#######################
################################################################################

ratio.dataset = merge(combined_data, pi.parquet, by = "metagenome_id")
head(ratio.dataset)

dt <- as.data.table(ratio.dataset)[, .(pI, aridity)]

dt[, aridity_cat := fifelse(
  aridity <= 0.2, "Arid",
  fifelse(aridity <= 0.65, "Semiarid", "Humid")
)]

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

tiff(file.path(OUT_DIR, "Fig1D_density_plot_metaproteomes.tiff"), 
     units="cm", width=12, height=10, res=300)
dp
dev.off()


############################################################
# Map of Australia - Aridity Index + metagenome points
############################################################
rm(list = ls())

# LIBRARIES
library(terra)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(rnaturalearth)    # coastlines / country borders
library(rnaturalearthdata)

############################################################
# CONFIGURATION
############################################################
POINTS_FILE <- METAGENOME_METADATA_FILE
OUT_PNG   <- file.path(OUT_DIR, "map_metagenomes.tif")

# Bounding box for Australia
AUS_EXT <- terra::ext(110, 155, -45, -10)

ai_rast <- terra::rast(AI_FILE) |>
  terra::crop(AUS_EXT) |>
  terra::aggregate(fact = 3, fun = "mean", na.rm = TRUE)

ai_rast <- terra::rast(AI_FILE) |>
  terra::crop(AUS_EXT)

# Mask oceans using Australia polygon
aus_vect <- terra::vect(rnaturalearth::ne_countries(
  country = "Australia", scale = "medium", returnclass = "sf"))

ai_rast <- terra::mask(ai_rast, aus_vect) |>
  terra::aggregate(fact = 3, fun = "mean", na.rm = TRUE)

# Convert to dataframe (NOW oceans are NA and removed)
ai_df <- as.data.frame(ai_rast, xy = TRUE, na.rm = TRUE) |>
  dplyr::rename(ai = 3) |>
  dplyr::mutate(
    ai = ai / 10000,
    category = dplyr::case_when(
      is.infinite(ai) | ai >= 0.65 ~ "Humid",
      ai >= 0.50 & ai < 0.65       ~ "Dry sub-humid",
      ai >= 0.20 & ai < 0.50       ~ "Semi-arid",
      ai >= 0.05 & ai < 0.20       ~ "Arid",
      ai <  0.05                   ~ "Hyper-arid",
      TRUE                         ~ NA_character_
    ),
    category = factor(category,
                      levels = c("Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid", "Humid"))
  ) |>
  dplyr::filter(!is.na(category))

############################################################
# SAMPLE POINTS
############################################################
points_df <- readr::read_tsv(POINTS_FILE, show_col_types = FALSE) |>
  tidyr::separate(lat_lon, into = c("lat", "lon"), sep = "\\s+", convert = TRUE)

############################################################
# COASTLINE / BORDERS
############################################################
aus_coast <- rnaturalearth::ne_countries(
  country = "Australia", scale = "medium", returnclass = "sf")

aus_states <- rnaturalearth::ne_states(
  country = "Australia", returnclass = "sf")

############################################################
# MAP
############################################################
aridity_colors <- c(
  "Hyper-arid"    = "#e31a1c",
  "Arid"          = "#fd8d3c",
  "Semi-arid"     = "#fecc5c",
  "Dry sub-humid" = "#ffffb2",
  "Humid"         = "#666666"
)

p_map <- ggplot() +
  # Aridity raster background
  geom_tile(data = ai_df, aes(x = x, y = y, fill = category)) +
  scale_fill_manual(
    values = aridity_colors,
    name   = "Aridity class",
    guide  = guide_legend(reverse = FALSE), na.value = NA
  ) +
  # State borders (thin, dashed)
  geom_sf(data = aus_states, fill = NA, colour = "grey50",
          linewidth = 0.25, linetype = "dashed", inherit.aes = FALSE) +
  # Country coastline on top (thicker)
  geom_sf(data = aus_coast, fill = NA, colour = "grey20",
          linewidth = 0.5, inherit.aes = FALSE) +
  # Metagenome sample points
  geom_jitter(data = points_df, aes(x = lon, y = lat),
              shape = 21, fill = "white", colour = "black",
              size = 3, width = 0.5, height = 1) +
  # Map extent
  coord_sf(xlim = c(110, 155), ylim = c(-45, -10), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    legend.title       = element_text(face = "bold", size = 10),
    legend.text        = element_text(size = 9),
    axis.title         = element_text(size = 12),
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.background   = element_rect(fill = "white")
  )
p_map

ggsave(
  filename = file.path(OUT_DIR, "Fig1A_map_metagenomes.tiff"),
  plot     = p_map,
  width    = 12, height = 10,
  units    = "cm", dpi = 300,
  compression = "lzw"
)
