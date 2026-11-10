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

library(data.table)
library(dplyr)
library(arrow)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(stringr)

##########################################################################################
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

combined_data$source_file <- gsub("\\.protparam\\.tsv$", "", combined_data$source_file)
names(combined_data)[names(combined_data) == "source_file"] <- "metagenome_id" 
head(combined_data)

# Read the aridity data and metadata of the metaproteomes
pi.parquet = read_parquet(METAPROTEOME_PI_FILE)
head(pi.parquet)

# Merge the datasets by metagenome id
ratio.dataset = merge(combined_data, pi.parquet, by = "metagenome_id")
rm(combined_data)
head(ratio.dataset)

# Plot
dt <- as.data.table(ratio.dataset)[, .(pI, aridity)]

dt[, aridity_cat := fifelse(
  aridity <= 0.2, "Arid",
  fifelse(aridity <= 0.65, "Semiarid", "Humid")
)]

set.seed(1)
dt_sample <- dt[, .SD[sample(.N, min(.N, 1e15))], by = aridity_cat]
head(dt_sample)

# Density estimates on a common grid
d_arid <- density(dt_sample$pI[dt_sample$aridity_cat == "Arid"],
                  from = min(dt_sample$pI),
                  to   = max(dt_sample$pI),
                  n = 1024)

d_humid <- density(dt_sample$pI[dt_sample$aridity_cat == "Humid"],
                   from = min(dt_sample$pI),
                   to   = max(dt_sample$pI),
                   n = 1024)

delta_df <- data.frame(
  pI = d_arid$x,
  delta = d_arid$y - d_humid$y
)

delta.pi = ggplot(delta_df, aes(pI, delta)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(color = "#e31a1c", linewidth = 1) +
  labs(
    x = "pI",
    y = expression(Delta * Density ~ "(Arid - Humid)")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

tiff(file.path(OUT_DIR, "Supplementary_Delta_Density_metaproteomes.tiff"), 
     units="cm", width=12, height=10, res=300)
delta.pi
dev.off()

library(pracma)

# Positive peaks
pos_peaks <- findpeaks(
  delta_df$delta,
  minpeakheight = 0.01,
  minpeakdistance = 20
)

pos_peaks <- data.frame(
  height = pos_peaks[,1],
  index  = pos_peaks[,2],
  pI_peak = delta_df$pI[pos_peaks[,2]]
)

# Negative peaks
neg_peaks <- findpeaks(
  -delta_df$delta,
  minpeakheight = 0.01,
  minpeakdistance = 20
)

neg_peaks <- data.frame(
  height = neg_peaks[,1],
  index  = neg_peaks[,2],
  pI_peak = delta_df$pI[neg_peaks[,2]]
)

selected_pos_peaks = pos_peaks

selected_neg_peaks <-
  neg_peaks %>%
  mutate(
    mode = case_when(
      pI_peak < 7    ~ "acidic",
      pI_peak < 9    ~ "neutral",
      TRUE           ~ "basic"
    )
  ) %>%
  group_by(mode) %>%
  slice_max(height, n = 1, with_ties = FALSE) %>%
  arrange(mode, desc(height))

selected_neg_peaks

selected_peaks = rbind(selected_neg_peaks, selected_pos_peaks)
selected_peaks

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

tiff(file.path(OUT_DIR, "Density_plot_metaproteomes.tiff"), 
     units="cm", width=12, height=10, res=300)
dp
dev.off()
