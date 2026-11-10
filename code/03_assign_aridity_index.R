## clean/reset environment 
rm(list=ls())

############################################################
# CONFIGURATION: edit these paths before running
############################################################
AI_FILE        <- "/YOURDIRECTORY/Global-AI_ET0__annual_v3_1/ai_v31_yr.tif"
MAGS_COORDS    <- "/YOURDIRECTORY/MAGs_with_lat.tsv"  # accession + ncbi_lat_lon (from soil_metadata.tsv)
MAGS_TABLE     <- "/YOURDIRECTORY/MAGs.tsv"           # MAG table to which aridity is added
MAGS_TABLE_OUT <- "/YOURDIRECTORY/MAGs_aridity.tsv"   # output

# Load required packages
library(raster)
library(dplyr)
library(stringr)

# Function to convert coordinate strings to decimal degrees
parse_coords <- function(coord_str) {
  parts <- str_split(coord_str, " ", simplify = TRUE)
  lat_val <- as.numeric(parts[1])
  lat_dir <- parts[2]
  lon_val <- as.numeric(parts[3])
  lon_dir <- parts[4]
  
  lat <- ifelse(lat_dir == "N", lat_val, -lat_val)
  lon <- ifelse(lon_dir == "E", lon_val, -lon_val)
  
  return(data.frame(latitude = lat, longitude = lon))
}

# Load Global Aridity Database raster (download first from: https://cgiarcsi.community/data/global-aridity-index/)
# Replace with your actual file path
gad_raster <- raster(AI_FILE)

# Example data frame with coordinates
# Replace this with your actual data
MAGs_with_lat <- read.delim(MAGS_COORDS)
data <- MAGs_with_lat %>%
  transmute(
    id    = accession,
    coords = ncbi_lat_lon
  )

# Process coordinates and extract aridity values
result <- data %>%
  mutate(coord_data = lapply(coords, parse_coords)) %>%
  tidyr::unnest(coord_data) %>%
  mutate(
    aridity = raster::extract(gad_raster, cbind(longitude, latitude))
  )

# Save aridity and ids
mags <- read.delim(MAGS_TABLE)
out  <- merge(mags, result[c("id","aridity")], by.x = "accession", 
              by.y = "id", all.x = TRUE)
write.table(out, MAGS_TABLE_OUT, sep = "\t", quote = FALSE, 
            row.names = FALSE)

