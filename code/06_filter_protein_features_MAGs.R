############################################################
# TSV QUALITY FILTERING BY VALUE RANGES
############################################################
# Reads all .tsv files from IN_DIR, removes rows containing
# NA values, and applies minimum quality/range filters to
# each variable. Filtered .tsv files are then written to
# OUT_DIR using the same filenames.
############################################################

############################################################
# LIBRARIES
############################################################

library(dplyr)
library(readr)
library(tools)

############################################################
# COMMAND-LINE PARAMETERS
############################################################

# Expected usage:
# Rscript param_qc.R IN_DIR OUT_DIR

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2)
  stop("Usage: Rscript param_qc.R IN_DIR OUT_DIR", call. = FALSE)

# Store input and output paths from command-line arguments
IN_DIR <- args[1]
OUT_DIR <- args[2]

# Create output directory if it does not already exist
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Find all .tsv files inside IN_DIR. Stop if none are found.
files <- list.files(
  IN_DIR,
  pattern = "\\.tsv$",
  full.names = TRUE
)

if (length(files) == 0)
  stop("No .tsv files found in IN_DIR")

# Global counters for the final report
total_in <- 0L      # Total number of rows read
total_kept <- 0L    # Total number of rows retained
total_drop <- 0L    # Total number of rows discarded

############################################################
# MINIMUM QUALITY RANGES
############################################################

# These constants define what is considered a valid value.
# Any row that does NOT satisfy all conditions will be removed.

LEN_MIN <- 30L       # Minimum allowed sequence length
LEN_MAX <- 50000L    # Maximum allowed sequence length

GR_MIN <- -4.5       # Minimum GRAVY value
GR_MAX <- 4.5        # Maximum GRAVY value

PI_MIN <- 0          # Minimum allowed pI
PI_MAX <- 14         # Maximum allowed pI

############################################################
# Main loop: process each TSV file
############################################################

for (f in files) {
  
  # FILE READING
  #
  # read_tsv:
  #   Reads a tab-delimited file
  #   - col_guess() attempts to infer column types
  #
  # tryCatch:
  #   If reading fails, logs the error and continues
  #   with the next file.
  
  df <- tryCatch(
    read_tsv(
      f,
      col_types = cols(
        .default = col_double(),   # Numeric columns
        seq_id = col_character()   # Sequence IDs are strings
      ),
      na = c("", "NA", "NaN", ".")
    ),
    error = function(e) {
      message(
        sprintf(
          "[%-30s] READ ERROR: %s",
          basename(f),
          e$message
        )
      )
      return(NULL)
    }
  )
  
  if (is.null(df))
    next
  
  # Detect column groups by prefix:
  #
  # percent_* : amino acid composition columns
  # dis_*     : disorder-related columns
  
  pct_cols <- df %>%
    select(starts_with("percent_")) %>%
    names()
  
  dis_cols <- df %>%
    select(starts_with("dis_")) %>%
    names()
  
  dis_cols <- setdiff(
    dis_cols,
    "dis_max_run_0_5"
  )
  
  # File-level statistics
  
  n0 <- nrow(df)
  
  total_in <- total_in + n0
  
  if (n0 == 0L) {
    
    # If the file contains no data rows,
    # write it unchanged and continue.
    
    write_tsv(
      df,
      file.path(
        OUT_DIR,
        basename(f)
      )
    )
    
    next
  }
  
  ##########################################################
  # FILTER 1: REMOVE ROWS CONTAINING NA VALUES
  ##########################################################
  
  df <- df %>%
    filter(
      if_all(
        all_of(
          c(
            "seq_id",
            "length",
            "gravy",
            "pI",
            "dis_max_run_0_5",
            pct_cols,
            dis_cols
          )
        ),
        ~ !is.na(.)
      )
    )
  
  ##########################################################
  # FILTER 2: APPLY RANGE-BASED QUALITY RULES
  ##########################################################
  
  # Create logical columns *_ok indicating whether
  # each rule is satisfied.
  #
  # All rules must be TRUE for a row to be retained.
  
  keep <- df %>%
    mutate(
      
      # length within [30, 50,000]
      len_ok =
        length >= LEN_MIN &
        length <= LEN_MAX,
      
      # GRAVY within [-4.5, 4.5]
      gr_ok =
        gravy >= GR_MIN &
        gravy <= GR_MAX,
      
      # pI within [0, 14]
      pi_ok =
        pI >= PI_MIN &
        pI <= PI_MAX,
      
      # Maximum disorder run within [0, length]
      run_ok =
        dis_max_run_0_5 >= 0 &
        dis_max_run_0_5 <= length,
      
      # All percent_* columns within [0, 100]
      pct_bounds_ok =
        if (length(pct_cols)) {
          if_all(
            all_of(pct_cols),
            ~ .x >= 0 & .x <= 100
          )
        } else TRUE,
      
      # All disorder metrics except run within [0, 1]
      dis_bounds_ok =
        if (length(dis_cols)) {
          if_all(
            all_of(dis_cols),
            ~ .x >= 0 & .x <= 1
          )
        } else TRUE
      
    ) %>%
    transmute(
      ok =
        len_ok &
        gr_ok &
        pi_ok &
        run_ok &
        pct_bounds_ok &
        dis_bounds_ok
    ) %>%
    pull(ok)
  
  # Rows that pass all filters
  
  kept <- df[
    keep,
    ,
    drop = FALSE
  ]
  
  # File and global statistics
  
  n_kept <- nrow(kept)
  n_drop <- n0 - n_kept
  
  total_kept <- total_kept + n_kept
  total_drop <- total_drop + n_drop
  
  # Write filtered file to OUT_DIR
  # preserving the original filename
  
  write_tsv(
    kept,
    file.path(
      OUT_DIR,
      basename(f)
    )
  )
  
  # Per-file summary
  
  message(
    sprintf(
      "[%-30s] in:%6d kept:%6d drop:%6d",
      basename(f),
      n0,
      n_kept,
      n_drop
    )
  )
}

############################################################
# GLOBAL SUMMARY
############################################################

message(
  "------------------------------------------------------------"
)

message(
  sprintf(
    "TOTAL in:%d kept:%d drop:%d (%.6f%% discarded)",
    total_in,
    total_kept,
    total_drop,
    ifelse(
      total_in > 0,
      100 * total_drop / total_in,
      0
    )
  )
)