############################################################
# TSV RANGE-BASED QUALITY FILTERING
# Removes rows (proteins) with missing values
# or values outside valid ranges
############################################################

library(dplyr)
library(tools)

############################################################
# COMMAND-LINE PARAMETERS
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2)
  stop(
    "Usage: Rscript param_qc.R IN_DIR OUT_DIR",
    call. = FALSE
  )

IN_DIR <- args[1]
OUT_DIR <- args[2]

dir.create(
  OUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# MINIMUM QUALITY RANGES
############################################################

# Length
LEN_MIN <- 30L
LEN_MAX <- 50000L

# Isoelectric point
PI_MIN <- 0
PI_MAX <- 14

# acid_aa / basic_aa ??? percentages in [0, 100]
PCT_VARS <- c(
  "acid_aa",
  "basic_aa"
)

# Variables that must be within [0, 1]
FRAC_VARS <- c(
  "dis_mean",
  "dis_frac_0_5",
  "dis_top10_mean",
  "dis_N_term30",
  "dis_C_term30"
)

############################################################
# FILES TO PROCESS
############################################################

files <- list.files(
  IN_DIR,
  pattern = "\\.tsv$",
  full.names = TRUE
)

if (length(files) == 0)
  stop("No .tsv files found in IN_DIR")

total_in <- 0L
total_kept <- 0L
total_drop <- 0L

############################################################
# MAIN LOOP
############################################################

for (f in files) {
  
  # Read file
  df <- tryCatch(
    {
      df <- read.delim(
        f,
        sep = "\t",
        header = TRUE,
        na.strings = c("", "NA", "NaN", "."),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      
      # seq_id as character
      df$seq_id <- as.character(df$seq_id)
      
      # Convert all remaining columns to numeric
      num_cols <- setdiff(
        names(df),
        "seq_id"
      )
      
      df[num_cols] <- lapply(
        df[num_cols],
        as.numeric
      )
      
      df
    },
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
  
  n0 <- nrow(df)
  total_in <- total_in + n0
  
  if (n0 == 0L) {
    
    write.table(
      df,
      file = file.path(
        OUT_DIR,
        basename(f)
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    next
  }
  
  ##########################################################
  # FILTER 1: REMOVE ROWS WITH NAs
  ##########################################################
  
  if (length(num_cols) > 0)
    df <- df[
      complete.cases(
        df[, num_cols, drop = FALSE]
      ),
      ,
      drop = FALSE
    ]
  
  if (nrow(df) == 0L) {
    
    write.table(
      df,
      file = file.path(
        OUT_DIR,
        basename(f)
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    message(
      sprintf(
        "[%-30s] in:%6d kept:%6d drop:%6d",
        basename(f),
        n0,
        0L,
        n0
      )
    )
    
    total_drop <- total_drop + n0
    
    next
  }
  
  ##########################################################
  # FILTER 2: QUALITY RANGE FILTERS
  ##########################################################
  
  # length within [LEN_MIN, LEN_MAX]
  len_ok <-
    df$length >= LEN_MIN &
    df$length <= LEN_MAX
  
  # pI within [0, 14]
  pi_ok <-
    df$pI >= PI_MIN &
    df$pI <= PI_MAX
  
  # dis_max_run_0_5 within [0, length]
  run_ok <-
    df$dis_max_run_0_5 >= 0 &
    df$dis_max_run_0_5 <= df$length
  
  # acid_aa / basic_aa within [0, 100]
  pct_ok <- rep(
    TRUE,
    nrow(df)
  )
  
  for (v in PCT_VARS) {
    
    if (v %in% names(df)) {
      
      pct_ok <- pct_ok &
        df[[v]] >= 0 &
        df[[v]] <= 100
    }
  }
  
  # Disorder variables within [0, 1]
  frac_ok <- rep(
    TRUE,
    nrow(df)
  )
  
  for (v in FRAC_VARS) {
    
    frac_ok <- frac_ok &
      df[[v]] >= 0 &
      df[[v]] <= 1
  }
  
  # Final filter
  keep <-
    len_ok &
    pi_ok &
    run_ok &
    frac_ok &
    pct_ok
  
  kept <- df[
    keep,
    ,
    drop = FALSE
  ]
  
  n_kept <- nrow(kept)
  n_drop <- n0 - n_kept
  
  total_kept <- total_kept + n_kept
  total_drop <- total_drop + n_drop
  
  ##########################################################
  # WRITE OUTPUT
  ##########################################################
  
  write.table(
    kept,
    file = file.path(
      OUT_DIR,
      basename(f)
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
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
# FINAL SUMMARY
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