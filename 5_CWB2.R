# ==============================================================================
# SCRIPT: Calculate CWB2 (Climatic Water Balance) - Annual Totals
# DESCRIPTION: 
#   Calculates CWB2 = Precipitation (P) - Evapotranspiration (EV).
#   1. Converts "monthly mean daily rates" -> "monthly totals" (mm/month).
#   2. Aggregates to Annual Totals (mm/year).
#   3. Computes Long-term Stats (Mean, Min, Max) over the specified period.
#
# INPUTS: NetCDF files for Precipitation and Evapotranspiration
# OUTPUTS: GeoTIFFs for Mean, Min, and Max Annual CWB, P, and EV.
# ==============================================================================

library(terra)
library(ncdf4)

# ----------------------------
# 1. Configuration (User Inputs)
# ----------------------------
# Define your period of interest
START_YEAR <- 1995L
END_YEAR   <- 2024L

# Variable Names in the NetCDF files (Check these using nc_open if unsure)
VAR_NAME_R  <- "r"          # Precipitation
VAR_NAME_EV <- "evspsbl"    # Evapotranspiration

# ----------------------------
# 2. File Paths
# ----------------------------
# Input Files
in_file_r <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/!NewMthPrecipitation/r_ERA5-Land_no-expt_mon_19500101-20241201.nc"
in_file_e <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/!newMthEvapotransiration/evspsbl_ERA5-Land_no-expt_mon_19500101-20241201.nc"

# Output Directory
out_dir   <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/!newMthEvapotransiration/stats_CWB2_1995_2024"

# ----------------------------
# 3. Processing Flags (Logic Control)
# ----------------------------
# TRUE: Multiplies data by days-in-month (required if data is "mm/day" or "mean daily")
# FALSE: Leaves data as-is (use if data is already "mm/month")
FORCE_MEAN_DAILY_TO_TOTAL  <- TRUE 

# TRUE: Automatically flips sign if EV is stored as negative values (common in ERA5)
AUTO_FLIP_EVAP_IF_NEGATIVE <- TRUE

# ----------------------------
# 4. Environment & Safety Checks
# ----------------------------
# Validate inputs exist before running heavy code
stopifnot(file.exists(in_file_r), file.exists(in_file_e))

# Create output directories
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tmp_dir <- file.path(out_dir, "terra_tmp")
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

# Configure Terra Options
# memfrac: 0.6 uses 60% of RAM
# tempdir: Redirects temp files to your external drive (Y:) to save C: drive space
terraOptions(tempdir = tmp_dir, progress = 1, memfrac = 0.6)

message("Setup Complete. Processing ", START_YEAR, "-", END_YEAR, "...")

# ----------------------------
# Helpers
# ----------------------------
# [Standard NetCDF attribute extraction helpers omitted for brevity...]
get_att <- function(nc_file, var, att) {
  nc <- nc_open(nc_file)
  on.exit(nc_close(nc), add = TRUE)
  v <- ncatt_get(nc, var, att)$value
  if (is.null(v) || is.na(v)) "" else as.character(v)
}
get_units <- function(nc_file, var) get_att(nc_file, var, "units")
get_long_name <- function(nc_file, var) get_att(nc_file, var, "long_name")

# TIME PARSING: Robust handling of CF-convention time units (days/hours since X)
get_cf_time_as_date <- function(nc) {
  tname <- if ("time" %in% names(nc$dim)) "time" else {
    cand <- names(nc$dim)[grepl("time", names(nc$dim), ignore.case = TRUE)]
    if (length(cand) == 0) stop("No time dimension found.") else cand[1]
  }
  tvals  <- ncvar_get(nc, tname)
  tunits <- nc$dim[[tname]]$units
  
  origin_str <- sub(".*since\\s*", "", tunits)
  origin_str <- sub("T.*", "", origin_str)
  origin_str <- sub("\\s.*$", "", origin_str)
  origin <- suppressWarnings(as.Date(origin_str))
  if (is.na(origin)) origin <- as.Date("1950-01-01")
  
  unit_part <- tolower(trimws(strsplit(tunits, "since", fixed = TRUE)[[1]][1]))
  d <- if (grepl("day", unit_part)) {
    origin + tvals
  } else if (grepl("hour", unit_part)) {
    origin + (tvals / 24)
  } else if (grepl("month", unit_part)) {
    seq(origin, by = "month", length.out = length(tvals))
  } else if (grepl("year", unit_part)) {
    seq(origin, by = "year", length.out = length(tvals))
  } else {
    stop("Unhandled CF time unit: ", tunits)
  }
  as.Date(format(d, "%Y-%m-01"))
}

ensure_time_monthly <- function(x, nc_file) {
  tt <- time(x)
  if (!is.null(tt) && !all(is.na(tt))) {
    time(x) <- as.Date(format(as.Date(tt), "%Y-%m-01"))
    return(x)
  }
  nc <- nc_open(nc_file)
  on.exit(nc_close(nc), add = TRUE)
  time(x) <- get_cf_time_as_date(nc)
  x
}

# *** CRITICAL FIX ***
# This function calculates days in month safely using POSIXlt.
# Previous versions failed because seq.Date() does not accept vectors for 'from'.
days_in_month <- function(dates_month_start) {
  d1 <- as.Date(format(as.Date(dates_month_start), "%Y-%m-01"))
  lt <- as.POSIXlt(d1)
  lt$mon <- lt$mon + 1 
  d2 <- as.Date(lt)
  as.integer(d2 - d1)
}

seconds_in_month <- function(dates_month_start) as.numeric(days_in_month(dates_month_start)) * 86400

# UNIT CONVERSION:
# Multiplies rate by seconds or days in month to get TOTAL accumulation.
to_monthly_total_if_rate <- function(x, units_str) {
  u <- tolower(units_str)
  is_rate <- grepl("kg\\s*m-2\\s*s-1", u) || grepl("kg m\\-2 s\\-1", u) ||
    grepl("m\\s*s-1", u) || grepl("m s\\-1", u)
  if (!is_rate) return(x)
  
  d <- time(x)
  if (is.null(d) || all(is.na(d))) stop("Time not set; cannot convert rate to monthly total.")
  sec <- seconds_in_month(as.Date(format(as.Date(d), "%Y-%m-01")))
  
  mult_mm <- if (grepl("m\\s*s-1", u) || grepl("m s\\-1", u)) 1000 else 1
  x * (sec * mult_mm)
}

to_monthly_total <- function(x, units_str, force_mean_daily = TRUE) {
  x2 <- to_monthly_total_if_rate(x, units_str)
  if (!identical(x2, x)) return(x2)
  
  if (!force_mean_daily) return(x2)
  
  u <- tolower(units_str)
  d <- time(x2)
  if (is.null(d) || all(is.na(d))) stop("Time not set; cannot convert mean-daily to monthly total.")
  nd <- days_in_month(as.Date(format(as.Date(d), "%Y-%m-01")))
  
  mult_mm <- if (grepl("^\\s*m\\b", u) && !grepl("mm", u)) 1000 else 1
  x2 * (nd * mult_mm)
}

# STATISTICAL HELPERS:
# These handle NA values appropriately during the running update.
nan_to_na <- function(x) ifel(is.nan(x), NA, x)
cellwise_min <- function(a, b) {
  ifel(is.na(a), b, ifel(is.na(b), a, ifel(b < a, b, a)))
}
cellwise_max <- function(a, b) {
  ifel(is.na(a), b, ifel(is.na(b), a, ifel(b > a, b, a)))
}
sum_add_na0 <- function(a, b) {
  ifel(is.na(a), 0, a) + ifel(is.na(b), 0, b)
}
count_non_na <- function(x) ifel(!is.na(x), 1, 0)

# ----------------------------
# Read rasters
# ----------------------------
r <- rast(sprintf("NETCDF:%s:%s", r_file, var_r))
e <- rast(sprintf("NETCDF:%s:%s", e_file, var_e))

r <- ensure_time_monthly(r, r_file)
e <- ensure_time_monthly(e, e_file)

if (is.na(crs(r))) crs(r) <- "EPSG:4326"
if (is.na(crs(e))) crs(e) <- "EPSG:4326"

# ALIGNMENT: 
# Crucial step. If grids don't match exactly (even by rounding errors), math will fail or produce NAs.
# Here we force E to match R.
if (!compareGeom(r, e, stopOnError = FALSE)) {
  message("Grid differs; resampling evspsbl to precipitation grid...")
  e <- resample(e, r, method = "bilinear")
  time(e) <- time(r)
}

# TIME MATCHING:
# Ensures we only process months where both Precipitation and Evapotranspiration exist.
tr <- as.Date(time(r))
te <- as.Date(time(e))
key_r <- format(tr, "%Y-%m-01")
key_e <- format(te, "%Y-%m-01")

common <- intersect(key_r, key_e)
if (length(common) == 0) stop("No overlapping months between r and evspsbl.")
common <- sort(common)

ir <- match(common, key_r)
ie <- match(common, key_e)

yrs_common <- as.integer(substr(common, 1, 4))
years <- START_YEAR:END_YEAR

if (!all(years %in% yrs_common)) {
  missing_years <- years[!(years %in% yrs_common)]
  stop("These years are missing: ", paste(missing_years, collapse = ", "))
}

units_r <- get_units(r_file, var_r)
units_e <- get_units(e_file, var_e)

# ----------------------------
# Running stats initialization
# ----------------------------
# These variables will hold the cumulative state across the loop
sumP <- sumEV <- sumCWB <- NULL
nP    <- nEV    <- nCWB    <- NULL
minP <- minEV <- minCWB <- NULL
maxP <- maxEV <- maxCWB <- NULL

flip_e <- FALSE
years_used <- integer(0)

# ----------------------------
# Year-by-year Loop
# ----------------------------
for (y in years) {
  # 1. IDENTIFY: Find the 12 months belonging to year 'y'
  pos <- which(yrs_common == y)
  if (length(pos) != 12) stop("Year ", y, " incomplete.")
  
  # 2. LOAD: Read only those 12 layers into memory
  r_y <- r[[ ir[pos] ]]
  e_y <- e[[ ie[pos] ]]
  
  # 3. CONVERT: Apply the days-in-month multiplication here
  r_y <- to_monthly_total(r_y, units_r, FORCE_MEAN_DAILY_TO_TOTAL)
  e_y <- to_monthly_total(e_y, units_e, FORCE_MEAN_DAILY_TO_TOTAL)
  
  # 4. FLIP SIGN: Handle negative EV if necessary
  if (AUTO_FLIP_EVAP_IF_NEGATIVE && isFALSE(flip_e) && length(years_used) == 0) {
    m_e <- as.numeric(global(e_y[[1]], "mean", na.rm = TRUE)[1, 1])
    if (is.finite(m_e) && m_e < 0) {
      message("Flipping negative evapotranspiration.")
      flip_e <- TRUE
    }
  }
  if (flip_e) e_y <- -e_y
  
  # 5. AGGREGATE: Calculate ANNUAL TOTALS (sum of 12 months)
  # *** FIX: Using generic sum() allows Terra to work correctly ***
  P_ann  <- sum(r_y, na.rm = TRUE)
  EV_ann <- sum(e_y, na.rm = TRUE)
  
  # 6. CALCULATE CWB2: (Total Annual P) - (Total Annual EV)
  CWB_ann <- P_ann - EV_ann
  
  P_ann   <- nan_to_na(P_ann)
  EV_ann  <- nan_to_na(EV_ann)
  CWB_ann <- nan_to_na(CWB_ann)
  
  # 7. UPDATE STATS: Add current year's Annual Total to the long-term history
  if (is.null(sumP)) {
    # First year initializes the stats
    sumP   <- ifel(is.na(P_ann), 0, P_ann)
    sumEV  <- ifel(is.na(EV_ann), 0, EV_ann)
    sumCWB <- ifel(is.na(CWB_ann), 0, CWB_ann)
    
    nP    <- count_non_na(P_ann)
    nEV   <- count_non_na(EV_ann)
    nCWB  <- count_non_na(CWB_ann)
    
    minP   <- P_ann; minEV  <- EV_ann; minCWB <- CWB_ann
    maxP   <- P_ann; maxEV  <- EV_ann; maxCWB <- CWB_ann
  } else {
    # Subsequent years update the running sum/min/max
    sumP   <- sum_add_na0(sumP, P_ann)
    sumEV  <- sum_add_na0(sumEV, EV_ann)
    sumCWB <- sum_add_na0(sumCWB, CWB_ann)
    
    nP    <- nP + count_non_na(P_ann)
    nEV   <- nEV + count_non_na(EV_ann)
    nCWB  <- nCWB + count_non_na(CWB_ann)
    
    minP   <- cellwise_min(minP, P_ann)
    minEV  <- cellwise_min(minEV, EV_ann)
    minCWB <- cellwise_min(minCWB, CWB_ann)
    
    maxP   <- cellwise_max(maxP, P_ann)
    maxEV  <- cellwise_max(maxEV, EV_ann)
    maxCWB <- cellwise_max(maxCWB, CWB_ann)
  }
  
  years_used <- c(years_used, y)
  
  # CLEANUP: Crucial for keeping RAM usage low
  rm(r_y, e_y, P_ann, EV_ann, CWB_ann)
  gc()
  suppressWarnings(try(tmpFiles(remove = TRUE), silent = TRUE))
  message("Processed year: ", y)
}

# ----------------------------
# Final Calculations & Output
# ----------------------------
# Calculate final means (Sum / Count)
meanP   <- ifel(nP   == 0, NA, sumP   / nP)
meanEV  <- ifel(nEV  == 0, NA, sumEV  / nEV)
meanCWB <- ifel(nCWB == 0, NA, sumCWB / nCWB)

wopt <- list(overwrite = TRUE, gdal = c("COMPRESS=LZW"))

write_r <- function(r, name, opts) {
  fname <- file.path(out_dir, sprintf("%s_%d_%d.tif", name, START_YEAR, END_YEAR))
  args <- c(list(x = r, filename = fname), opts)
  do.call(writeRaster, args)
}

# Saving 9 files: Mean, Min, and Max of the Annual Totals
write_r(meanP,   "P_ann_mean",    wopt)
write_r(minP,    "P_ann_min",     wopt)
write_r(maxP,    "P_ann_max",     wopt)

write_r(meanEV,  "EV_ann_mean",   wopt)
write_r(minEV,   "EV_ann_min",    wopt)
write_r(maxEV,   "EV_ann_max",    wopt)

write_r(meanCWB, "CWB2_ann_mean", wopt)
write_r(minCWB,  "CWB2_ann_min",  wopt)
write_r(maxCWB,  "CWB2_ann_max",  wopt)

message("Done. Years processed: ", paste(range(years_used), collapse = "-"))