# ==============================================================================
# SCRIPT: Calculate CWB1 (Climatic Water Balance) - Annual Means
# DESCRIPTION: 
#   Calculates CWB1 = Precipitation (P) - Potential Evapotranspiration (PET).
#   
#   LOGIC:
#   1. Convert "monthly mean daily rates" -> "monthly totals" (mm/month).
#   2. Calculate Monthly CWB = P - PET.
#   3. Statistic per Year: MEAN of the 12 monthly CWB values [mm/month].
#   4. Long-term Stats: MEAN, MIN, MAX of those Annual Means over the period.
#
# INPUTS: NetCDF files for Precipitation (r) and Potential Evapotranspiration (pet)
# OUTPUTS: GeoTIFFs for Mean, Min, and Max Annual Mean CWB.
# ==============================================================================

library(terra)
library(ncdf4)

# ----------------------------
# 1. Configuration (User Inputs)
# ----------------------------
# Define your period of interest
START_YEAR <- 1995L
END_YEAR   <- 2024L

# Variable Names in the NetCDF files
VAR_NAME_R   <- "r"     # Precipitation
VAR_NAME_PET <- "pet"   # Potential Evapotranspiration

# ----------------------------
# 2. File Paths
# ----------------------------
# Input Files
in_file_r   <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/!NewMthPrecipitation/r_ERA5-Land_no-expt_mon_19500101-20241201.nc"
in_file_pet <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/MthDailyAccumulatedPotentialEvapotranspiration/pet_ERA5-Land_no-expt_mon_19500101-20241201.nc"

# Output Directory
out_dir     <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/!NewMthPrecipitation/stats_CWB1_1995_2024"

# ----------------------------
# 3. Processing Flags
# ----------------------------
# TRUE: Multiplies daily rates by days-in-month to get total mm/month.
# Essential for CWB1 so that P and PET are comparable (mm vs mm) before subtraction.
FORCE_MEAN_DAILY_TO_TOTAL <- TRUE

# TRUE: Automatically flips sign if PET is stored as negative values.
AUTO_FLIP_PET_IF_NEGATIVE <- TRUE

# ----------------------------
# 4. Environment & Safety Checks
# ----------------------------
stopifnot(file.exists(in_file_r), file.exists(in_file_pet))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tmp_dir <- file.path(out_dir, "terra_tmp")
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

# optimize for memory (60% RAM usage) and redirect temp files to external drive
terraOptions(tempdir = tmp_dir, progress = 1, memfrac = 0.6)

message("Setup Complete. Processing CWB1 (Annual Means) for ", START_YEAR, "-", END_YEAR, "...")


# ==============================================================================
# 5. HELPERS
# ==============================================================================

# NetCDF Attribute Helpers
get_att <- function(nc_file, var, att) {
  nc <- nc_open(nc_file)
  on.exit(nc_close(nc), add = TRUE)
  v <- ncatt_get(nc, var, att)$value
  if (is.null(v) || is.na(v)) "" else as.character(v)
}
get_units <- function(nc_file, var) get_att(nc_file, var, "units")

# Time Parsing (Handles 'days since', 'hours since', etc.)
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

# *** VECTORIZED DATE FIX ***
# safely calculates days in month for a whole vector of dates
days_in_month <- function(dates_month_start) {
  d1 <- as.Date(format(as.Date(dates_month_start), "%Y-%m-01"))
  lt <- as.POSIXlt(d1)
  lt$mon <- lt$mon + 1
  d2 <- as.Date(lt)
  as.integer(d2 - d1)
}

seconds_in_month <- function(dates_month_start) as.numeric(days_in_month(dates_month_start)) * 86400

# Unit Conversion: Rates -> Monthly Totals
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

# Running Stats Helpers
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


# ==============================================================================
# 6. MAIN EXECUTION
# ==============================================================================

# 1. Read Files (Lazy Load)
r   <- rast(sprintf("NETCDF:%s:%s", in_file_r, VAR_NAME_R))
pet <- rast(sprintf("NETCDF:%s:%s", in_file_pet, VAR_NAME_PET))

r   <- ensure_time_monthly(r, in_file_r)
pet <- ensure_time_monthly(pet, in_file_pet)

if (is.na(crs(r)))   crs(r)   <- "EPSG:4326"
if (is.na(crs(pet))) crs(pet) <- "EPSG:4326"

# 2. Align Grids
if (!compareGeom(r, pet, stopOnError = FALSE)) {
  message("Grid differs; resampling PET to precipitation grid...")
  pet <- resample(pet, r, method = "bilinear")
  time(pet) <- time(r)
}

# 3. Time Matching
tr   <- as.Date(time(r))
tpet <- as.Date(time(pet))
key_r   <- format(tr, "%Y-%m-01")
key_pet <- format(tpet, "%Y-%m-01")

common <- intersect(key_r, key_pet)
if (length(common) == 0) stop("No overlapping months between r and pet.")
common <- sort(common)

ir   <- match(common, key_r)
ipet <- match(common, key_pet)

yrs_common <- as.integer(substr(common, 1, 4))
years <- START_YEAR:END_YEAR

if (!all(years %in% yrs_common)) {
  missing_years <- years[!(years %in% yrs_common)]
  stop("Missing years in overlapping period: ", paste(missing_years, collapse = ", "))
}

units_r   <- get_units(in_file_r, VAR_NAME_R)
units_pet <- get_units(in_file_pet, VAR_NAME_PET)

# 4. Initialize Running Statistics
sumCWB <- NULL
nCWB   <- NULL
minCWB <- NULL
maxCWB <- NULL

flip_pet <- FALSE
years_used <- integer(0)

# 5. Year-by-Year Processing Loop
for (y in years) {
  pos <- which(yrs_common == y)
  if (length(pos) != 12) stop("Year ", y, " has ", length(pos), " months (expected 12).")
  
  # Load 12 months for year Y
  r_y   <- r[[ ir[pos] ]]
  pet_y <- pet[[ ipet[pos] ]]
  
  # Convert rates to monthly totals (mm)
  # 
  
  [Image of Calendar]
  ensuring correct days per month (28, 29, 30, 31) are used
  r_y   <- to_monthly_total(r_y, units_r, FORCE_MEAN_DAILY_TO_TOTAL)
  pet_y <- to_monthly_total(pet_y, units_pet, FORCE_MEAN_DAILY_TO_TOTAL)
  
  # Auto-flip sign if PET is negative
  if (AUTO_FLIP_PET_IF_NEGATIVE && isFALSE(flip_pet) && length(years_used) == 0) {
    m_p <- as.numeric(global(pet_y[[1]], "mean", na.rm = TRUE)[1, 1])
    if (is.finite(m_p) && m_p < 0) {
      message("PET appears negative; flipping sign to make it positive.")
      flip_pet <- TRUE
    }
  }
  if (flip_pet) pet_y <- -pet_y
  
  # Calculate Monthly CWB
  cwb_y <- r_y - pet_y
  
  # --- CWB1 SPECIFIC LOGIC ---
  # Statistic: Mean of the 12 monthly values for this year
  ann_stat_y <- mean(cwb_y, na.rm = TRUE)
  ann_stat_y <- nan_to_na(ann_stat_y)
  
  # Update Long-term Stats (Mean of Means, Min of Means, Max of Means)
  if (is.null(sumCWB)) {
    sumCWB <- ifel(is.na(ann_stat_y), 0, ann_stat_y)
    nCWB   <- count_non_na(ann_stat_y)
    minCWB <- ann_stat_y
    maxCWB <- ann_stat_y
  } else {
    sumCWB <- sum_add_na0(sumCWB, ann_stat_y)
    nCWB   <- nCWB + count_non_na(ann_stat_y)
    minCWB <- cellwise_min(minCWB, ann_stat_y)
    maxCWB <- cellwise_max(maxCWB, ann_stat_y)
  }
  
  years_used <- c(years_used, y)
  
  # Cleanup RAM
  rm(r_y, pet_y, cwb_y, ann_stat_y)
  gc()
  suppressWarnings(try(tmpFiles(remove = TRUE), silent = TRUE))
  message("Processed year: ", y)
}

if (length(years_used) == 0) stop("No years processed.")

# Final Calculation
meanCWB <- ifel(nCWB == 0, NA, sumCWB / nCWB)

# ----------------------------
# 7. Write Outputs
# ----------------------------
wopt <- list(overwrite = TRUE, gdal = c("COMPRESS=LZW"))

write_r <- function(r, name, opts) {
  fname <- file.path(out_dir, sprintf("%s_%d_%d.tif", name, START_YEAR, END_YEAR))
  args <- c(list(x = r, filename = fname), opts)
  do.call(writeRaster, args)
}

# Save files: These represent the Mean/Min/Max of "Annual Mean CWB"
write_r(meanCWB, "CWB1_ann_mean_mean", wopt)
write_r(minCWB,  "CWB1_ann_mean_min",  wopt)
write_r(maxCWB,  "CWB1_ann_mean_max",  wopt)

message("Done. Years processed: ", paste(range(years_used), collapse = "-"), 
        ". Values are 'Annual Mean Monthly CWB' (mm/month).")