# CWB1 = Monthly precipitation (r) - Monthly daily accumulated potential evapotranspiration (pet)
# Statistic:
#   1) For each year: MEAN over the 12 monthly CWB1 values
#   2) For 1994–2024: MEAN of those annual means
#
# Output: GeoTIFF raster

library(terra)
library(ncdf4)

# ----------------------------
# Files
# ----------------------------
r_file   <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/!NewMthPrecipitation/r_ERA5-Land_no-expt_mon_19500101-20241201.nc"
pet_file <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/MthDailyAccumulatedPotentialEvapotranspiration/pet_ERA5-Land_no-expt_mon_19500101-20241201.nc"
stopifnot(file.exists(r_file), file.exists(pet_file))

out_dir <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/!NewMthPrecipitation/stats"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tmp_dir <- file.path(out_dir, "terra_tmp")
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

# Use a temp directory on a fast drive with space
terraOptions(tempdir = tmp_dir, progress = 1, memfrac = 0.6)

# ----------------------------
# Helpers
# ----------------------------
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

get_units <- function(nc_file, var) {
  nc <- nc_open(nc_file)
  on.exit(nc_close(nc), add = TRUE)
  u <- ncatt_get(nc, var, "units")$value
  if (is.null(u) || is.na(u)) "" else as.character(u)
}

seconds_in_month <- function(dates_month_start) {
  d1 <- as.Date(format(dates_month_start, "%Y-%m-01"))
  d2 <- seq(d1, by = "month", length.out = length(d1) + 1)[-1]
  as.numeric(d2 - d1) * 86400
}

to_monthly_total_if_rate <- function(x, units_str) {
  u <- tolower(units_str)
  is_rate <- grepl("kg\\s*m-2\\s*s-1", u) || grepl("kg m\\-2 s\\-1", u) ||
    grepl("m\\s*s-1", u)       || grepl("m s\\-1", u)
  if (!is_rate) return(x)
  
  d <- time(x)
  if (is.null(d) || all(is.na(d))) stop("Time not set; cannot convert rate to monthly total.")
  sec <- seconds_in_month(as.Date(format(as.Date(d), "%Y-%m-01")))
  
  mult <- if (grepl("m\\s*s-1", u) || grepl("m s\\-1", u)) 1000 else 1
  x * (sec * mult)  # kg m-2 s-1 -> mm/month ; m s-1 -> mm/month
}

mean_na <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)

# ----------------------------
# Read rasters (do NOT subset big blocks yet)
# ----------------------------
r   <- rast(sprintf("NETCDF:%s:%s", r_file, "r"))
pet <- rast(sprintf("NETCDF:%s:%s", pet_file, "pet"))

# ERA5-Land lon/lat
crs(r) <- "EPSG:4326"
crs(pet) <- "EPSG:4326"

r   <- ensure_time_monthly(r, r_file)
pet <- ensure_time_monthly(pet, pet_file)

# Align grids once (if needed)
if (!compareGeom(r, pet, stopOnError = FALSE)) {
  message("Grid differs; resampling PET to precipitation grid...")
  pet <- resample(pet, r, method = "bilinear")
  time(pet) <- time(r)
}

# Match months by YYYY-MM-01
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
years <- 1994:2024
years <- years[years %in% yrs_common]
if (length(years) == 0) stop("No overlapping years in 1994–2024 after matching time axes.")

# Units (optional conversion if variables are rates)
units_r   <- get_units(r_file, "r")
units_pet <- get_units(pet_file, "pet")
message("Units r   : ", units_r)
message("Units pet : ", units_pet)

# Directory for annual mean rasters
ann_dir <- file.path(out_dir, "CWB1_annual_means")
dir.create(ann_dir, showWarnings = FALSE, recursive = TRUE)

# If PET is negative by convention, flip it automatically (based on first matched month)
AUTO_FLIP_PET_IF_NEGATIVE <- TRUE

# ----------------------------
# Year-by-year processing (memory-safe)
# ----------------------------
annual_files <- character(0)

for (y in years) {
  pos <- which(yrs_common == y)          # positions within 'common'
  if (length(pos) == 0) next
  
  # Read only the 12 months for that year (or however many exist)
  r_y   <- r[[ ir[pos] ]]
  pet_y <- pet[[ ipet[pos] ]]
  
  # Convert to monthly totals if stored as rates
  r_y   <- to_monthly_total_if_rate(r_y, units_r)
  pet_y <- to_monthly_total_if_rate(pet_y, units_pet)
  
  if (AUTO_FLIP_PET_IF_NEGATIVE && y == years[1]) {
    m_pet <- as.numeric(global(pet_y[[1]], "mean", na.rm = TRUE)[1, 1])
    if (is.finite(m_pet) && m_pet < 0) {
      message("PET appears negative; flipping sign to make it positive (first layer mean = ", signif(m_pet, 5), ").")
      AUTO_FLIP_PET_IF_NEGATIVE <- FALSE  # only decide once
      flip_pet <- TRUE
    } else {
      flip_pet <- FALSE
    }
  }
  if (exists("flip_pet") && isTRUE(flip_pet)) pet_y <- -pet_y
  
  # Monthly CWB1 then annual mean
  cwb_y <- r_y - pet_y
  ann_mean_y <- app(cwb_y, fun = mean_na)
  
  out_y <- file.path(ann_dir, sprintf("CWB1_annual_mean_%d.tif", y))
  writeRaster(ann_mean_y, out_y, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  annual_files <- c(annual_files, out_y)
  
  rm(r_y, pet_y, cwb_y, ann_mean_y)
  gc()
}

if (length(annual_files) == 0) stop("No annual mean files were created.")

# ----------------------------
# Mean of annual means (1994–2024)
# ----------------------------
ann_stack <- rast(annual_files)
result <- app(ann_stack, fun = mean_na)
names(result) <- "CWB1_mean_of_annual_mean_1994_2024"

out_tif <- file.path(out_dir, "CWB1_mean_annual_mean_1994_2024.tif")
writeRaster(result, out_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
message("Saved: ", out_tif)

# Sanity checks
print(c(
  annual_layers_created = nlyr(ann_stack),
  first_year = as.integer(sub(".*_(\\d{4})\\.tif$", "\\1", annual_files[1])),
  last_year  = as.integer(sub(".*_(\\d{4})\\.tif$", "\\1", annual_files[length(annual_files)]))
))
