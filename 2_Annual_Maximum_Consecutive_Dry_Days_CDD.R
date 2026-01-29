# Annual Maximum Consecutive Dry Days (CDD)
# Task: average (mean) for the 1995–2024 period
#
# Input : NetCDF with annual CDD (days)
# Output: GeoTIFF raster (mean CDD over 1995–2024)

library(terra)
library(ncdf4)

# ----------------------------
# Paths
# ----------------------------
base_dir <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/AnnualMaximumConsecutiveDryDays"
out_dir  <- file.path(base_dir, "stats")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# If you want to hardcode the file:
nc_file <- file.path(base_dir, "cdd_ERA5-Land_no-expt_yr_19500101-20240101.nc")
if (!file.exists(nc_file)) {
  # fallback: find any nc in folder
  nc_files <- list.files(base_dir, recursive = TRUE, full.names = TRUE, pattern = "\\.nc$")
  stopifnot(length(nc_files) > 0)
  nc_file <- nc_files[1]
}
message("Using NetCDF: ", nc_file)

# ----------------------------
# Detect variable name (prefer 'cdd')
# ----------------------------
nc <- nc_open(nc_file)
on.exit(nc_close(nc), add = TRUE)

var_names <- names(nc$var)
var_name <- var_names[grepl("^cdd$|cdd", var_names, ignore.case = TRUE)][1]
if (is.na(var_name) || !nzchar(var_name)) {
  stop("Could not detect a CDD variable. Variables found:\n  ",
       paste(var_names, collapse = ", "),
       "\nSet var_name manually.")
}
message("Using variable: ", var_name)

# ----------------------------
# Read raster
# ----------------------------
r <- rast(sprintf("NETCDF:%s:%s", nc_file, var_name))

# ERA5-Land is lon/lat; set CRS so output GeoTIFF has correct metadata
crs(r) <- "EPSG:4326"

# ----------------------------
# Ensure time is set (annual)
# ----------------------------
get_cf_time_yearly <- function(nc, n_layers) {
  tname <- if ("time" %in% names(nc$dim)) "time" else {
    cand <- names(nc$dim)[grepl("time", names(nc$dim), ignore.case = TRUE)]
    if (length(cand) == 0) NA_character_ else cand[1]
  }
  if (is.na(tname)) return(seq(as.Date("1950-01-01"), by = "year", length.out = n_layers))
  
  tunits <- nc$dim[[tname]]$units
  origin_str <- sub(".*since\\s*", "", tunits)
  origin_str <- sub("T.*", "", origin_str)
  origin_str <- sub("\\s.*$", "", origin_str)
  origin <- suppressWarnings(as.Date(origin_str))
  if (is.na(origin)) origin <- as.Date("1950-01-01")
  
  seq(origin, by = "year", length.out = n_layers)
}

tt <- time(r)
if (is.null(tt) || all(is.na(tt))) {
  time(r) <- get_cf_time_yearly(nc, nlyr(r))
}

# ----------------------------
# Mean for 1995–2024 (use available years in file)
# ----------------------------
mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

dates <- time(r)
yrs   <- as.integer(format(dates, "%Y"))

sel <- which(yrs >= 1995 & yrs <= 2024)
if (length(sel) == 0) stop("No layers found for 1995–2024. Check time parsing and file content.")

yrs_used <- sort(unique(yrs[sel]))
message("Years used: ", yrs_used[1], "–", yrs_used[length(yrs_used)],
        " (n=", length(yrs_used), ")")

result <- app(r[[sel]], fun = mean_na)
names(result) <- "cdd_mean_1995_2024"

# ----------------------------
# Save
# ----------------------------
out_tif <- file.path(out_dir, "cdd_mean_1995_2024.tif")
writeRaster(result, out_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
message("Saved: ", out_tif)

# Optional: single summary number (global mean over raster cells)
# global(result, "mean", na.
