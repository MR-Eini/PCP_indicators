# Monthly maximum 1-day precipitation (rx1day)
# Statistic:
#   1) For each year: MAX over the 12 monthly rx1day values
#   2) For 1995–2024: MEAN of those annual maxima
#
# Input : NetCDF with monthly rx1day (typically mm)
# Output: GeoTIFF raster

library(terra)
library(ncdf4)

# ----------------------------
# Paths
# ----------------------------
base_dir <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/!NewMthMax1-dayPrecipitation"
out_dir  <- file.path(base_dir, "stats")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

nc_file <- file.path(base_dir, "rx1day_ERA5-Land_no-expt_mon_19500101-20241201.nc")
if (!file.exists(nc_file)) stop("File not found: ", nc_file)

var_name <- "rx1day"  # as per product name; script will verify

# ----------------------------
# Verify variable exists + open nc (for time)
# ----------------------------
nc <- nc_open(nc_file)
on.exit(nc_close(nc), add = TRUE)

var_names <- names(nc$var)
if (!(var_name %in% var_names)) {
  # fallback: auto-detect
  hit <- var_names[grepl("^rx1day$|rx1day|rx1", var_names, ignore.case = TRUE)]
  if (length(hit) == 0) {
    stop("Could not find rx1day in variables. Found:\n  ", paste(var_names, collapse = ", "))
  }
  var_name <- hit[1]
}
message("Using variable: ", var_name)

# ----------------------------
# Read raster
# ----------------------------
r <- rast(sprintf("NETCDF:%s:%s", nc_file, var_name))

# ERA5-Land is lon/lat; set CRS so output GeoTIFF has correct metadata
crs(r) <- "EPSG:4326"

# ----------------------------
# Ensure time is set (monthly)
# ----------------------------
get_cf_time_monthly <- function(nc, n_layers) {
  tname <- if ("time" %in% names(nc$dim)) "time" else {
    cand <- names(nc$dim)[grepl("time", names(nc$dim), ignore.case = TRUE)]
    if (length(cand) == 0) NA_character_ else cand[1]
  }
  if (is.na(tname)) return(seq(as.Date("1950-01-01"), by = "month", length.out = n_layers))
  
  tunits <- nc$dim[[tname]]$units
  origin_str <- sub(".*since\\s*", "", tunits)
  origin_str <- sub("T.*", "", origin_str)
  origin_str <- sub("\\s.*$", "", origin_str)
  origin <- suppressWarnings(as.Date(origin_str))
  if (is.na(origin)) origin <- as.Date("1950-01-01")
  
  seq(origin, by = "month", length.out = n_layers)
}

tt <- time(r)
if (is.null(tt) || all(is.na(tt))) {
  time(r) <- get_cf_time_monthly(nc, nlyr(r))
}

# ----------------------------
# Compute: mean( annual_max( monthly ) ) for 1995–2024
# ----------------------------
max_na  <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

dates <- time(r)
yrs   <- as.integer(format(dates, "%Y"))

sel <- which(yrs >= 1995 & yrs <= 2024)
if (length(sel) == 0) stop("No layers found for 1995–2024. Check time parsing and file content.")

r_sel   <- r[[sel]]
yrs_sel <- yrs[sel]  # tapp() index cannot contain NA

# 1) annual max over months
annual_max <- tapp(r_sel, index = yrs_sel, fun = max_na)

# 2) mean over years
result <- app(annual_max, fun = mean_na)
names(result) <- "rx1day_mean_of_annual_monthly_max_1995_2024"

# ----------------------------
# Save
# ----------------------------
out_tif <- file.path(out_dir, "rx1day_mean_annual_max_monthly_1995_2024.tif")
writeRaster(result, out_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
message("Saved: ", out_tif)

# ----------------------------
# Sanity checks
# ----------------------------
years_used <- sort(unique(yrs_sel))
print(c(
  total_layers      = nlyr(r),
  layers_1995_2024  = length(sel),
  years_covered     = length(years_used),
  first_year        = years_used[1],
  last_year         = years_used[length(years_used)]
))

