# Monthly very heavy precipitation days (r20)
# Statistic:
#   1) For each year: MAX over the 12 monthly values
#   2) For 1995–2024: MEAN of those annual maxima
#
# Input : NetCDF with monthly r20 (days)
# Output: GeoTIFF raster

library(terra)
library(ncdf4)

# ----------------------------
# Paths
# ----------------------------
base_dir <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/MthVeryHeavyPrecipitationDays"
out_dir  <- file.path(base_dir, "stats")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Find NetCDF
# ----------------------------
nc_files <- list.files(base_dir, recursive = TRUE, full.names = TRUE, pattern = "\\.nc$")
cand <- nc_files[grepl("(^|[^a-z0-9])r20([^a-z0-9]|$)|very[_ -]?heavy.*precip", nc_files, ignore.case = TRUE)]
if (length(cand) == 0) stop("No candidate NetCDF found for r20 under: ", base_dir)

nc_file <- cand[1]
message("Using NetCDF: ", nc_file)

# Variable name (as you confirmed)
var_name <- "r20"

# ----------------------------
# Read raster
# ----------------------------
r <- rast(sprintf("NETCDF:%s:%s", nc_file, var_name))

# ERA5-Land is lon/lat; set CRS so output GeoTIFF has correct metadata
# (If your product is not lon/lat, remove this line.)
crs(r) <- "EPSG:4326"

# ----------------------------
# Ensure time is set (monthly)
# ----------------------------
nc <- nc_open(nc_file)
on.exit(nc_close(nc), add = TRUE)

get_cf_time_monthly <- function(nc, n_layers) {
  tname <- if ("time" %in% names(nc$dim)) "time" else {
    cand <- names(nc$dim)[grepl("time", names(nc$dim), ignore.case = TRUE)]
    if (length(cand) == 0) NA_character_ else cand[1]
  }
  
  # fallback if time metadata is missing
  if (is.na(tname)) return(seq(as.Date("1950-01-01"), by = "month", length.out = n_layers))
  
  tunits <- nc$dim[[tname]]$units
  
  # parse origin from CF units like "days since 1950-01-01 ..."
  origin_str <- sub(".*since\\s*", "", tunits)
  origin_str <- sub("T.*", "", origin_str)
  origin_str <- sub("\\s.*$", "", origin_str)
  origin <- suppressWarnings(as.Date(origin_str))
  if (is.na(origin)) origin <- as.Date("1950-01-01")
  
  # for monthly data, a monthly sequence is sufficient for year grouping
  seq(origin, by = "month", length.out = n_layers)
}

tt <- time(r)
if (is.null(tt) || all(is.na(tt))) {
  time(r) <- get_cf_time_monthly(nc, nlyr(r))
}

# ----------------------------
# Compute statistic for 1995–2024
# ----------------------------
max_na  <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

dates <- time(r)
yrs   <- as.integer(format(dates, "%Y"))

sel <- which(yrs >= 1995 & yrs <= 2024)
if (length(sel) == 0) stop("No layers found for 1995–2024. Check time parsing and file content.")

r_sel   <- r[[sel]]
yrs_sel <- yrs[sel]  # no NA allowed in tapp() index

# 1) annual max across the 12 monthly layers for each year
annual_max <- tapp(r_sel, index = yrs_sel, fun = max_na)

# 2) mean of annual maxima across 1995–2024
result <- app(annual_max, fun = mean_na)
names(result) <- "r20_mean_of_annual_monthly_max_1995_2024"

# ----------------------------
# Save
# ----------------------------
out_tif <- file.path(out_dir, "r20_mean_annual_max_monthly_1995_2024.tif")
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
