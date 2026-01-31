#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ncdf4)
  library(terra)
})

# Avoid mixing raster + terra I/O
if ("package:raster" %in% search()) {
  detach("package:raster", unload = TRUE, character.only = TRUE)
}

# ----------------------------
# USER SETTINGS
# ----------------------------
PERIOD_START <- 1995L
PERIOD_END   <- 2024L

MODE <- "CWB2"  # "CWB1" (r - pet) or "CWB2" (r - evspsbl)

r_file <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/!NewMthPrecipitation/r_ERA5-Land_no-expt_mon_19500101-20241201.nc"
pet_file <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/MthDailyAccumulatedPotentialEvapotranspiration/pet_ERA5-Land_no-expt_mon_19500101-20241201.nc"
evspsbl_file <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/!newMthEvapotransiration/evspsbl_ERA5-Land_no-expt_mon_19500101-20241201.nc"

# Copy final outputs to Y: (set FALSE if you only want local outputs)
COPY_TO_Y <- TRUE
out_base_Y <- "Y:/Maps/climate/CDS_CopernicusInteractiveClimateAtlas/_CWB_stats"

SAVE_ANNUAL_COMPONENTS <- TRUE
CLAMP_NEGATIVE_EV_TO_ZERO <- TRUE
READ_FAIL_ACTION <- "stop"   # "stop" | "skip_month" | "set_zero"

# Local work (recommended)
LOCAL_WORK_DIR  <- "C:/temp/CWB_work"
LOCAL_NC_CACHE  <- file.path(LOCAL_WORK_DIR, "nc_cache")
LOCAL_OUT_BASE  <- file.path(LOCAL_WORK_DIR, "out")
LOCAL_TERRA_TMP <- file.path(LOCAL_WORK_DIR, "terra_tmp")
dir.create(LOCAL_NC_CACHE, recursive = TRUE, showWarnings = FALSE)
dir.create(LOCAL_OUT_BASE, recursive = TRUE, showWarnings = FALSE)
dir.create(LOCAL_TERRA_TMP, recursive = TRUE, showWarnings = FALSE)

terraOptions(tempdir = LOCAL_TERRA_TMP, progress = 1, memfrac = 0.6)

stopifnot(file.exists(r_file))
if (MODE == "CWB1") stopifnot(file.exists(pet_file))
if (MODE == "CWB2") stopifnot(file.exists(evspsbl_file))

RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ----------------------------
# Helpers
# ----------------------------
localize_nc <- function(src) {
  src_n <- normalizePath(src, winslash="/", mustWork=TRUE)
  dst   <- file.path(LOCAL_NC_CACHE, basename(src_n))
  if (!file.exists(dst) || file.info(dst)$size != file.info(src_n)$size) {
    ok <- file.copy(src_n, dst, overwrite = TRUE)
    if (!ok) stop("Failed to copy NetCDF to local cache: ", dst)
  }
  normalizePath(dst, winslash="/", mustWork=TRUE)
}

cf_time_to_month_start <- function(nc, fallback_origin = as.Date("1950-01-01")) {
  tname <- if ("time" %in% names(nc$dim)) "time" else {
    cand <- names(nc$dim)[grepl("time", names(nc$dim), ignore.case = TRUE)]
    if (length(cand)==0) stop("No time dimension found.") else cand[1]
  }
  tvals  <- ncvar_get(nc, tname)
  tunits <- nc$dim[[tname]]$units
  
  origin_str <- sub(".*since\\s*", "", tunits)
  origin_str <- sub("T.*", "", origin_str)
  origin_str <- sub("\\s.*$", "", origin_str)
  origin <- suppressWarnings(as.Date(origin_str))
  if (is.na(origin)) origin <- fallback_origin
  
  unit_part <- tolower(trimws(strsplit(tunits, "since", fixed=TRUE)[[1]][1]))
  d <- if (grepl("day", unit_part)) {
    origin + tvals
  } else if (grepl("hour", unit_part)) {
    origin + (tvals/24)
  } else if (grepl("month", unit_part)) {
    seq(origin, by="month", length.out=length(tvals))
  } else if (grepl("year", unit_part)) {
    seq(origin, by="year", length.out=length(tvals))
  } else {
    seq(fallback_origin, by="month", length.out=length(tvals))
  }
  d <- as.Date(format(d, "%Y-%m-01"))
  if (any(is.na(d))) d <- seq(fallback_origin, by="month", length.out=length(tvals))
  d
}

days_in_month_vec <- function(dates_month_start) {
  d1 <- as.Date(format(as.Date(dates_month_start), "%Y-%m-01"))
  y  <- as.integer(format(d1, "%Y"))
  m  <- as.integer(format(d1, "%m"))
  y2 <- y + (m == 12)
  m2 <- ifelse(m == 12, 1L, m + 1L)
  d2 <- as.Date(sprintf("%04d-%02d-01", y2, m2))
  as.integer(d2 - d1)
}

seconds_in_month_vec <- function(dates_month_start) {
  d1 <- as.Date(format(as.Date(dates_month_start), "%Y-%m-01"))
  d2 <- seq(d1, by="month", length.out=length(d1)+1)[-1]
  as.numeric(d2 - d1) * 86400
}

get_dim_info <- function(nc, var) {
  dims <- nc$var[[var]]$dim
  dn <- tolower(vapply(dims, `[[`, "", "name"))
  lon_id  <- which(dn %in% c("lon","longitude","x","rlon","long"))
  lat_id  <- which(dn %in% c("lat","latitude","y","rlat"))
  time_id <- which(grepl("time", dn))
  if (length(lon_id)!=1 || length(lat_id)!=1 || length(time_id)!=1) {
    stop("Cannot identify lon/lat/time dims for ", var, " dims=", paste(dn, collapse=","))
  }
  list(
    dn=dn, lon_id=lon_id, lat_id=lat_id, time_id=time_id,
    lon_name=dims[[lon_id]]$name, lat_name=dims[[lat_id]]$name, time_name=dims[[time_id]]$name,
    lon_n=dims[[lon_id]]$len, lat_n=dims[[lat_id]]$len, time_n=dims[[time_id]]$len
  )
}

get_att <- function(nc, var, att) {
  v <- ncatt_get(nc, var, att)$value
  if (is.null(v) || is.na(v)) "" else as.character(v)
}
get_units <- function(nc, var) get_att(nc, var, "units")
get_long_name <- function(nc, var) get_att(nc, var, "long_name")

get_fill_value <- function(nc, var) {
  fv <- ncatt_get(nc, var, "_FillValue")$value
  if (is.null(fv) || is.na(fv)) fv <- ncatt_get(nc, var, "missing_value")$value
  if (is.null(fv) || is.na(fv)) NA_real_ else as.numeric(fv)
}

monthly_multiplier <- function(units_str, long_name_str, month_dates) {
  u  <- tolower(units_str)
  ln <- tolower(long_name_str)
  
  is_rate <- grepl("kg\\s*m-2\\s*s-1", u) || grepl("kg m\\-2 s\\-1", u) ||
    grepl("m\\s*s-1", u)         || grepl("m s\\-1", u)
  
  if (is_rate) {
    sec <- seconds_in_month_vec(month_dates)
    mult <- sec
    if (grepl("m\\s*s-1", u) || grepl("m s\\-1", u)) mult <- mult * 1000
    return(mult)
  }
  
  if (grepl("monthly mean of daily accumulated", ln)) {
    return(days_in_month_vec(month_dates))
  }
  
  rep(1, length(month_dates))
}

make_template_from_lonlat <- function(lon, lat) {
  lon_sort <- sort(lon); lat_sort <- sort(lat)
  resx <- median(diff(lon_sort)); resy <- median(diff(lat_sort))
  rast(
    ncols = length(lon), nrows = length(lat),
    ext = ext(min(lon_sort)-resx/2, max(lon_sort)+resx/2,
              min(lat_sort)-resy/2, max(lat_sort)+resy/2),
    crs = "EPSG:4326"
  )
}

nc_read_block <- function(nc, var, diminfo,
                          row0_terra, nrb_terra, nlon,
                          lon_asc, lat_asc, fillvalue,
                          time_index) {
  
  nlat <- diminfo$lat_n
  lat_start <- if (lat_asc) (nlat - row0_terra - nrb_terra + 2) else row0_terra
  lat_count <- nrb_terra
  
  lon_start <- 1L
  lon_count <- nlon
  
  nd <- length(diminfo$dn)
  start <- rep(1L, nd)
  count <- vapply(nc$var[[var]]$dim, `[[`, 1L, "len")
  
  start[diminfo$lon_id]  <- lon_start;  count[diminfo$lon_id]  <- lon_count
  start[diminfo$lat_id]  <- lat_start;  count[diminfo$lat_id]  <- lat_count
  start[diminfo$time_id] <- time_index; count[diminfo$time_id] <- 1L
  
  arr <- ncvar_get(nc, var, start=start, count=count, collapse_degen=FALSE)
  perm <- c(diminfo$lat_id, diminfo$lon_id, diminfo$time_id)
  m <- aperm(arr, perm)[,,1]
  
  if (!is.na(fillvalue)) m[abs(m - fillvalue) < 1e-12] <- NA_real_
  if (lat_asc)  m <- m[nrow(m):1, , drop=FALSE]
  if (!lon_asc) m <- m[, ncol(m):1, drop=FALSE]
  m
}

on_read_fail <- function(action, var_tag, date_str, idx, nrb, nlon, msg) {
  if (action == "stop") stop(sprintf("[%s] Read failed at %s (time index %d): %s", var_tag, date_str, idx, msg))
  if (action == "skip_month") return(matrix(NA_real_, nrow=nrb, ncol=nlon))
  if (action == "set_zero")   return(matrix(0, nrow=nrb, ncol=nlon))
  stop("Unknown READ_FAIL_ACTION: ", action)
}

safe_copy <- function(src, dst) {
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(src, dst, overwrite = TRUE)
  if (!ok) stop("Copy failed: ", dst)
}

# ----------------------------
# Localize NetCDFs (avoid Y: read problems)
# ----------------------------
r_local  <- localize_nc(r_file)
pet_local <- if (MODE=="CWB1") localize_nc(pet_file) else NA_character_
ev_local  <- if (MODE=="CWB2") localize_nc(evspsbl_file) else NA_character_

# ----------------------------
# Open NetCDFs via ncdf4
# ----------------------------
ncP <- nc_open(r_local); on.exit(nc_close(ncP), add=TRUE)
varP <- "r"
dimP <- get_dim_info(ncP, varP)

if (MODE=="CWB1") {
  ncE <- nc_open(pet_local); on.exit(nc_close(ncE), add=TRUE)
  varE <- "pet"
} else {
  ncE <- nc_open(ev_local);  on.exit(nc_close(ncE), add=TRUE)
  varE <- "evspsbl"
}
dimE <- get_dim_info(ncE, varE)

lonP <- ncvar_get(ncP, dimP$lon_name); latP <- ncvar_get(ncP, dimP$lat_name)
lonE <- ncvar_get(ncE, dimE$lon_name); latE <- ncvar_get(ncE, dimE$lat_name)

if (length(lonP)!=length(lonE) || length(latP)!=length(latE) || any(lonP!=lonE) || any(latP!=latE)) {
  stop("P and EV grids differ. Streamed backend requires identical lon/lat vectors.")
}

lon <- lonP; lat <- latP
lon_asc <- lon[1] < lon[length(lon)]
lat_asc <- lat[1] < lat[length(lat)]

tP <- cf_time_to_month_start(ncP, as.Date("1950-01-01"))
tE <- cf_time_to_month_start(ncE, as.Date("1950-01-01"))
kP <- format(tP, "%Y-%m-01")
kE <- format(tE, "%Y-%m-01")

common <- sort(intersect(kP, kE))
yrs_common <- as.integer(substr(common, 1, 4))
common <- common[yrs_common >= PERIOD_START & yrs_common <= PERIOD_END]
if (length(common)==0) stop("No overlapping months in requested period.")

idxP_all <- match(common, kP)
idxE_all <- match(common, kE)
dates_common <- as.Date(common)

years <- PERIOD_START:PERIOD_END
years <- years[years %in% as.integer(format(dates_common, "%Y"))]
if (length(years)==0) stop("No overlapping years in requested period.")

uP <- get_units(ncP, varP); lnP <- get_long_name(ncP, varP)
uE <- get_units(ncE, varE); lnE <- get_long_name(ncE, varE)
message("P  long_name: ", lnP)
message("P  units    : ", uP)
message("EV long_name: ", lnE)
message("EV units    : ", uE)

multP <- monthly_multiplier(uP, lnP, dates_common)
multE <- monthly_multiplier(uE, lnE, dates_common)

fillP <- get_fill_value(ncP, varP)
fillE <- get_fill_value(ncE, varE)

tmpl <- make_template_from_lonlat(lon, lat)
nlon <- ncol(tmpl)

bs <- terra::blocks(tmpl)
row_vec  <- if (is.data.frame(bs)) bs$row   else bs$row
nrow_vec <- if (is.data.frame(bs)) bs$nrows else bs$nrows
n_blocks <- length(row_vec)

# ----------------------------
# Output dirs (NEW run folder)
# ----------------------------
tag <- if (MODE=="CWB1") "CWB1_r_minus_pet" else "CWB2_r_minus_evspsbl"
sub_tag <- sprintf("%d_%d_RUN_%s", PERIOD_START, PERIOD_END, RUN_ID)

local_out_dir <- file.path(LOCAL_OUT_BASE, tag, sub_tag)
dir.create(local_out_dir, recursive=TRUE, showWarnings=FALSE)

local_annP_dir <- file.path(local_out_dir, "annual_P_total_mm_yr")
local_annE_dir <- file.path(local_out_dir, "annual_EV_total_mm_yr")
local_annD_dir <- file.path(local_out_dir, "annual_CWB_mm_yr")
dir.create(local_annP_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(local_annE_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(local_annD_dir, recursive=TRUE, showWarnings=FALSE)

Y_out_dir <- if (COPY_TO_Y) file.path(out_base_Y, tag, sub_tag) else NA_character_

annual_CWB_local <- character(0)

# ----------------------------
# Main loop
# ----------------------------
for (y in years) {
  message("Year: ", y)
  pos <- which(as.integer(format(dates_common, "%Y")) == y)
  if (length(pos)==0) next
  n_months_y <- length(pos)
  
  fP <- file.path(local_annP_dir, sprintf("P_total_%d.tif", y))
  fE <- file.path(local_annE_dir, sprintf("EV_total_%d.tif", y))
  fD <- file.path(local_annD_dir, sprintf("CWB_%d.tif", y))
  
  # IMPORTANT: do not reuse/overwrite old names from previous failed run; these are new per RUN_ID anyway
  if (file.exists(fP) || file.exists(fE) || file.exists(fD)) {
    stop("Unexpected existing output in run folder: ", y)
  }
  
  # Start streaming writers (local disk)
  rP_out <- terra::writeStart(tmpl, filename=fP, overwrite=TRUE,
                              wopt=list(datatype="FLT4S", gdal=c("COMPRESS=LZW")))
  rE_out <- terra::writeStart(tmpl, filename=fE, overwrite=TRUE,
                              wopt=list(datatype="FLT4S", gdal=c("COMPRESS=LZW")))
  rD_out <- terra::writeStart(tmpl, filename=fD, overwrite=TRUE,
                              wopt=list(datatype="FLT4S", gdal=c("COMPRESS=LZW")))
  
  if (!inherits(rP_out, "SpatRaster") || !inherits(rE_out, "SpatRaster") || !inherits(rD_out, "SpatRaster")) {
    stop("terra::writeStart failed. This indicates an unreleased GDAL handle in the session. Restart R and rerun.")
  }
  
  on.exit({
    try(terra::writeStop(rP_out), silent=TRUE)
    try(terra::writeStop(rE_out), silent=TRUE)
    try(terra::writeStop(rD_out), silent=TRUE)
  }, add=TRUE)
  
  for (b in seq_len(n_blocks)) {
    row0 <- as.integer(row_vec[b])
    nrb  <- as.integer(nrow_vec[b])
    
    sumP <- matrix(0,  nrow=nrb, ncol=nlon)
    sumE <- matrix(0,  nrow=nrb, ncol=nlon)
    cntP <- matrix(0L, nrow=nrb, ncol=nlon)
    cntE <- matrix(0L, nrow=nrb, ncol=nlon)
    
    for (pp in pos) {
      date_str <- as.character(dates_common[pp])
      
      tiP <- idxP_all[pp]
      mP <- tryCatch(
        nc_read_block(ncP, varP, dimP, row0, nrb, nlon, lon_asc, lat_asc, fillP, tiP),
        error=function(e) on_read_fail(READ_FAIL_ACTION, "P", date_str, tiP, nrb, nlon, conditionMessage(e))
      )
      
      tiE <- idxE_all[pp]
      mE <- tryCatch(
        nc_read_block(ncE, varE, dimE, row0, nrb, nlon, lon_asc, lat_asc, fillE, tiE),
        error=function(e) on_read_fail(READ_FAIL_ACTION, "EV", date_str, tiE, nrb, nlon, conditionMessage(e))
      )
      
      mP <- mP * multP[pp]
      mE <- mE * multE[pp]
      if (MODE=="CWB2" && CLAMP_NEGATIVE_EV_TO_ZERO) mE[mE < 0] <- 0
      
      vP <- !is.na(mP); vE <- !is.na(mE)
      if (any(vP)) { tmp <- mP; tmp[!vP] <- 0; sumP <- sumP + tmp; cntP <- cntP + (vP*1L) }
      if (any(vE)) { tmp <- mE; tmp[!vE] <- 0; sumE <- sumE + tmp; cntE <- cntE + (vE*1L) }
    }
    
    P_ann <- sumP; E_ann <- sumE
    P_ann[cntP < n_months_y] <- NA_real_
    E_ann[cntE < n_months_y] <- NA_real_
    CWB <- P_ann - E_ann
    CWB[is.na(P_ann) | is.na(E_ann)] <- NA_real_
    
    start_cell <- terra::cellFromRowCol(tmpl, row=row0, col=1)
    rP_out <- terra::writeValues(rP_out, as.vector(t(P_ann)), start_cell)
    rE_out <- terra::writeValues(rE_out, as.vector(t(E_ann)), start_cell)
    rD_out <- terra::writeValues(rD_out, as.vector(t(CWB)),  start_cell)
  }
  
  rP_out <- terra::writeStop(rP_out)
  rE_out <- terra::writeStop(rE_out)
  rD_out <- terra::writeStop(rD_out)
  
  annual_CWB_local <- c(annual_CWB_local, fD)
  
  if (COPY_TO_Y) {
    dir.create(file.path(Y_out_dir, "annual_P_total_mm_yr"), recursive=TRUE, showWarnings=FALSE)
    dir.create(file.path(Y_out_dir, "annual_EV_total_mm_yr"), recursive=TRUE, showWarnings=FALSE)
    dir.create(file.path(Y_out_dir, "annual_CWB_mm_yr"), recursive=TRUE, showWarnings=FALSE)
    safe_copy(fP, file.path(Y_out_dir, "annual_P_total_mm_yr", basename(fP)))
    safe_copy(fE, file.path(Y_out_dir, "annual_EV_total_mm_yr", basename(fE)))
    safe_copy(fD, file.path(Y_out_dir, "annual_CWB_mm_yr", basename(fD)))
  }
}

if (length(annual_CWB_local)==0) stop("No annual CWB rasters created.")

# Multi-year stats (local)
ann_stack <- rast(annual_CWB_local)
cwb_mean <- mean(ann_stack, na.rm=TRUE)
cwb_min  <- min(ann_stack,  na.rm=TRUE)
cwb_max  <- max(ann_stack,  na.rm=TRUE)

out_mean <- file.path(local_out_dir, sprintf("%s_CWB_mean_%d_%d_mm_yr.tif", tag, PERIOD_START, PERIOD_END))
out_min  <- file.path(local_out_dir, sprintf("%s_CWB_min_%d_%d_mm_yr.tif",  tag, PERIOD_START, PERIOD_END))
out_max  <- file.path(local_out_dir, sprintf("%s_CWB_max_%d_%d_mm_yr.tif",  tag, PERIOD_START, PERIOD_END))

writeRaster(cwb_mean, out_mean, overwrite=TRUE, gdal=c("COMPRESS=LZW"))
writeRaster(cwb_min,  out_min,  overwrite=TRUE, gdal=c("COMPRESS=LZW"))
writeRaster(cwb_max,  out_max,  overwrite=TRUE, gdal=c("COMPRESS=LZW"))

if (COPY_TO_Y) {
  safe_copy(out_mean, file.path(Y_out_dir, basename(out_mean)))
  safe_copy(out_min,  file.path(Y_out_dir, basename(out_min)))
  safe_copy(out_max,  file.path(Y_out_dir, basename(out_max)))
}

message("DONE")
message("Local output folder: ", local_out_dir)
if (COPY_TO_Y) message("Y: output folder: ", Y_out_dir)
message("Local NetCDF cache: ", LOCAL_NC_CACHE)
message("Terra tempdir: ", LOCAL_TERRA_TMP)
