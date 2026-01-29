# Copernicus Interactive Climate Atlas (C3S) — Climate Statistics from NetCDF (R / terra)

This repository contains **R scripts** to compute climate statistics from **monthly and annual NetCDF** products downloaded from the **Copernicus Climate Data Store (CDS) — Multi-origin C3S Atlas** and stored locally.

Data source (download portal):  
https://cds.climate.copernicus.eu/datasets/multi-origin-c3s-atlas?tab=download

All scripts use **terra** (and **ncdf4**) and are designed for large rasters (e.g., ERA5-Land global grids). For memory-heavy workflows, the scripts implement **year-by-year processing** (12 months at a time) to prevent R from terminating due to RAM exhaustion.

---

## Requirements

### R packages
- `terra`
- `ncdf4`

Install:
```r
install.packages(c("terra", "ncdf4"))
```

### Notes on CRS warnings
When reading some NetCDF files, GDAL may warn about axis units. For ERA5-Land lat/lon grids, scripts set:
```r
crs(r) <- "EPSG:4326"
```
to ensure GeoTIFF output has correct CRS metadata. The warning is not fatal for these time aggregations.

---

## Local data layout (example)

You keep datasets in:
```
Y:\Maps\climate\CDS_CopernicusInteractiveClimateAtlas\
```

Example subfolders used in the scripts:
- `MthVeryHeavyPrecipitationDays`
- `AnnualMaximumConsecutiveDryDays`
- `!NewMthMax1-dayPrecipitation`
- `!NewMthPrecipitation`
- `MthDailyAccumulatedPotentialEvapotranspiration`
- `!newMthEvapotransiration`

Each script writes outputs into a `stats/` folder within the corresponding dataset folder.

---

## Implemented statistics

### 1) Monthly very heavy precipitation days (R20 / `r20`)
**Definition:** Monthly count of days with daily accumulated precipitation (liquid water equivalent) **> 20 mm** (units: days).  
**Statistic computed:**
1. For each year: **MAX** across the 12 monthly values  
2. For **1995–2024**: **MEAN** of those annual maxima

**Input example**
```
...\MthVeryHeavyPrecipitationDays\r20_ERA5-Land_no-expt_mon_19500101-20241201.nc
var = r20
```

**Output**
- `stats/r20_mean_annual_max_monthly_1995_2024.tif`

---

### 2) Annual maximum consecutive dry days (CDD / `cdd`)
**Definition:** Annual maximum number of consecutive dry days (units: days).  
**Statistic computed:** Mean over **1995–2024**.

**Input example**
```
...\AnnualMaximumConsecutiveDryDays\cdd_ERA5-Land_no-expt_yr_19500101-20240101.nc
var = cdd
```

**Output**
- `stats/cdd_mean_1995_2024.tif`

---

### 3) Monthly maximum 1-day precipitation (RX1day / `rx1day`)
**Definition:** Monthly maximum 1-day precipitation amount (often mm).  
**Statistic computed:**
1. For each year: **MAX** across 12 monthly values  
2. For **1995–2024**: **MEAN** of those annual maxima

**Input example**
```
...\!NewMthMax1-dayPrecipitation\rx1day_ERA5-Land_no-expt_mon_19500101-20241201.nc
var = rx1day
```

**Output**
- `stats/rx1day_mean_annual_max_monthly_1995_2024.tif`

---

### 4) CWB1 = Monthly precipitation − Monthly daily accumulated potential evapotranspiration
**Definition:**
```
CWB1_month = r_month − pet_month
```
**Statistic computed:**
1. For each year: **MEAN** across 12 monthly CWB1 values  
2. For **1994–2024**: **MEAN** of those annual means

**Inputs**
```
r   : ...\!NewMthPrecipitation\r_ERA5-Land_no-expt_mon_19500101-20241201.nc
pet : ...\MthDailyAccumulatedPotentialEvapotranspiration\pet_ERA5-Land_no-expt_mon_19500101-20241201.nc
```

**Output**
- `stats/CWB1_mean_annual_mean_1994_2024.tif`

**Memory note:** This computation can be RAM-heavy. The repository includes a **memory-safe** script that processes **one year at a time** and writes annual rasters to disk, then averages them.

---

### 5) CWB2 = Monthly precipitation − Monthly evapotranspiration
**Definition:**
```
CWB2_month = r_month − evspsbl_month
```
**Statistic computed:**
1. For each year: **MEAN** across 12 monthly CWB2 values  
2. For **1994–2024**: **MEAN** of those annual means

**Inputs**
```
r       : ...\!NewMthPrecipitation\r_ERA5-Land_no-expt_mon_19500101-20241201.nc
evspsbl : ...\!newMthEvapotransiration\evspsbl_ERA5-Land_no-expt_mon_19500101-20241201.nc
```

**Output**
- `stats/CWB2_mean_annual_mean_1994_2024.tif`

**Units/sign note:** Some evap variables may be stored as rates (e.g., `kg m-2 s-1`) or use negative sign conventions. The memory-safe script:
- optionally converts rates to monthly totals using month length
- optionally flips sign if evap is negative on average in the first processed layer

---

## How to run

If you keep scripts under `scripts/`, run from R/RStudio:

```r
source("scripts/01_r20_mean_annual_max_monthly_1995_2024.R")
source("scripts/02_cdd_mean_1995_2024.R")
source("scripts/03_rx1day_mean_annual_max_monthly_1995_2024.R")
source("scripts/04_CWB1_mean_annual_mean_1994_2024_memory_safe.R")
source("scripts/05_CWB2_mean_annual_mean_1994_2024_memory_safe.R")
```

Outputs will be written to the relevant `stats/` folders.

---

## Performance / stability tips

1. **Use a fast temp directory with space** (important for large grids)
   ```r
   terraOptions(tempdir = "<path>", memfrac = 0.6, progress = 1)
   ```
2. Prefer the **memory-safe year-by-year** scripts for CWB1/CWB2.
3. If grids differ between files, scripts resample the second dataset to match precipitation.

---

## License
MIT
