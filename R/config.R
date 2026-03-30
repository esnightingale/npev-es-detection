################################################################################
# Analysis configuration - set country and time period before running workflow
# Source this file at the start of each script, or run interactively to set params
################################################################################

# Country to analyse: "Pakistan", "Madagascar" or "Nigeria"
country <- "Nigeria"

# Time period (inclusive)
tp <- lubridate::ymd(c("2021-01-01", "2024-12-01"))
years <- lubridate::year(tp)
years_lab <- paste0(years, collapse = "-")

# Derived parameters (do not edit unless adding new country)
iso_code <- switch(country,
  "Pakistan" = "PAK",
  "Madagascar" = "MDG",
  "Nigeria" = "NGA",
  stop("Unknown country. Add mapping in R/config.R")
)

# Local projection for spatial operations (UTM)
proj_local <- switch(country,
  "Pakistan" = "EPSG:32642",  # UTM 42N
  "Madagascar" = "EPSG:32738", # UTM 38S
  "Nigeria" = "EPSG:32632",   # UTM 32N
  stop("Unknown country. Add mapping in R/config.R")
)

# WorldPop raster prefix (pak, nga, etc.)
raster_prefix <- tolower(iso_code)

# Paths
dir <- file.path("data", country)
# Dropbox location for raw POLIS exports (manually download per country)
dropbox_polis <- "/Users/phpuenig/Dropbox/Polio-SPEC/Analysis/Data/polis/csv"

################################################################################
