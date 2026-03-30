################################################################################
# Read and tidy raw POLIS downloads for AFP cases and ES samples
# + Select / setup vars of interest
# + virus_type_ord_1..5 and virus_type_n_unique from virus_type_s (AFP + ES)
# + Impute missing coordinates with district centroids
# + Define additional variables (classify sites by sampling regularity)
# + No subsetting by time/country
################################################################################

library(logr)
library(here)
library(parsedate)
library(tidyverse)
library(sf)

source(here("R/config.R"))

# In/out directories
datadir <- file.path(dropbox_polis, country)
outdir <- here("data",country)

# Choose to include or exclude AFP cases with missing district
impute_miss_adm <- FALSE 

log_open(here(outdir, country, "./Cleaning/clean_polis.log"))

log_print(paste("POLIS data cleaning log for", country))

# ---------------------------------------------------------------------------- #
# Additional helper functions

period_to_date = function(x,format='%Y-%b'){
  yr = floor(x)
  mo = round(12*(x-yr),0)+1
  d = '1'
  format(ymd(paste(yr,mo,d,sep='-')),format)
}

date_to_period = function(x){
  year(x) + (month(x)-1)/12
}

# ---------------------------------------------------------------------------- #
# virus_type_s: five logical indicators (fixed POLIS reporting order / groupings)
# 1 NPEV, 2 Sabin/vaccine (VACCINE*), 3 VDPV (cVDPV / aVDPV / VDPV*), 4 wild (WILD*),
# 5 other non-empty tokens not in 1–4. Missing/empty virus_type_s -> all FALSE.
# Plus virus_type_n_unique = count of distinct trimmed tokens (0 if missing/empty).

parse_virus_type_tokens <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(character(0))
  }
  x <- x[1]
  if (is.na(x)) {
    return(character(0))
  }
  x <- str_trim(as.character(x))
  if (!nzchar(x)) {
    return(character(0))
  }
  parts <- unlist(strsplit(x, ",", fixed = TRUE))
  parts <- str_trim(parts)
  parts[nzchar(parts)]
}

virus_type_token_class <- function(t) {
  t <- str_trim(t)
  if (!nzchar(t) || is.na(t)) {
    return(NA_integer_)
  }
  if (t == "NPEV") {
    return(1L)
  }
  if (grepl("WILD", t, fixed = TRUE)) {
    return(4L)
  }
  if (grepl("VACCINE", t, fixed = TRUE)) {
    return(2L)
  }
  if (grepl("cVDPV", t, fixed = TRUE) || grepl("aVDPV", t, fixed = TRUE) ||
      grepl("VDPV", t, fixed = TRUE)) {
    return(3L)
  }
  5L
}

add_virus_type_order_indicators <- function(d) {
  if (!"virus_type_s" %in% names(d)) {
    return(d)
  }
  tok_lists <- lapply(seq_len(nrow(d)), function(i) parse_virus_type_tokens(d$virus_type_s[i]))
  cls_mat <- matrix(FALSE, nrow = nrow(d), ncol = 5L)
  for (i in seq_along(tok_lists)) {
    toks <- tok_lists[[i]]
    if (!length(toks)) {
      next
    }
    for (t in toks) {
      k <- virus_type_token_class(t)
      if (!is.na(k) && k >= 1L && k <= 5L) {
        cls_mat[i, k] <- TRUE
      }
    }
  }
  d$virus_type_ord_1 <- cls_mat[, 1L]
  d$virus_type_ord_2 <- cls_mat[, 2L]
  d$virus_type_ord_3 <- cls_mat[, 3L]
  d$virus_type_ord_4 <- cls_mat[, 4L]
  d$virus_type_ord_5 <- cls_mat[, 5L]
  d$virus_type_n_unique <- vapply(tok_lists, function(toks) length(unique(toks)), integer(1))
  d
}

# ---------------------------------------------------------------------------- #
# District shapefiles

log_print("Reading district shapefiles")

shape2 <- readRDS(here(outdir,"shape2.rds"))

#---------------------------------------------------------------------------- #
# NPAFP rate (for population denominator)

log_print("Reading NPAFP data")

periods = date_to_period(seq(tp[1],tp[2], by = "month"))

npafp_raw <- readxl::read_xlsx(file.path(datadir, 
                                       list.files(datadir,pattern = "Indicators"))) |> 
  clean_names() %>%
  mutate(guid = paste0("{",toupper(admin_2_guid),"}"),
         period = date_to_period(dmy(paste0("1",period))))

row_fill <- npafp_raw |>
  expand(guid = unique(shape2$guid_old), period = periods) 
  
# Expand time points for each guid
npafp <- npafp_raw |>
  right_join(row_fill) |>
  mutate(month = ymd(paste0(period_to_date(period), "-01")),
         year = floor(period),
         polis_npafp_rate = value,
         missing_denom = is.na(denominator)) |> 
  select(guid, period, month, year, polis_npafp_rate, numerator, denominator, missing_denom) 

names(npafp)
summary(npafp)

tabyl(npafp, missing_denom, year)
tabyl(npafp, year)

log_print(paste0("Imputing ",
                 sum (npafp$missing_denom)," (",round(mean(npafp$missing_denom),2)*100,
                 "%) missing denominators, for ",
                 unique(npafp$guid[npafp$missing_denom]) %>% length, " unique districts."))

# Impute with nearest non-missing value in time
npafp <- npafp |> 
  group_by(guid) |> 
  fill(denominator, .direction = "downup") |> 
  ungroup() |> 
  rename(imputed_denom = missing_denom)

# shape2 <- left_join(shape2, select(npafp, guid, pop_u15_polisdenom))

#---------------------------------------------------------------------------- #
# AFP linelist

log_print("Reading AFP data")

# Prefer AFP case file over Specimen file when both exist
human_files <- list.files(datadir, pattern = "Human")
afp_file <- human_files[grepl("AFP_Only", human_files) |
                         (grepl("Detailed", human_files) & !grepl("Specimen", human_files))]
if (length(afp_file) == 0) afp_file <- human_files[1]

afp_cols <- names(readxl::read_excel(file.path(datadir, afp_file), n_max = 0))
afp_types <- rep("text", length(afp_cols))
afp_types[grepl("Date", afp_cols)] <- "date"
afp_types[grepl("Age", afp_cols) | grepl("Number", afp_cols)] <- "numeric"
afp_types[which(afp_cols %in% c("X", "Y"))] <- "numeric"
afp_types[grepl("PoNS", afp_cols) | grepl("Exact", afp_cols) | grepl("Paralysis", afp_cols) | grepl("Contact", afp_cols)] <- "skip"

afp_raw <- readxl::read_xlsx(file.path(datadir, afp_file), col_types = afp_types) %>%
  clean_names()

names(afp_raw)
summary(afp_raw)

log_print("Set up AFP variables and drop 'Not an AFP' cases")

afp_raw %>%
  rename(iso_3_code = country_iso3,
         date = case_date,
         onset_date = date_onset,
         report_wk = afp_reporting_week) %>%
  mutate(guid = paste0("{", toupper(admin_2_guid), "}"),
         period = date_to_period(date),
         month = floor_date(date, "month"),
         year = year(date),
         age_m = calculated_age_months) %>%
  select(epid, place_admin_1, place_admin_2, guid, period, month, year, date,
         onset_date, notification_date, investigation_date, stool_date_sent_to_lab, report_wk,
         age_m, classification,
         vaccine_1:npev, stool_adequacy,
         virus_type_s, virus_cluster_s, virus_is_orphan, emergence_group_s,
         iso_3_code, x, y) %>%
  filter(classification != "Not an AFP") -> afp

log_print("Checking ages")

# Check age distribution
summary(afp$age_m)

# Remove erroneous age > 120
afp$age_m[afp$age_m > 120] <- NA

log_print("Checking missing coordinates")

# Check and impute missingness in coordinates if necessary
summary(select(afp, x, y))

if (any(is.na(afp$x), is.na(afp$y))){
  afp %>%
    mutate(coord_imp = (is.na(x) | is.na(y))) %>%
    left_join(select(shape2, guid, center_lat, center_lon) %>%
                st_drop_geometry()) %>%
    mutate(x = coalesce(x, center_lon),
           y = coalesce(y, center_lat)) -> afp
  
  summary(select(afp, x, y))
  summary(afp$coord_imp)
}

log_print("Checking missingness in admin_2_guid")

table(afp$guid == "{NA}")

afp_miss_guid <- filter(afp, guid == "{NA}")
afp_miss_guid %>% 
  select(epid, place_admin_1, place_admin_2, guid, x, y) 

write_csv(afp_miss_guid, here(outdir, country, "Cleaning/afp_missing_adm2.csv"))

# => get guid from shapefile in which point falls if coordinates avaliable?

afp_miss_guid_sf <- filter(afp, guid == "{NA}") %>% 
  st_as_sf(coords = c("x","y"), crs = st_crs(4326))

# Plot to check
ggplot() +
  geom_sf(data = filter(shape2, iso_3_code == iso_code)) +
  geom_sf(data = filter(afp_miss_guid_sf, iso_3_code == iso_code), cex = 0.7) + 
  coord_sf()
ggsave(here(outdir, country, "Cleaning/afp_missing_adm2.png"), height = 7, width = 8)

#### ----- INCLUDE OR EXCLUDE CASES WITH MISSING GUID ------ ####

log_print(paste("Impute missing district:",impute_miss_adm))

if (impute_miss_adm){
  
  log_print(paste("INCLUDE and IMPUTE missing admin_2_guid for", nrow(afp_miss_guid), "cases"))

    # Fill in guid according to district in which coordinates fall
  afp %>% 
    mutate(guid_imp = guid == "{NA}") %>%
    st_as_sf(coords = c("x","y"), crs = st_crs(4326), remove = FALSE) %>% 
    st_join(select(shape2, guid), join = st_within, left = TRUE) %>%
    st_drop_geometry() %>% 
    group_by(epid) %>%
    # Remove duplicated where point falls across multiple districts
    slice_head(n = 1) %>%
    mutate(guid.x = na_if(guid.x, "{NA}"),
           guid = coalesce(guid.x, guid.y)) %>% 
    ungroup() %>% 
    select(-guid.x, -guid.y) -> tmp
  
  # Check uniqueness of original records
  tmp %>% group_by(epid) %>% tally() 
  
  table(is.na(tmp$guid),tmp$guid_imp)
  
  # Set this as main dataset
  afp <- tmp
  
}else{
  log_print(paste("EXCLUDE", nrow(afp_miss_guid), "cases with missing admin_2_guid"))
  afp <- filter(afp, guid != "{NA}")
}

# ---------------------------------------------------------------------------- #
# ES linelist

log_print("Reading ES data")

es_cols <- names(readxl::read_excel(file.path(datadir, 
                                      list.files(datadir,pattern = "Env")), n_max = 0))
es_types <- rep("text",length(es_cols))
es_types[grepl("Date",es_cols)] <- "date"
es_types[which(es_cols %in% c("X","Y"))] <- "numeric"

es_raw <- readxl::read_xlsx(file.path(datadir, 
                                      list.files(datadir,pattern = "Env")),
                            col_types = es_types) %>%
  clean_names()

log_print("Set up ES variables")
# Nigeria ES has both env_sample_id and sample_id; Pakistan has env_sample_id only
es_raw <- es_raw %>%
  rename(iso_3_code = country_iso3,
         x = y,
         y = x)
if ("sample_id" %in% names(es_raw) && "env_sample_id" %in% names(es_raw)) {
  es_raw <- select(es_raw, -env_sample_id)  # keep sample_id, drop duplicate
} else if ("env_sample_id" %in% names(es_raw)) {
  es_raw <- rename(es_raw, sample_id = env_sample_id)
}
es_raw <- es_raw %>%
  separate(environmental_site, into = c("site_id","site_name"), sep = " -") %>% 
  mutate(guid = paste0("{",toupper(admin_2_guid),"}"),
         # across(c(collection_date, date_received_in_lab),
         #        lubridate::dmy),
         month = floor_date(ymd(collection_date), "month"),
         year = year(collection_date),
         travel_time = difftime(date_received_in_lab, collection_date, "days"),
        across(site_id:site_name, trimws))  -> es

log_print("Checking missing coordinates")
summary(select(es, x, y))
table(es$guid == "{NA}")

log_print("Impute missing site coordinates with district centroid")
es %>% 
  mutate(coord_imp = (is.na(x) | is.na(y))) %>% 
  left_join(select(shape2, guid, center_lat, center_lon) %>% 
              st_drop_geometry()) %>%
  mutate(x = coalesce(x, center_lon),
         y = coalesce(y, center_lat)) -> es

summary(select(es, x, y))
tabyl(es$coord_imp)

log_print("Check uniqueness of locations per site id")
es %>% 
  group_by(site_id) %>% 
  summarise(x = n_distinct(x), y = n_distinct(y)) %>%
  summary()

# ---------------------------------------------------------------------------- #
# Save cleaned datasets

log_print("Derive virus_type_ord_1..5 and virus_type_n_unique from virus_type_s")
afp <- add_virus_type_order_indicators(afp)
es <- add_virus_type_order_indicators(es)

log_print(paste("Saving cleaned datasets in", outdir))

write_csv(afp, here(outdir, "linelist_afp_clean.csv"))
write_csv(es, here(outdir, "linelist_es_clean.csv")) 
write_csv(npafp, here(outdir, "npafp_clean.csv"))

################################################################################
################################################################################