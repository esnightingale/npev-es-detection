################################################################################
# Read and tidy raw POLIS downloads for AFP cases and ES samples
# + Select / setup vars of interest
# + Impute missing coordinates with district centroids
# + Define additional variables (classify sites by sampling regularity)
# + No subsetting by time/country
################################################################################

library(logr)
library(here)
library(parsedate)
library(tidyverse)
library(sf)

# In/out directories
subdir <- "Pakistan"
datadir <- paste0("/Users/phpuenig/Dropbox/Polio-SPEC/Analysis/Data/polis/csv/",subdir) 
outdir <- "data"

# Choose to include or exclude AFP cases with missing district
impute_miss_adm <- FALSE 

log_open(here(outdir,subdir,"./Cleaning/clean_polis.log"))

log_print(paste("POLIS data cleaning log for",subdir))

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
# District shapefiles

log_print("Reading district shapefiles")

source("../polio-spec-models/R/utils/IDM/polis_shapes_query_sf.R")
filename_shapes = paste0("/Users/phpuenig/Dropbox/Polio-SPEC/Analysis/Data/",
                         "WHO_POLIO_GLOBAL_GEODATABASE.gdb")

shape2 <- st_read_geodb(filename_shapes,2) %>% 
  st_make_valid() |> 
  filter(enddate > "2024-01-01", iso_3_code == "PAK") 
shape1 <- read_sf_geodb(filename_shapes,1) |> 
  filter(enddate > "2024-01-01", iso_3_code == "PAK")
shape0 <- read_sf_geodb(filename_shapes,0) |> 
  filter(enddate > "2024-01-01", iso_3_code == "PAK")

#---------------------------------------------------------------------------- #
# NPAFP rate (for population denominator)

log_print("Reading NPAFP data")

periods = date_to_period(seq(ymd("2020-01-01"),ymd("2024-12-01"), by = "month"))

npafp_raw <- readxl::read_xlsx(file.path(datadir, 
                                       list.files(datadir,pattern = "Indicators"))) |> 
  clean_names() %>%
  mutate(guid = paste0("{",toupper(admin_2_guid),"}"),
         period = date_to_period(dmy(paste0("1",period))))

row_fill <- npafp_raw |>
  expand(guid = unique(shape2$guid), period = periods) 
  
# Expand time points for each guid
npafp <- npafp_raw |>
  right_join(row_fill) |>
  mutate(month = ymd(paste0(period_to_date(period), "-01")),
         year = floor(period),
         missing_denom = is.na(denominator)) |> 
  select(guid, period, month, year, numerator, denominator, missing_denom) 

names(npafp)
summary(npafp)

tabyl(npafp, missing_denom, year)
tabyl(npafp, year)
# All missing for 2021?

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

afp_cols <- names(readxl::read_excel(file.path(datadir, 
                                              list.files(datadir,pattern = "Human")), n_max = 0))
afp_types <- rep("text",length(afp_cols))
afp_types[grepl("Date",afp_cols)] <- "date"
afp_types[grepl("Age",afp_cols) | grepl("Number",afp_cols)] <- "numeric"
afp_types[which(afp_cols %in% c("X","Y"))] <- "numeric"
afp_types[grepl("PoNS",afp_cols) | grepl("Exact",afp_cols) | grepl("Paralysis",afp_cols)| grepl("Contact",afp_cols)] <- "skip"

afp_raw <- readxl::read_xlsx(file.path(datadir, 
                                       list.files(datadir, pattern = "Human")),
                             col_types = afp_types) %>%
  clean_names()

names(afp_raw)
summary(afp_raw)

tabyl(afp_raw$classification, exclude = NULL)

log_print("Set up AFP variables and drop 'Not an AFP' cases")

afp_raw %>% 
  rename(iso_3_code = country_iso3,
         # Case date is defined as date of paralysis onset - take this as general date for FFI analysis
         date = case_date,
         onset_date = date_onset,
         report_wk = afp_reporting_week) %>% 
  mutate(guid = paste0("{",toupper(admin_2_guid),"}"),
         period = date_to_period(date),
         month = floor_date(date, "month"),
         year = year(date),
         final_class = case_when(grepl("WILD1",virus_type_s) ~ "WPV1",
                           grepl("WILD3",virus_type_s) ~ "WPV3",
                           grepl("cVDPV2", virus_type_s) ~ "cVDPV2",
                           grepl("aVDPV2", virus_type_s) ~ "aVDPV2",
                           grepl("VDPV", virus_type_s) ~ "VDPV",
                           grepl("NPEV", virus_type_s) ~ "NPEV"),
         # Calculated age has the least missingness out of all age vars
         age_m = calculated_age_months) %>% 
  select(epid, place_admin_1, place_admin_2, guid:year, date,
         onset_date, notification_date, investigation_date, stool_date_sent_to_lab, report_wk,
         age_m,
         final_class, classification,
         vaccine_1:npev, stool_adequacy,
         virus_type_s, virus_cluster_s, virus_is_orphan,emergence_group_s,
         iso_3_code,x,y) %>% 
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

write_csv(afp_miss_guid, here(outdir,subdir,"Cleaning/afp_missing_adm2.csv"))

# => get guid from shapefile in which point falls if coordinates avaliable?

afp_miss_guid_sf <- filter(afp, guid == "{NA}") %>% 
  st_as_sf(coords = c("x","y"), crs = st_crs(4326))

# Plot to check
ggplot() +
  geom_sf(data = filter(shape2, iso_3_code == "PAK")) +
  geom_sf(data= filter(afp_miss_guid_sf, iso_3_code == "PAK"), cex = 0.7) + 
  coord_sf()
ggsave(here(outdir,subdir,"Cleaning/afp_missing_adm2.png"), height = 7, width = 8)

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
es_raw %>% 
  rename(iso_3_code = country_iso3,
         sample_id = env_sample_id,
         x = y,
         y = x) %>% 
  separate(environmental_site, into = c("site_id","site_name"), sep = " -") %>% 
  mutate(guid = paste0("{",toupper(admin_2_guid),"}"),
         # across(c(collection_date, date_received_in_lab),
         #        lubridate::dmy),
         month = floor_date(ymd(collection_date), "month"),
         year = year(collection_date),
         final_class = case_when(grepl("WILD1",virus_type_s) ~ "WPV1",
                           grepl("WILD3",virus_type_s) ~ "WPV3",
                           grepl("cVDPV2", virus_type_s) ~ "cVDPV2",
                           grepl("aVDPV2", virus_type_s) ~ "aVDPV2",
                           grepl("VDPV", virus_type_s) ~ "VDPV",
                           grepl("VACCINE", virus_type_s) ~ "Vaccine",
                           grepl("NPEV", virus_type_s) ~ "NPEV"),
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

log_print(paste("Saving cleaned datasets in", here(outdir,subdir)))

saveRDS(shape2, here(outdir,subdir,"shape2.rds"))
saveRDS(shape1, here(outdir,subdir,"shape1.rds"))
saveRDS(shape0, here(outdir,subdir,"shape0.rds"))

write_csv(afp, here(outdir,subdir,"linelist_afp_clean.csv"))
write_csv(es, here(outdir,subdir,"linelist_es_clean.csv")) 
write_csv(npafp, here(outdir,subdir,"npafp_clean.csv"))

################################################################################
################################################################################