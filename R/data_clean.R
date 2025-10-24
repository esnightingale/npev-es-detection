################################################################################
# Read and tidy raw POLIS downloads for AFP cases and ES samples
# + Select / setup vars of interest
# + Impute missing coordinates with district centroids
# + Define additional variables (classify sites by sampling regularity)
# + No subsetting by time/country
################################################################################

library(here)
library(parsedate)
library(tidyverse)

# In/out directories
country <- "Pakistan"
datadir <- paste0("/Users/phpuenig/Dropbox/Polio-SPEC/Analysis/Data/polis/csv",country) 
outdir <- "inputs"

source(here("R/utils/polis_shapes_query_sf.R"))

#---------------------------------------------------------------------------- #
# Surveillance indicators

country <- "Pakistan"
ind_afp <- read_csv(list.files(datadir,
                          "polis/csv/","Indicators_Detailed_Dataset_may_contain_sensitive_data_12-04-2024_13-08-29.csv")) %>%
  clean_names() %>%
  mutate(pop_u15_polisdenom = denominator,
         guid = paste0("{",admin_2_guid,"}"))

# View(ind_afp %>%
#        filter(indicator_id == 3527))

summary(ind_afp$denominator)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 3773   54379  101295  163528  154966 5311744    7528

shape2 <- left_join(shape2, select(ind_afp, guid, pop_u15_polisdenom))
#---------------------------------------------------------------------------- #
# AFP linelist

# Read data for Afghanistan, Pakistan and Nigeria
path <- paste0(datadir,"polis/csv/")
files <- list.files(path, pattern = "Human")
afp_raw <- lapply(files, 
                  function(f) read_csv(paste0(path,f))) %>% 
  bind_rows() %>% 
  clean_names() 

# afp_raw <- read_csv(paste0(datadir,
#                            "cvd2_master_2023-05-12-030357_e44c0c9_results/results/linelist_afp_raw.csv"))

table(afp_raw$classification, exclude = NULL)
# Compatible Confirmed (wild)        Discarded       Not an AFP          Pending             VDPV 
#        284              575           256058             1409             2557             1154 
table(is.na(afp_raw$admin_2_guid))
#  FALSE   TRUE 
# 261351    686 

afp_raw_miss_guid <- filter(afp_raw, is.na(admin_2_guid))
write_csv(afp_raw_miss_guid, here(path, "afp_missing_adm2.csv"))

afp_raw %>% 
  rename(iso_3_code = country_iso3,
         # y = exact_longitude,
         # x = exact_latitude
         # Case date is defined as date of paralysis onset - take this as general date for FFI analysis
         date = case_date,
         onset_date = date_onset,
         report_wk = afp_reporting_week) %>% 
  mutate(guid = paste0("{",toupper(admin_2_guid),"}"),
         across(c(onset_date, date, notification_date, investigation_date, stool_date_sent_to_lab),
                lubridate::dmy),
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
  select(epid, guid:month, date,
         onset_date, notification_date, investigation_date, stool_date_sent_to_lab, report_wk,
         age_m,
         final_class, classification,
         vaccine_1:npev, stool_adequacy,
         virus_type_s, virus_cluster_s, virus_is_orphan,emergence_group_s,
         iso_3_code,x,y
         # po_ns_id:po_ns_wild_cluster_name
         ) %>% 
  filter(classification != "Not an AFP") -> afp

# Check age
summary(afp$age_m)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# 0.00    19.00    36.00    47.66    62.00 22841.00      422 

# Remove erroneous age > 1000y
afp$age_m[afp$age_m > 100*12] <- NA

# Check missingness in guid
table(afp$guid == "{NA}",is.na(afp$x))
#        FALSE
# FALSE 261351
# TRUE     686
# => all records with missing guid have known coordinates
#    get guid from shapefile in which point falls

afp_missguid <- filter(afp, guid == "{NA}") %>% 
  st_as_sf(coords = c("x","y"), crs = st_crs(4326))

# Plot for Pakistan to check
ggplot() +
  geom_sf(data = filter(shape2, iso_3_code == "PAK")) +
  geom_sf(data= filter(afp_missguid, iso_3_code == "PAK"), cex = 0.7) + 
  coord_sf()

# Fill in guid according to district in which coordinates fall
afp %>% 
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
# tmp %>% group_by(epid) %>% tally() %>% View()

table(is.na(tmp$guid))

# Set this as main dataset
# afp <- tmp

# Check and impute missingness in coordinates if necessary
# summary(select(afp, x, y))

# afp %>%  
#   mutate(coord_imp = (is.na(x) | is.na(y))) %>% 
#   left_join(select(shape2, guid, center_lat, center_lon) %>% 
#               st_drop_geometry()) %>%
#   mutate(x = coalesce(x, center_lon),
#          y = coalesce(y, center_lat)) -> afp

# summary(select(afp, x, y))

# Exclude cases with missing admin 3
afp <- filter(afp, guid != "{NA}")

# ---------------------------------------------------------------------------- #
# ES linelist

# Read data for Afghanistan, Pakistan and Nigeria
es_raw <- read_csv(paste0(datadir,
                          "polis/csv/EnvSamples_Summary_Dataset_may_contain_sensitive_data_20-10-2023_12-02-32.csv"
                          # "EnvSamples_Detailed_Dataset_may_contain_sensitive_data_04-08-2023_14-32-23.csv"
                          )) %>% 
  clean_names()

es_raw %>% 
  rename(iso_3_code = country_iso3,
         sample_id = env_sample_id,
         x = y,
         y = x) %>% 
  separate(environmental_site, into = c("site_id","site_name"), sep = " -") %>% 
  mutate(guid = paste0("{",toupper(admin_2_guid),"}"),
         across(c(collection_date, date_received_in_lab),
                lubridate::dmy),
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

# Check and impute missingness in coords
summary(select(es, x, y))
# x                y         
# Min.   : 3.185   Min.   : 4.723  
# 1st Qu.: 7.461   1st Qu.:10.217  
# Median :12.455   Median :12.778  
# Mean   :34.283   Mean   :18.945  
# 3rd Qu.:68.766   3rd Qu.:31.456  
# Max.   :74.895   Max.   :36.743  
# NA's   :930      NA's   :930    

table(is.na(es$guid))
# FALSE 
# 25609 

es %>% 
  mutate(coord_imp = (is.na(x) | is.na(y))) %>% 
  left_join(select(shape2, guid, center_lat, center_lon) %>% 
              st_drop_geometry()) %>%
  mutate(x = coalesce(x, center_lon),
         y = coalesce(y, center_lat)) -> es

summary(select(es, x, y))

# Check uniqueness of locations per site id
# es %>% group_by(site_id) %>% summarise(x = n_distinct(x), y = n_distinct(y)) %>% View()

# ---------------------------------------------------------------------------- #
# Define consistency class

per <- range(es$month)
st <- per[1]
tdy <- ceiling_date(per[2],"month")

# Summarise by month of collection
ids <- unique(es$site_id)
months <- seq(st, tdy, by='months')
months <- data.frame(site_id = rep(ids, each = length(months)),
                     month = rep(months, times = length(ids)))

es %>% 
  group_by(site_id, month) %>% 
  count() %>%
  ungroup() %>% 
  right_join(months) %>%
  mutate(year = year(month),
         n = replace_na(n, 0)) -> es_mth

# By year of collection
es_mth %>%   
  group_by(site_id, year) %>% 
  summarise(monthly = all(n > 0),
            avg_monthly = sum(n) >= n_distinct(month),
            yr_total_samples = sum(n)) -> es_yr

es_yr %>% 
  # select(-yr_total_samples, -samples_ge12) %>% 
  dplyr::select(-yr_total_samples, -monthly) %>% 
  pivot_wider(names_from = year, 
              values_from = avg_monthly, 
              names_prefix = "monthly_") -> temp
es_yr %>% 
  dplyr::select(-monthly, -avg_monthly) %>%
  pivot_wider(names_from = year, 
              values_from = yr_total_samples, 
              names_prefix = "n_") %>% 
  full_join(temp) -> es_yr_reg

es_yr_reg %>% 
  rowwise() %>% 
  mutate(monthly_last5 = all(c_across(`monthly_2019`:`monthly_2023`)),
         monthly_gt5 = (monthly_last5 & 
                          any(c_across(`monthly_2014`:`monthly_2018`)))) %>%
  mutate(site_class = as.factor(case_when(monthly_gt5 ~ "Regular",
                                          monthly_last5 ~ "New",
                                          TRUE ~ "Sporadic"))) %>% 
  dplyr::select(site_id, site_class) -> es_yr_class

table(es_yr_class$site_class)
# Sporadic  Regular      New 
#      724       37        2 

# To do:
# Investigate actual frequency of testing versus monthly assumption

# Add class to main ES data
es <- es %>% full_join(es_yr_class)

# ---------------------------------------------------------------------------- #
# Save cleaned datasets

saveRDS(shape2, here(outdir,"shape2.rds"))
saveRDS(shape1, here(outdir,"shape1.rds"))
saveRDS(shape0, here(outdir,"shape0.rds"))

write_csv(data, here(outdir,"tsir_data.csv"))
write_csv(afp, here(outdir,"linelist_afp_clean.csv"))
write_csv(es, here(outdir,"linelist_es_clean.csv")) 

################################################################################
################################################################################