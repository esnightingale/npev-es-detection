################################################################################
################################################################################
# Analysis data setup

# -> Save outputs for analysis
################################################################################

pacman::p_load(tidyverse, sf, janitor, units, here, lubridate, terra)

rm(list = ls(all = TRUE))  # clear all data

################################################################################

# Set parameters (from config)
source(here("R/config.R"))

filename_afp = here(dir,"linelist_afp_clean.csv")
filename_es = here(dir,"linelist_es_clean.csv")
filename_npafp = here(dir,"npafp_clean.csv")

# ----- #

# Read input data

## Surveillance data
### AFP case records
afp_ = read_csv(filename_afp) |> 
  filter(between(month, tp[1], tp[2]))

### Environmental surveillance sample records 
es_ = read_csv(filename_es) |> 
  filter(between(month, tp[1], tp[2]))

### NPAFP indicator for denominators
npafp_ind <- read_csv(filename_npafp) |> 
  filter(between(month, tp[1], tp[2]))

### Spatial data
shape2 <- readRDS(here(dir,"shape2.rds")) 
poprast_proj <- rast(here(dir,"WorldPop", "proj_pop_rast.tif"))
rain_proj <- rast(here(dir,"WorldPop", "proj_precip_rast.tif"))
built_proj <- rast(here(dir,"WorldPop", "proj_built_rast.tif"))

### Novel-T watersheds (run R/setup/setup_watersheds.Rmd first; for Nigeria create data/Nigeria/Novel-T/ and run with country="Nigeria")
watershed_path <- here(dir, "analysis/detection", "watershed.rds")
if (file.exists(watershed_path)) {
  catch_wsh <- readRDS(watershed_path)
} else {
  message("Watershed file not found at ", watershed_path, 
          ". Using 5km buffers only. Run R/setup/setup_watersheds.Rmd for watershed catchments.")
  # Empty sf without site_id to avoid duplicate column names in st_join
  catch_wsh <- st_sf(id = integer(0), population = numeric(0), geometry = st_sfc(crs = 4326))
}

################################################################################
# Set up data for analysis ------------------------------------------------

# AFP data ----------------------------------------------------------------

afp <- afp_ |> 
  # Flag 
  # - *any* EV detection (poliovirus, vaccine, or NPEV) 
  # - NPEV (with or without PV)
  # - PV only (PV or vaccine, without NPEV)
  # Distinguish polio vs non-polio for denominator/rate calculations
  mutate(month_of_year = month(month, label = TRUE, abbr = TRUE),
         ev = !is.na(virus_type_s),
         npev = grepl("NPEV", virus_type_s),
         pv_only = ev & !npev,
         npafp = !pv_only) 

# # Check AFP data has all districts/months:
# tmp <- expand.grid(
#     district = levels(npafp_ind$guid),
#     month_index = 1:12
#   )

afp_agg <- afp |>
  # Aggregate to district and month
  group_by(guid, month) |> 
  summarise(n_afp = n(),
            n_npafp = sum(npafp),
            n_ev = sum(ev),
            n_npev = sum(npev, na.rm = T),
            n_pv_only = sum(pv_only)) |> 
  ungroup() |> 
  # Join with U15 population denominator, expanding to all districts/months
  right_join(npafp_ind) |> 
  # Define rate/prevalence variables of interest and time indicators
  mutate(npafp_rate = n_npafp*1e5*12/denominator,
         ev_afp_p = n_ev/n_afp) |>
  rename(guid_old = guid) |> 
  # Add population density
  full_join(shape2 |> 
              st_drop_geometry() |> 
              select(adm1_name, guid, guid_old, guid_total_pop, guid_pop_dens)) 

# Now fill in zero-counts for district/months with no afp notifications, and
# add time variables
afp_agg <- mutate(afp_agg, across(c(n_afp:ev_afp_p), ~replace_na(., 0)),
                  t = as.numeric(as.factor(period)),
                  month_of_year = month(month, label = TRUE, abbr = TRUE)) 

summary(afp_agg)

# ES data -----------------------------------------------------------------

es = es_ %>% 
  dplyr::select(iso_3_code,admin_1,admin_2,guid,
                site_id, site_name, x, y, coord_imp,
                collection_date, travel_time, sample_condition,
                npev, virus_type_s, virus_cluster_s, emergence_group_s) |>
  # Define additional time variables, PV/EV positivity and site type
  mutate(site_name = toupper(site_name),
         month = floor_date(collection_date, "month"),
         month_of_year = month(collection_date, label = T),
         year = year(month),
         ev = !is.na(virus_type_s),
         npev = grepl("NPEV", virus_type_s),
         pv_only = ev & !npev,
         # Define factor for EV positivity
         ev_f = case_when(!ev ~ "EV-",
                            ev ~ "EV+"), 
         across(c(iso_3_code:site_id,sample_condition, ev_f), as.factor)) 

summary(es)

# ES sites 
sites <- es |> 
  mutate(type = case_when(grepl("TREATMENT",site_name) | grepl("PLANT", site_name) |
                             (grepl("PUMP", site_name) & grepl("STATION", site_name)) ~ "WWTP or Pumping station",
                          TRUE ~ "Other") |> 
      factor(levels = c("WWTP or Pumping station","Other"))) |> 
  group_by(guid, site_id, site_name, x, y, coord_imp, type) |> 
  summarise(first = min(collection_date),
            first_y = year(first),
            total = n(),
            p_pv_only = mean(pv_only),
            p_npev = mean(npev),
            p_ev = mean(ev)) |> 
  ungroup() |> 
  # Exclude samples from unmapped sites
  filter(!is.na(x)) |> 
  st_as_sf(coords = c("x","y"), crs = 4326, remove = F) 

tabyl(sites$type)
# sites$type   n   percent
# WWTP or Pumping station  80 0.1251956
# Other 559 0.8748044

# Define consistency class ------------------------------------------------

# Consider:
# + Ad-hoc sites at least 6m old but only sampled for a short period
# + Established regularly sampled sites, at least 1 sample/month for the time period
# + Newer regularly sampled sites, at least 1 sample per month since 2023

# First plot out all sites in the dataset
ggplot(es, aes(collection_date, site_id, col = as.numeric(site_id))) +
  geom_point() + 
  theme(axis.text.y = element_blank()) + 
  scale_colour_viridis_c() + 
  guides(col = "none") +
  labs(x = "Sample collection date", y = "Site")

per <- range(es$month)
st <- per[1]
tdy <- ceiling_date(per[2],"month")

# Summarise by month of collection
ids <- unique(es$site_id)
months <- seq(st, tdy, by='months')
months <- data.frame(site_id = rep(ids, each = length(months)),
                     month = rep(months, times = length(ids)))

# By month of collection
es |> 
  group_by(site_id, month) %>% 
  count() |> 
  ungroup() %>% 
  right_join(months) %>% 
  mutate(year = year(month),
         n = replace_na(n, 0)) |> 
  right_join(select(sites, site_id, first, total)) |> 
  group_by(site_id) |> 
  mutate(
    # Flag as "active" after first collection
    active = month >= floor_date(first,"month"),
    # Average sampling rate since first collection (whole 4y period)
    avg_rate_active = mean(n[active])) -> es_mth

# By year of collection
es_mth %>%   
  group_by(site_id, first, total, year, avg_rate_active) %>% 
  summarise(
    # No. active months
    active_mths = sum(active),
    # Total samples 
    yr_total = sum(n),
    # No. months with any sample
    sampled_mths = sum(n>0),
    # No. months with any sample since first collection
    sampled_mths_active = sum(n[active]>0),
    # Average sampling rate since first collection
    avg_rate_active = mean(n[active])) |> 
  ungroup() |> 
  mutate(
    # Indicator - sample for every month during this period?
    all_mths = sampled_mths == 12,
    # Indicator - sample for every month since first collection?
    all_mths_active = sampled_mths_active == active_mths,
    # Indicator - sample for at least half of months during this period?
    semi_reg = sampled_mths >= 6,
    # Indicator - sample for at least half of observed months since first collection?
    semi_reg_active = sampled_mths_active >= floor(active_mths/2)) -> es_yr

es_yr |> 
  dplyr::select(site_id:year, all_mths) |> 
  pivot_wider(names_from = year, 
              values_from = all_mths, 
              names_prefix = "all_mths") -> temp
es_yr |> 
  dplyr::select(site_id:year,all_mths_active) |> 
  pivot_wider(names_from = year, 
              values_from = all_mths_active, 
              names_prefix = "all_mths_active") -> temp2
es_yr |> 
  dplyr::select(site_id:year,semi_reg) |> 
  pivot_wider(names_from = year, 
              values_from = semi_reg, 
              names_prefix = "semi_reg") -> temp3
es_yr |> 
  dplyr::select(site_id:year,semi_reg_active) |> 
  pivot_wider(names_from = year, 
              values_from = semi_reg_active, 
              names_prefix = "semi_reg_active") -> temp4

es_yr |> 
  dplyr::select(site_id:year,yr_total) |> 
  pivot_wider(names_from = year, 
              values_from = yr_total, 
              names_prefix = "n") |> 
  full_join(temp) |> full_join(temp2) |> full_join(temp3) |> full_join(temp4) -> es_yr_reg

# Define regularity of sampling
es_yr_reg %>% 
  rowwise() %>% 
  mutate(
    # All years at least monthly sampling
    monthly_all = all(`all_mths2021`:`all_mths2024`),
    # At least monthly sampling since first collection
    monthly_active = all(`all_mths_active2021`:`all_mths_active2024`),
    # All years semi-regular
    semireg = all(`semi_reg2021`:`semi_reg2024`),
    # Semi-regular since first collection
    semireg_active = all(`semi_reg_active2021`:`semi_reg_active2024`),
    # Monthly for any year
    # monthly_st = any(`all_mths2021`:`all_mths2024`),
    # Sporadic/short-term sampling - at least 12m active but less than 10 samples in total, or semi-regular for max two years
    sporadic = (between(total, 2, 12) & first < as.POSIXct("2024-01-01")),
    short_term = sum(`semi_reg_active2021`:`semi_reg_active2024`) == 1,
    singleton = total == 1) |> 
  mutate(site_class = factor(case_when(
    singleton ~ "Singleton",
    monthly_all ~ "Monthly (whole period)",
    monthly_active ~ "Monthly (since established)",
    # monthly_st ~ "Monthly for limited time",
    semireg ~ "Semi-regular (whole period)",
    semireg_active ~ "Semi-regular (since established)",
    sporadic | short_term ~ "Sporadic or short-term",
    TRUE ~ "Other"),
    levels = c("Monthly (whole period)",
               "Semi-regular (whole period)",
               "Monthly (since established)",
               "Semi-regular (since established)",
               "Sporadic or short-term",
               "Singleton",
               "Other"))) |> 
  select(starts_with("site")) -> es_yr_class

table(es_yr_class$site_class)
# Monthly (whole period)      Semi-regular (whole period)      Monthly (since established) Semi-regular (since established)           Sporadic or short-term                        Singleton                            Other 
# 54                                9                               59                                3                              180                              334                                0 

# Add class to main ES data and site data
es <- full_join(es, es_yr_class)
sites <- left_join(sites, es_yr_class)

# Catchments: Radial (5km) ------------------------------------------------

# Catchment size
catch_5k <- sites %>% 
  select(site_id, x, y, geometry) |> 
  st_transform(proj_local) %>% 
  st_buffer(dist = units::set_units(5, "km"))

catch_5k$catchment_pop <- exactextractr::exact_extract(poprast_proj,
                                                        catch_5k, 
                                                        fun = "sum")  

# Also categorise by radial catchment size
catch_5k$catchment_cat <- factor(catch_5k$catchment_pop > 1e5, 
                                 levels = c(FALSE,TRUE),
                                    labels = c("<100,000", ">100,000"))

tabyl(catch_5k$catchment_cat)
# catch_5k$catchment_cat   n   percent
# <100,000 525 0.8215962
# >100,000 114 0.1784038

# Surrounding geography
## Precipitation
catch_5k$precip21 <- exactextractr::exact_extract(rain_proj, 
                                                  catch_5k, 
                                                  fun = "mean")  
## Urbanicity - % built surface area
tmp <- terra::extract(built_proj, catch_5k, fun = sum, na.rm = T)[,-1]
names(tmp) <- paste0("built",years[1]:years[2])

area_5k <- st_area(catch_5k)[1]
catch_5k <- cbind(catch_5k, tmp) |> 
  mutate(across(starts_with("built"), \(x) x/area_5k))

# Return buffers to unprojected for consistency
catch_5k <- st_transform(catch_5k, 4326)

# Catchments: Watershed ---------------------------------------------------
# Currently bugged for Pakistan watersheds - leaving for now as not currently using

if (nrow(catch_wsh) > 0) {
  # Rename novel-t population estimate for consistency, then also extract estimate from worldpop raster
  catch_wsh <- rename(catch_wsh, catchment_pop_nt = population) |> 
    st_transform(proj_local)
  catch_wsh$catchment_pop_wpop <- exactextractr::exact_extract(poprast_proj,
                                                            catch_wsh, 
                                                            fun = "sum")  
  catch_wsh <- st_transform(catch_wsh, 4326)
  ggplot(catch_wsh |> pivot_longer(starts_with("catchment_pop_")), aes(value, fill = name)) + 
    geom_density(alpha = 0.5) + 
    scale_x_continuous(trans = "log")
  plot(log(catch_wsh$catchment_pop_nt), log(catch_wsh$catchment_pop_wpop))
  abline(0,1)
}

# Connect sites/5k buffer to associated watersheds (where available)
sites_wsh <- st_join(sites, catch_wsh, 
                       st_is_within_distance, dist = set_units(100,"m"))
if (nrow(catch_wsh) > 0) {
  sites_wsh <- sites_wsh |> group_by(site_id) |> slice_head(n = 1) |> ungroup()
}
sites_wsh <- sites_wsh |> 
  left_join(select(st_drop_geometry(catch_5k), site_id, catchment_pop))

if (nrow(sites_wsh) > 0 && "catchment_pop_nt" %in% names(sites_wsh)) {
  plot(log(sites_wsh$catchment_pop), log(sites_wsh$catchment_pop_nt))
  abline(0,1)
}

wsh_sites <- if (nrow(catch_wsh) > 0 && "site" %in% names(catch_wsh)) {
  catch_wsh |>
    filter(!is.na(site)) |> 
    right_join(st_drop_geometry(sites_wsh |> select(site_id, id)))
} else {
  st_sf(geometry = st_sfc(crs = 4326))
}

adm1_example <- unique(shape2$adm1_name)[1]
site_pattern <- paste0(iso_code, "/")
if (nrow(wsh_sites) > 0 && nrow(sites_wsh) > 0) {
  ggplot() +
    geom_sf(data = shape2 |> filter(adm1_name == adm1_example)) +
    geom_sf(data = wsh_sites |> filter(grepl(site_pattern, site_id)),
            aes(geometry = geometry), fill = "blue") +
    geom_sf(data = sites_wsh |> filter(grepl(site_pattern, site_id)),
            aes(geometry = geometry), col = "red")
}
# 
# catch_wsh_dist <- st_intersection(catch_wsh, select(shape2, guid))
# 
# catch_dist_max <- catch_wsh_dist |> 
#   mutate(area = st_area(catch_dist)) |> 
#   st_drop_geometry() |> 
#   group_by(site_id) |> 
#   mutate(area_p = as.numeric(area/sum(area))) |>
#   slice_max(order_by = area_p) 
# 
# catch_dist_max |> 
#   mutate(area_p = cut(area_p, c(0,.5,.6,.7,.8,.9,1))) |> 
#   janitor::tabyl(area_p)

# Save analysis datasets --------------------------------------------------

outdir <- paste0(dir,"/analysis")
dir.create(here(outdir, "prevalence"), recursive = TRUE, showWarnings = FALSE)
dir.create(here(outdir, "detection"), recursive = TRUE, showWarnings = FALSE)

# AFP cases
saveRDS(afp, here(outdir, "prevalence/afp_linelist.rds"))
saveRDS(afp_agg, here(outdir, "prevalence/afp_analysis.rds"))

# ES 
saveRDS(es, here(outdir, "detection/es_analysis.rds"))
saveRDS(sites, here(outdir, "detection/sites_analysis.rds"))

# Catchment shapefiles and population sizes
saveRDS(catch_5k, here(outdir,"detection/catch_5k.rds"))
saveRDS(catch_wsh, here(outdir,"detection/catch_wsh.rds"))

################################################################################
################################################################################