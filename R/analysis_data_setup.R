################################################################################
################################################################################
# Analysis data setup

# -> Save outputs for analysis
################################################################################

pacman::p_load(tidyverse, sf, janitor, units, here, lubridate, terra)

rm(list = ls(all = TRUE))  # clear all data

################################################################################

# Set parameters
## Input file paths
dir <- "data/Pakistan"
filename_afp = here(dir,"linelist_afp_clean.csv")
filename_es = here(dir,"linelist_es_clean.csv")
filename_npafp = here(dir,"npafp_clean.csv")

## Time period of interest
tp <- ymd(c("2021-01-01","2024-12-01"))

## Local projection
proj_local <- "EPSG:32642"

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

## Spatial data
### District shapefiles
shape2 <- readRDS(here(dir,"shape2.rds")) 

### Population raster
poprast <- rast(here(dir,"WorldPop","pak_ppp_2020_UNadj_constrained.tif")) 

### Novel-T watersheds
catch_wsh <- readRDS(here(dir,"analysis", "watershed.rds"))

### Annual mean precipitation https://hub.worldpop.org/doi/10.5258/SOTON/WP00772 
rain <- rast(here(dir,"WorldPop","pak_ppt_2021_yravg_tc_100m_v1.tif"))

### Built up surface https://hub.worldpop.org/doi/10.5258/SOTON/WP00772 
f <- list.files(here(dir, "WorldPop"), 
                full.names = T,
                pattern = "pak_built_S_GHS_U_wFGW_100m_v1_")
built <- rast(f)

################################################################################
# Set up data for analysis ------------------------------------------------

# Shapefiles --------------------------------------------------------------

## Projected data
shape2_proj <- st_transform(shape2, proj_local) 
poprast_proj <- project(poprast, proj_local)
rain_proj <- project(rain, proj_local)
built_proj <- project(built, proj_local) 

# Define district area in km2
shape2$shape_area_km2  <- st_area(shape2_proj) |> units::set_units("km^2")
hist(shape2$shape_area_km2)
summary(shape2$shape_area_km2)

## Extract total district populations from raster
shape2$guid_total_pop <- exactextractr::exact_extract(poprast_proj,
                                                 shape2_proj, 
                                                 fun = "sum") 

# Define district population density
shape2$guid_pop_dens <- shape2$guid_total_pop/shape2$shape_area_km2
hist(shape2$guid_pop_dens)  
summary(shape2$guid_pop_dens)

# ggplot(shape2) +
#   geom_sf(aes(fill = as.numeric(guid_pop_dens))) +
#   scale_fill_viridis_c(name = paste0("Pop Density (", 
#                                      deparse_unit(shape2$guid_pop_dens), ")"),
#                        trans = "log10",
#                        labels = scales::comma) +
#   theme_void()

# AFP data ----------------------------------------------------------------

afp <- afp_ |> 
  # Distinguish polio and non-polio cases
  mutate(npafp = !(grepl("WILD", virus_type_s) | grepl("VDPV",virus_type_s))) |> 
  # exclude AFP among over-15s 
  filter(age_m < 15*12) 

# # Check NPAFP data has all districts/months:
# tmp <- expand.grid(
#     district = levels(npafp_ind$guid),
#     month_index = 1:12
#   )

afp_agg <- afp |>
  # Aggregate to district and month
  group_by(guid, month) |> 
  summarise(n_afp = n(),
            n_npafp = sum(npafp, na.rm = T),
            n_npev = sum(npev == "Yes", na.rm = T)) |> 
  ungroup() |> 
  # Join with U15 population denominator, expanding to all districts/months
  right_join(npafp_ind) |>
  # Define rate/prevalence variables of interest and time indicators
  mutate(npafp_r = n_npafp*1e5/denominator,
         npev_npafp_p = n_npev/n_npafp,
         npev_afp_p = n_npev/n_afp,
         t = as.numeric(as.factor(period)),
         month_of_year = month(month, label = TRUE, abbr = TRUE)) |> 
  # Add population density
  full_join(shape2 |> 
              st_drop_geometry() |> 
              select(adm1_name, guid, guid_total_pop, guid_pop_dens), by = "guid") 

# Note: Some denominators are imputed (all for 2021)

# Check agreement of npafp numerator and tallied npafp case records:
afp |> 
  mutate(check = n_npafp == numerator,
         diff = n_npafp - numerator) -> tmp

tabyl(tmp$check)
hist(tmp$diff)

# => Most differences are small, generally numerator > tallied records 
#    - Could be slight differences in the way we've defined polio cases?
#    - Only excluded 6 due to missing guid so that's not it
  
# Now fill in zero-counts for district/months with no afp notifications
afp <- mutate(afp, across(c(n_afp:npev_afp_p), ~replace_na(., 0))) 

summary(afp)

# ES data -----------------------------------------------------------------

es = es_ %>% 
  dplyr::select(iso_3_code,admin_1,admin_2,guid,
                site_id, site_name, x, y, coord_imp,
                collection_date, travel_time, sample_condition, final_class,
                npev, virus_type_s, virus_cluster_s, emergence_group_s) |>
  # Define additional time variables, PV/NPEV positivity and site type
  mutate(month = floor_date(collection_date, "month"),
         month_of_year = month(collection_date, label = T),
         year = year(month),
         # Flag any samples in which poliovirus was detected to later remove
         pv = grepl("WILD",virus_type_s) | grepl("PV", virus_type_s),
         # Define factor for NPEV positivity
         npev_f = case_when(is.na(npev) ~ "NPEV-",
                            !is.na(npev) ~ "NPEV+"), 
         site_type = case_when((grepl("PUMPING",site_name)|
                                  (grepl("PUMP", site_name) & grepl("STATION", site_name))|
                                  grepl("PLANT", site_name)) ~ "WWTP or Pumping station",
                               TRUE ~ "Other") |> 
           factor(levels = c("WWTP or Pumping station","Other")),
         across(c(iso_3_code:site_id,sample_condition, npev, final_class), as.factor)) |> 
  # Filter to sites with at least min_ss samples
  group_by(site_id) |> 
  mutate(site_total = n()) |> 
  ungroup() 

summary(es)
tabyl(es$site_type)

# ES sites 

sites <- es |> 
  group_by(guid, site_id, x, y) |> 
  summarise(first_collection = min(collection_date),
            first_collection_y = year(first_collection),
            total_collections = n(),
            p_wpv = mean(final_class == "WPV1"),
            p_npev = mean(!is.na(npev))) |> 
  ungroup() |> 
  st_as_sf(coords = c("x","y"), crs = 4326, remove = F) 

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
es %>% 
  group_by(site_id) |> 
  summarise(site_first = min(collection_date),
            total = n()) |> 
  ungroup() -> site_first

# By month of collection
es |> 
  group_by(site_id, month) %>% 
  count() |> 
  ungroup() %>% 
  right_join(months) %>%
  mutate(year = year(month),
         n = replace_na(n, 0)) |> 
  right_join(site_first) |> 
  group_by(site_id) |> 
  mutate(
    # Flag as "active" after first collection
    active = month >= floor_date(site_first,"month"),
    # Average sampling rate since first collection (whole 4y period)
    avg_rate_active = mean(n[active])) -> es_mth

# By year of collection
es_mth %>%   
  group_by(site_id, site_first, total, year, avg_rate_active) %>% 
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
    sporadic = (between(total, 2, 12) & site_first < as.POSIXct("2024-01-01")),
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
sites <- full_join(sites, select(es_yr_class, site_id, site_class))


# Catchments: Radial (5km) ------------------------------------------------

# Catchment size
catch_5k <- sites %>% 
  select(site_id, x, y, geometry) |> 
  st_transform(proj_local) %>% 
  st_buffer(dist = units::set_units(5, "km"))

sites$catchment_pop_5k <- catch_5k$catchment_pop_5k <- exactextractr::exact_extract(poprast_proj,
                                                        catch_5k, 
                                                        fun = "sum")  

# Also add to ES dataset
es <- left_join(es, select(sites, site_id, catchment_pop_5k)) |> 
  # Categorise by radial catchment size
  mutate(catchment_cat_5k = factor(catchment_pop_5k > 1e5, 
                            labels = c("<100,000", ">100,000")))

tabyl(es$catchment_cat_5k)

# Surrounding geography
## Precipitation
sites$precip21 <- catch_5k$precip21 <- exactextractr::exact_extract(rain_proj, 
                                                                    catch_5k, 
                                                                    fun = "mean")  
## Urbanicity - % built surface area
tmp <- terra::extract(built_proj, catch_5k, fun = sum)[,-1]
names(tmp) <- paste0("built",21:24)

area_5k <- st_area(catch_5k)[1]
sites <- cbind(sites, tmp) |> 
  mutate(across(starts_with("built"), \(x) x/area_5k))
catch_5k <- cbind(catch_5k, tmp) |> 
  mutate(across(starts_with("built"), \(x) x/area_5k))

# Return buffers to unprojected for consistency
catch_5k <- st_transform(catch_5k, 4326)

# Catchments: Watershed ---------------------------------------------------

# Rename novel-t population estimate for consistency, then also extract estimate from worldpop raster
catch_wsh <- rename(catch_wsh, catchment_pop_nt = population) |> 
  st_transform(proj_local)

catch_wsh$catchment_pop_wpop <- exactextractr::exact_extract(poprast_proj,
                                                          catch_wsh, 
                                                          fun = "sum")  
# Consists of vary many watersheds (partitioning the whole area in extracted tiles, not just attributed to sites?)
# Can't directly identify which link to sites
# Filter these before extracting other spatial covariates

# # Surrounding geography
# ## Precipitation
# catch_wsh$precip21 <- exactextractr::exact_extract(rain_proj, 
#                                                    catch_wsh, 
#                                                    fun = "mean")  
# ## Urbanicity - % built surface area
# tmp <- extract(built_proj, catch_wsh, fun = mean)[,-1]
# names(tmp) <- paste0("built",21:24)
# 
# catch_wsh <- cbind(catch_wsh, tmp) |> 
#   mutate(across(starts_with("built"), \(x) x/(1000^2*area)))

# Return to unprojected for consistency
catch_wsh <- st_transform(catch_wsh, 4326)

ggplot(catch_wsh |> pivot_longer(starts_with("catchment_pop_")), aes(value, fill = name)) + 
  geom_density(alpha = 0.5) + 
  scale_x_continuous(trans = "log")

plot(log(catch_wsh$catchment_pop_nt), log(catch_wsh$catchment_pop_wpop))
abline(0,1)
# Novel-T populations somewhat smaller than populations extracted from wpop raster to the same polygons, and less varied.

# Connect sites to associated watersheds (where available)
sites_wsh <- st_join(sites, catch_wsh, 
                       st_is_within_distance, dist = set_units(100,"m")) |> 
  group_by(site_id) |> 
  slice_head(n = 1)

plot(log(sites_wsh$catchment_pop_5k), log(sites_wsh$catchment_pop_nt))
abline(0,1)
# Novel-T catchments overall smaller in pop size than 5km radius

# tabyl(sites_wsh, site_id) |> View()

wsh_sites <- catch_wsh |>
  filter(!is.na(site)) |> 
  right_join(st_drop_geometry(sites_wsh |> select(site_id, id)))

ggplot() +
  geom_sf(data = shape2 |> filter(adm1_name == "AJK")) +
  geom_sf(data = wsh_sites |> filter(grepl("PAK/AJ",site_id)),
          aes(geometry = geometry), fill = "blue") +
  geom_sf(data = sites_wsh |> filter(grepl("PAK/AJ",site_id)),
          aes(geometry = geometry), col = "red")
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

# AFP cases
saveRDS(afp, here(outdir, "afp_analysis.rds"))

# ES 
saveRDS(es, here(outdir, "es_analysis.rds"))
saveRDS(sites, here(outdir, "sites_analysis.rds"))

# Catchment shapefiles and population sizes
saveRDS(catch_5k, here(outdir,"catch_5k.rds"))
saveRDS(catch_wsh, here(outdir,"catch_wsh.rds"))

# Shapefiles
saveRDS(shape2, here(outdir, "shape2.rds"))

# Projected rasters
writeRaster(poprast_proj, 
            here(dir,"WorldPop", "proj_pop_rast.tif"), 
                 overwrite=T)
writeRaster(rain_proj, 
            here(dir,"WorldPop", "proj_precip_rast.tif"), 
            overwrite=T)
writeRaster(built_proj, 
            here(dir,"WorldPop", "proj_built_rast.tif"), 
            overwrite=T)

################################################################################
################################################################################