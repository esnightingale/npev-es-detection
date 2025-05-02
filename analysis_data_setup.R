################################################################################
################################################################################
# Analysis data setup

# -> Save outputs for analysis
################################################################################

library(tidyverse)
library(sf)
library(janitor)
library(units)
library(here)
library(lubridate)
library(terra)

rm(list = ls(all = TRUE))  # clear all data

dir <- "data/Pakistan"

# Local projection
proj_local <- "EPSG:32642"

# Inputs
filename_afp = here(dir,"linelist_afp_clean.csv")
filename_es = here(dir,"linelist_es_clean.csv")
filename_npafp = here(dir,"npafp_clean.csv")

################################################################################

# Time period of interest
tp <- ymd(c("2021-01-01","2024-12-01"))

# AFP case records
afp_ = read_csv(filename_afp) |> 
  filter(between(month, tp[1], tp[2]))

# Environmental surveillance sample records 
es_ = read_csv(filename_es) |> 
  filter(between(month, tp[1], tp[2]))

# NPAFP indicator for denominators
npafp_ind <- read_csv(filename_npafp) |> 
  filter(between(month, tp[1], tp[2]))

# District shapefiles
shape2 <- readRDS(here(dir,"shape2.rds"))

# Population raster
poprast <- rast(here(dir,"pak_ppp_2020_UNadj_constrained.tif"))

# Set up data for analysis ------------------------------------------------

# Shapefiles --------------------------------------------------------------

shape2_proj <- st_transform(shape2, proj_local) 

shape2$shape_area_km2  <- st_area(shape2_proj) |> units::set_units("km^2")
hist(shape2$shape_area_km2)
summary(shape2$shape_area_km2)

# Extract total district populations from raster
shape2$guid_total_pop <- exactextractr::exact_extract(poprast,
                                                 shape2_proj, 
                                                 fun = "sum") 

shape2$guid_pop_dens <- shape2$guid_total_pop/shape2$shape_area_km2
hist(shape2$guid_pop_dens)  
summary(shape2$guid_pop_dens)

ggplot(shape2) +
  geom_sf(aes(fill = as.numeric(guid_pop_dens))) +
  scale_fill_viridis_c(name = paste0("Pop Density (", 
                                     deparse_unit(shape2$guid_pop_dens), ")"),
                       trans = "log10",
                       labels = scales::comma) +
  theme_void()

# Define spatial adjacency matrix
nb <- spdep::poly2nb(shape2)  # Create neighbors list 
W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE) |> as.matrix() # Convert to adjacency matrix
rownames(W) <- shape2$guid # Set rownames to district names

# AFP data ----------------------------------------------------------------

afp <- afp_ |> 
  # Distinguish polio and non-polio cases
  mutate(npafp = !(grepl("WPV", virus_type_s) | grepl("VDPV",virus_type_s))) |> 
  # exclude AFP among over-15s 
  filter(age_m < 15*12) |>
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
afp <- mutate(afp,
              across(c(n_afp:npev_afp_p), ~replace_na(., 0))) 

summary(afp)

# ES data -----------------------------------------------------------------

es = es_ %>% 
  dplyr::select(iso_3_code,admin_1,admin_2,guid,
                site_id, site_name, x, y, coord_imp,
                collection_date, travel_time, sample_condition, final_class,
                npev, virus_type_s, virus_cluster_s, emergence_group_s) |>
  # Define site type
  mutate(site_type = case_when((grepl("PUMPING",site_name)|
                                  (grepl("PUMP", site_name) & grepl("STATION", site_name))|
                                  grepl("PLANT", site_name)) ~ "WWTP or Pumping station",
                               TRUE ~ "Other") |> 
           factor(levels = c("WWTP or Pumping station","Other")),
         across(c(iso_3_code:site_id,sample_condition, npev, final_class), as.factor))

summary(es)
tabyl(es$site_type)

# ES catchments -----------------------------------------------------------

poprast <- project(poprast, proj_local)

# Unique sites
es %>% 
  select(iso_3_code:site_name, x, y, coord_imp, guid) %>% 
  distinct() -> es_sites

# Estimate 5km catchment population
es_catch5 <- es_sites %>% 
  st_as_sf(coords = c("x","y"), crs = 4326) %>% 
  st_transform(proj_local) %>% 
  st_buffer(dist = units::set_units(5, "km"))

es_sites$catchment_pop_5k <- exactextractr::exact_extract(poprast,
                                                        es_catch5, 
                                                        fun = "sum")  

es <- left_join(es, select(es_sites, site_id, guid, catchment_pop_5k))

# Save analysis datasets --------------------------------------------------

# AFP cases
saveRDS(afp, here(dir, "afp_analysis.rds"))

# ES 
saveRDS(es, here(dir, "es_analysis.rds"))
saveRDS(es_sites, here(dir, "es_sites.rds"))

# Shapefiles
saveRDS(shape2, here(dir,"shape2.rds"))
saveRDS(W, here(dir,"W.rds"))

# Projected population raster
writeRaster(poprast, 
            here(dir,"proj_pop_rast.tif"), 
                 overwrite=T)

################################################################################
################################################################################