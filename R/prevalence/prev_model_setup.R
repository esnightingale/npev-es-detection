################################################################################
################################################################################
# Set up data for EV prevalence model
################################################################################

source(here::here("R/config.R"))
dir <- file.path(dir, "analysis/prevalence")

# AFP analysis data
data <- readRDS(here(dir, "afp_analysis.rds"))

# Shapefiles 
shape2 <- readRDS(here(dir,"../shape2.rds")) |> st_transform(4326)

# Set up data
preddata <- mutate(data, 
                   # Scale time variables
                   time = t,
                   t = (time - min(time))/diff(range(time)),
                   m = as.numeric(month_of_year),
                   month_of_year = (m - min(m))/diff(range(m)),
                   # standardise pop density for stability
                   log_pop_dens = log(guid_pop_dens) |> as.numeric(),
                   log_pop_dens_s = (log_pop_dens - mean(log_pop_dens))/sd(log_pop_dens)) |> 
  rename(guid_old = guid) |> 
  left_join(shape2 |> st_drop_geometry() |> select(guid_old, guid))

# Exclude zero trial rows for fitting (will still predict for these later) 
fitdata <- filter(preddata, n_afp > 0)
write_rds(preddata, here(dir, "preddata.rds"))
write_rds(fitdata, here(dir, "fitdata.rds"))

# Define adm1-guid key for later matching:
prov_guid <- preddata |> 
  select(adm1_name, guid) |> 
  distinct()

################################################################################
################################################################################