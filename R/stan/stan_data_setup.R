pacman::p_load(tidyverse, sf, janitor, terra, exactextractr, mgcv, gstat, broom)
theme_set(theme_minimal())

# Inputs
indir <- "data/Pakistan"
proj_local <- "EPSG:32631"

# ES linelist
es <- readRDS(file.path(indir,"es.rds")) |> 
  mutate(nonpolio = is.na(final_class)|final_class %in% c("NPEV","VACCINE"))

# AFP linelist
afp <- readRDS(file.path(indir,"afp.rds"))

# Population raster
poprast  <- rast(file.path(indir,"proj_pop_rast.tif"))

# Administrative boundaries
shape2 <- readRDS(file.path(indir, "shape2.rds")) |> 
  st_transform(proj_local)

shape0 <- st_union(shape2)
