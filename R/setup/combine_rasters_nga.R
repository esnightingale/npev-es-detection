################################################################################
################################################################################
# Combine ESRI Sentinal-2 land cover tiles for Nigeria
# https://www.arcgis.com/home/item.html?id=cfcb7609de5f478eb7666240902d4d3d
################################################################################

library(terra)

# Nigeria country shape
shape0 <- readRDS("../polio-spec-models/inputs/Nigeria - manuscript/nga/shape2.rds") |> 
  st_transform(32631)

# Read individual tiles and create spatial raster collection
f <- paste0("data/",list.files("data", pattern = "20210101"))
rast_collect <- lapply(f, function(f) rast(f) |> project("EPSG:32631")) |> sprc()

# Define all wrt built area
built_fun <- function(x){
  x <- 
}

# Combine into one raster using mosaic
lc_nga <- mosaic(rast_collect, fun = Mode)

lc_nga_crop <- crop(lc_nga, shape0, snap="near")

plot(lc_nga_crop)

ggplot() +
  geom_raster(data = lc_nga_crop) + 
  geom_sf(data = shape0)

writeRaster(lc_nga, "data/esri_land_cover_nga.tif")

################################################################################
