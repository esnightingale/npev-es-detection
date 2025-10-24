################################################################################
# Aside: Divide into grid for Novel-T extraction
################################################################################

pacman::p_load(tidyverse, sf, janitor, terra, exactextractr, mgcv, gstat, broom)
theme_set(theme_minimal())

# Inputs
indir <- "data/Pakistan"
proj_local <- "EPSG:32631"

# ES linelist
es <- readRDS(file.path(indir,"es.rds"))

# Administrative boundaries
shape2 <- readRDS(file.path(indir, "shape2.rds")) |> 
  st_transform(proj_local)

shape0 <- st_union(shape2)

# Define sites ------------------------------------------------------------

es |> 
  select(guid, admin_0:admin_2, site_id, site_name, x, y) |> 
  unique() -> es_sites

glimpse(es_sites)

# Save unique sites for extraction of novel-t catchments
es_sf <- st_as_sf(es_sites, coords = c("x","y"), crs = 4326) |> mutate(id = row_number())
st_write(es_sf |> select(-site_id, -site_name), "data/Pakistan/es_sites_pak.shp", append = F)

# Option 1: Make grid ------------------------------------------------------

# Make grid from the extents
pk_grid <- st_make_grid(shape0 |> st_transform(4326), cellsize = c(0.9,1))

# Define a 10km buffer around sites (assuming catchment will extend no more than 10km from site location)
es_union_buff <- st_buffer(es_sf, units::set_units(10, "km")) |> 
  st_union() |> 
  st_make_valid()

ggplot() +
  geom_sf(data = shape0|> st_transform(4326), fill = "grey", color = NA) +
  geom_sf(data = pk_grid, fill = NA) +
  geom_sf(data = es_union_buff, fill = "blue", alpha = 0.2) +
  geom_sf(data = es_sf, col = "blue", cex = 0.5) +
  theme_void()

# Filter to only grid cells that intersect 50km buffer
pk_grid_filt <- pk_grid[lengths(st_intersects(pk_grid, es_union_buff)) > 0, ]

ggplot() +
  geom_sf(data = shape0|> st_transform(4326), fill = "grey", color = NA) +
  geom_sf(data = pk_grid_filt, fill = NA) +
  geom_sf(data = es_union_buff, fill = "blue", alpha = 0.2) +
  geom_sf(data = es_sf, col = "blue") +
  theme_void()

st_write(pk_grid_filt, "data/Pakistan/grid_partition.shp", append = F)

# Option 2: K-means clustering --------------------------------------------

points <- st_coordinates(es_sf)

# Apply K-Means clustering 
set.seed(123)
kmeans_result <- kmeans(points, centers = 90)

print(kmeans_result$cluster)

es_sf$cluster_km <- as.factor(kmeans_result$cluster)

ggplot() +
  geom_sf(data = shape0|> st_transform(4326), fill = "grey", color = NA) +
  geom_sf(data = es_sf, aes(col = cluster_km)) + 
  scale_colour_viridis_d(option= "turbo") + 
  guides(col = "none") + 
  theme_void()

tabyl(kmeans_result$cluster) |> arrange(desc(n))

# Define buffered bbox ----------------------------------------------------

es_buff <- st_buffer(es_sf, units::set_units(10, "km"))

# Compute bounding boxes for each cluster
# As sf object
bbox_sf <- do.call(rbind, lapply(split(es_buff, es_buff$cluster_km), function(cluster_data) {
  bbox <- st_bbox(cluster_data) |> st_as_sfc() |> st_as_sf()
  return(bbox)
}))

ggplot() +
  geom_sf(data = shape0|> st_transform(4326), fill = "grey", color = NA) +
  geom_sf(data = bbox_sf, fill = NA) +
  geom_sf(data = es_sf, aes(col = cluster_km)) + 
  scale_colour_viridis_d(option= "turbo") + 
  guides(col = "none") +
  theme_void()

clust_perim <- st_perimeter(bbox_sf)
summary(clust_perim/4)

# Write shapefile of cluster extents
st_write(bbox_sf, "data/Pakistan/km_cluster_bbox.shp", append = F)


bounding_boxes <- do.call(rbind, lapply(split(es_buff, es_buff$cluster_km), function(cluster_data) {
  print(cluster_data)
  bbox <- st_bbox(cluster_data)
  return(data.frame(
    cluster = unique(cluster_data$cluster),
    min_lon = bbox["xmin"], max_lon = bbox["xmax"],
    min_lat = bbox["ymin"], max_lat = bbox["ymax"]
  ))
}))

# Print bounding boxes
print(bounding_boxes)
