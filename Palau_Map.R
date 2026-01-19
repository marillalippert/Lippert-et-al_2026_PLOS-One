
library(tidyverse)
library(stars)
library(rnaturalearth)
library(ggspatial)
library(mapview)
library(cowplot)
library(ggnewscale)

in_situ_coords <- read.csv("Data/Nighttime_Full_InSitu.csv", stringsAsFactors = F)
in_situ_coords$Latitude <- as.numeric(in_situ_coords$Latitude)
in_situ_coords$Longitude <- as.numeric(in_situ_coords$Longitude)
sst_coords <- read.csv("Data/Regional_CRW_SST_Full.csv", stringsAsFactors = F)
colnames(sst_coords) <- c("long", "lat", "Date", "Temp", "Region")

# Convert colony coordinates to `sf` spatial points for plotting
colony_coords <- in_situ_coords |> filter(!is.na(Latitude)) |> 
  group_by(Region, Colony) |> 
  summarize(
    Latitude = unique(Latitude), 
    Longitude = unique(Longitude)) |> 
#  ) |> 
 # filter(!is.na(Latitude)) |> 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

# Calculate bounding box for SST coordinates, and use bounding
# boxes for each region to create `sf` spatial polygons for plotting
sst_coords_minmax <- sst_coords |> group_by(Region) |> 
  summarize(min_long = min(long) - 0.025,
            min_lat = min(lat) - 0.025,
            max_long = max(long) + 0.025, 
            max_lat = max(lat) + 0.025) |> 
  st_as_sf(coords = c("min_long", "min_lat"), remove = FALSE, crs = 4326) |> 
  rowwise() |> 
  mutate(
    geometry = 
      st_make_grid(setNames(unlist(c(min_long, min_lat, max_long, max_lat)), 
                            c("xmin","ymin","xmax","ymax")), 
                   n = 1, crs = 4326)
  )

# Read in bathymetry data, convert from elevation to depth by negating
bathy <- stars::read_stars("Data/Palau_bathymetry.tiff")
bathy$depth <- -bathy$Palau_bathymetry.tiff

# Create depth contour polygons from bathymetry data
# Set values above sea level to NA to use different color for land
depth_breaks <- c(0, 10, 50, 100, seq(2000, 8000, 1000))
depth_contours <- stars::st_contour(dplyr::select(bathy, depth), breaks = depth_breaks)
depth_contours$Min[is.infinite(depth_contours$Min)] <- NA

# Spatial buffer to make white background for colony coordinates
points_backdrop <- st_union(st_buffer(colony_coords, dist = 1.75e3))

# Make orthographic inset map
ortho_crs <- "+proj=ortho +lat_0=7.5 +lon_0=134.5"
world <- dplyr::select(ne_countries(scale = "large", returnclass = "sf"), geometry)

# Make border for globe
world_edge <- world |> st_bbox() |> st_make_grid(n = 100) |> st_as_sf() |> 
  st_cast("MULTILINESTRING") |> 
  st_cast("LINESTRING", do_split = TRUE) |> 
  st_transform(crs = ortho_crs) |> 
  mutate(npts = npts(x, by_feature = TRUE)) |> 
  filter(npts >= 4) |> 
  st_cast("POLYGON") |> st_union() |> st_convex_hull()

# Orthographic map of continents
world_ortho <- world |> 
  st_cast("MULTILINESTRING") |> 
  st_cast("LINESTRING", do_split = TRUE) |> 
  st_transform(crs = ortho_crs) |> 
  mutate(npts = npts(geometry, by_feature = TRUE)) |> 
  filter(npts >= 4) |> 
  st_cast("POLYGON")

# Point for Palau
palau_pt <- data.frame(x = 134.5, y = 7.5) |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  st_transform(ortho_crs)

# Inset map
world_inset <- ggplot() + 
  geom_sf(data = world_edge, fill = "white", color = "#012738", linewidth = 1) +
  geom_sf(data = world_ortho, fill = "#012738", color = "#012738") + 
  geom_sf(data = palau_pt, color = "red", size = 2.5) +
  theme_void()

palau_map <- ggplot() + 
  geom_sf(aes(fill = Min), data = depth_contours, linewidth = 0.2, color = NA) +
  geom_sf(aes(color = Min), data = depth_contours, linewidth = 0.2, fill = NA, show.legend = FALSE) +
  geom_sf(data = filter(depth_contours, is.na(Min)), color = colorspace::darken("#F5F5DC", 0.5), fill = NA) + 
  scale_fill_binned(breaks = depth_breaks, trans = "sqrt", na.value = "#F5F5DC",
                    type = "gradient", low = "#c2edff", high = "#012f42") + 
  scale_color_binned(breaks = depth_breaks, trans = "sqrt", 
                    type = "gradient", low = colorspace::darken("#c2edff", 0.2), 
                    high = colorspace::darken("#012f42", 0.3)) + 
  new_scale_color() +
  geom_sf(data = sst_coords_minmax, fill = NA, linewidth = 3, color = "white") + 
  geom_sf(data = points_backdrop, fill = "white", color = NA) + 
  geom_sf(aes(color = Region), data = sst_coords_minmax, fill = NA, linewidth = 1) + 
  geom_sf(aes(color = Region), data = colony_coords, size = 2) + 
  coord_sf(xlim = c(134.05, 134.95), ylim = c(6.9, 8.1)) + 
  annotation_scale(location = "br", text_col = "white", text_face = "bold") +
  annotation_north_arrow(location = "br", style = north_arrow_fancy_orienteering(text_col = "white"), 
                         height = unit(1, "cm"), pad_y = unit(0.75, "cm")) +
  scale_color_manual(values = colorspace::lighten(c("goldenrod3", "darkgreen", "magenta4"), 0.15)) +
  theme(
    axis.title = element_blank(), 
    legend.position.inside = c(0, 1), 
    legend.justification = c(0, 1), 
    legend.key = element_blank()
  ) + 
  labs(fill = "Depth (m)")

palau_map <- palau_map |> 
  ggdraw() + 
  draw_plot(world_inset, -0.015, 0.73, height = 0.25, width = 0.5)

ggsave("Figures/Main_Figs/Palau_map.tif", palau_map, height = 7, width = 7, units = "in", dpi = 300)

