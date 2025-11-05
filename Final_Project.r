library(sf)
library(spatstat)
library(tidyverse)

# Load block group and MSA files

blckgrp2010 <- st_read("nhgis0011_shape/nhgis0011_shapefile_tl2010_us_blck_grp_2010/US_blck_grp_2010.shp")
blckgrp2020 <- st_read("nhgis0011_shape/nhgis0011_shapefile_tl2020_us_blck_grp_2020/US_blck_grp_2020.shp")
msa2010 <- st_read("nhgis0011_shape/nhgis0011_shapefile_tl2010_us_cbsa_2010/US_cbsa_2010.shp")
msa2020 <- st_read("nhgis0011_shape/nhgis0011_shapefile_tl2020_us_cbsa_2020/US_cbsa_2020.shp")

# Filter by geos

denver2010 <- msa2010 %>%
  filter(GEOID10 == 19740)
denver2020 <- msa2020 %>%
  filter(GEOID == 19740)

denver2010 <- st_transform(denver2010, st_crs(blckgrp2010)) %>%
  st_union() %>%
  st_make_valid()

blckgrp2010.denver <- st_intersection(blckgrp2010, denver2010)

denver2020 <- st_transform(denver2020, st_crs(blckgrp2020)) %>%
  st_union() %>%
  st_make_valid()

blckgrp2020.denver <- st_intersection(blckgrp2020, denver2020)

plot(st_geometry(blckgrp2010.denver))
plot(st_geometry(blckgrp2020.denver))

# Load RTD files and filter only the commuter rail lines

RTD <- st_read("Local_transit/Local_transit.shp")

RTD.commuter <- RTD %>%
  filter(agency_id == "RTD" & route_shor %in% c("A","B","G","N"))

RTD.stations <- st_read("LightrailStations/LightrailStations.shp")

RTD.commuter.stations <- RTD.stations %>%
  # Keep only commuter lines
  filter(str_detect(RAIL_LINE, "A|B|G|N")) %>%
  
  # Separate multi-line entries into multiple rows
  separate_rows(RAIL_LINE, sep = "-") %>%
  mutate(RAIL_LINE = str_trim(RAIL_LINE)) %>%
  filter(RAIL_LINE %in% c("A", "B", "G", "N")) %>%
  
  # Add correct ref dates per line
  mutate(
    ref_date = case_when(
      RAIL_LINE == "A" ~ as.Date("2016-04-22"),
      RAIL_LINE == "B" ~ as.Date("2016-07-25"),
      RAIL_LINE == "G" ~ as.Date("2019-04-26"),
      RAIL_LINE == "N" ~ as.Date("2020-09-21"),
      TRUE ~ NA_Date_
    )
  )

# NNdist function for filtering/weights

rail.dist <- function(blckgrp, rail_stations, study_area) {
  # ---- CRS setup ----
  target_crs <- 26913  # UTM Zone 13N
  blckgrp <- st_transform(blckgrp, target_crs)
  rail_stations <- st_transform(rail_stations, target_crs)
  study_area <- st_transform(study_area, target_crs)
  
  # ---- Geometry prep ----
  center <- st_centroid(blckgrp)
  study_area_valid <- st_make_valid(st_union(study_area))
  rail_stations_valid <- st_make_valid(rail_stations)
  
  study_buffer <- st_buffer(study_area_valid, 200)
  win <- as.owin(st_union(study_buffer))
  
  # ---- Convert to planar points ----
  center.xy <- st_coordinates(center)
  rail.xy <- st_coordinates(rail_stations_valid)
  
  center.ppp <- ppp(center.xy[,1], center.xy[,2], window = win, check = FALSE)
  rail.ppp <- ppp(rail.xy[,1], rail.xy[,2], window = win, check = FALSE)
  
  # ---- Compute nearest station index and distance ----
  nn.df <- nncross(center.ppp, rail.ppp)  # returns both $which and $dist
  
  # ---- Attach route and distance ----
  rail_attrs <- st_drop_geometry(rail_stations_valid) %>%
    mutate(station_id = row_number())
  
  blckgrp.dist <- blckgrp %>%
    mutate(
      dist_to_rail = nn.df$dist,
      nearest_station = rail_attrs$NAME[nn.df$which],
      nearest_route = rail_attrs$RAIL_LINE[nn.df$which]
    )
  
  return(blckgrp.dist)
}

# Run nn dist function

denver.2010.dist <- rail.dist(blckgrp2010.denver,RTD.commuter.stations,denver2010)
denver.2020.dist <- rail.dist(blckgrp2020.denver,RTD.commuter.stations,denver2020)

# Sanity-check choropleth plot

ggplot(denver.2010.dist) +
  geom_sf(aes(fill = dist_to_rail), color = "white", size = 0.1) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  labs(
    title = "Denver Census Block Groups by Distance to Nearest Commuter Rail Station",
    fill = "Distance"
  ) +
  theme_minimal()
