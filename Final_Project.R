library(sf)
library(spatstat)
library(spdep)
library(spData)
library(tidyverse)
library(lubridate)
library(rlang)
library(fixest)
library(modelsummary)

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

blckgrp2010.denver <- st_intersection(blckgrp2010, denver2010) %>%
  mutate(GEOID10 = as.numeric(GEOID10))


denver2020 <- st_transform(denver2020, st_crs(blckgrp2020)) %>%
  st_union() %>%
  st_make_valid()

blckgrp2020.denver <- st_intersection(blckgrp2020, denver2020) %>%
  mutate(GEOID10 = as.numeric(GEOID))

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
  
  # Add correct ref dates per line
  mutate(
    ref_date = case_when(
      str_detect(RAIL_LINE, "A") ~ as.Date("2016-04-22"),
      str_detect(RAIL_LINE, "N") & !str_detect(RAIL_LINE, "A") ~ as.Date("2020-09-21"),
      PID == 228 ~ as.Date("2016-07-25"),
      str_detect(RAIL_LINE, "G") & !str_detect(RAIL_LINE, "A") ~ as.Date("2019-04-26"),
      TRUE ~ NA_Date_
    )
  )

# Create buffered treatment and control groups for event study

control.treatment <- function(blckgrp, rail_stations, study_area, data, lehd, year, stdevs=2) {
  # ---- CRS setup ----
  target_crs <- 26913  # UTM Zone 13N
  blckgrp <- st_transform(blckgrp, target_crs)
  rail_stations <- st_transform(rail_stations, target_crs)
  study_area <- st_transform(study_area, target_crs)
  
  # ---- Filter rail stations by ref date ----
  rail_stations.filtered <- rail_stations %>%
    filter(year(ref_date) <= year + 3)
  
  # ---- Geometry prep ----
  center <- st_centroid(blckgrp)
  study_area_valid <- st_make_valid(st_union(study_area))
  rail_stations_valid <- st_make_valid(rail_stations.filtered)
  rail_buffer <- st_buffer(rail_stations_valid, 1180) # buffer 1 standard deviation
  
  study_buffer <- st_buffer(study_area_valid, 200)
  win <- as.owin(st_union(study_buffer))
  
  # ---- Convert to planar points ----
  center.xy <- st_coordinates(center)
  rail.xy <- st_coordinates(rail_stations_valid)
  
  center.ppp <- ppp(center.xy[,1], center.xy[,2], window = win, check = FALSE)
  rail.ppp <- ppp(rail.xy[,1], rail.xy[,2], window = win, check = FALSE)
  
  # ---- Compute nearest station index and distance ----
  nn.df <- nncross(center.ppp, rail.ppp)  # returns $which and $dist
  
  # ---- Attach route, distance, and ref date ----
  rail_attrs <- st_drop_geometry(rail_stations_valid) %>%
    mutate(station_id = row_number())
  
  blckgrp.dist <- blckgrp %>%
    mutate(
      center_dist_to_rail = nn.df$dist,
      nearest_station = rail_attrs$NAME[nn.df$which],
      nearest_route = rail_attrs$RAIL_LINE[nn.df$which],
      ref_date = rail_attrs$ref_date[nn.df$which]
    )
  
  # ---- Join and return data set, assigning control and treatment groups ----
  blckgrp.data <- inner_join(blckgrp.dist, data, by = c("GISJOIN" = "GISJOIN")) %>%
    mutate(
      treat = case_when(
        center_dist_to_rail <= 805 + 375 * stdevs | lengths(st_intersects(blckgrp.dist, rail_buffer)) > 0 ~ 1,
        TRUE ~ 0
      ))
  
  blckgrp.lehd.data <- left_join(blckgrp.data,lehd, by = c("GEOID10" = "id"))
  
  return(blckgrp.lehd.data)
}

# Load ACS data
data.2013 <- read.csv("nhgis0011_csv/nhgis0011_ds201_20135_blck_grp.csv")
data.2014 <- read.csv("nhgis0011_csv/nhgis0011_ds206_20145_blck_grp.csv")
data.2015 <- read.csv("nhgis0011_csv/nhgis0011_ds215_20155_blck_grp.csv")
data.2016 <- read.csv("nhgis0011_csv/nhgis0011_ds225_20165_blck_grp.csv")
data.2017 <- read.csv("nhgis0011_csv/nhgis0011_ds233_20175_blck_grp.csv")
data.2018 <- read.csv("nhgis0011_csv/nhgis0011_ds239_20185_blck_grp.csv")
data.2019 <- read.csv("nhgis0011_csv/nhgis0011_ds244_20195_blck_grp.csv")
data.2020 <- read.csv("nhgis0011_csv/nhgis0011_ds249_20205_blck_grp.csv")
data.2021 <- read.csv("nhgis0011_csv/nhgis0011_ds254_20215_blck_grp.csv")
data.2022 <- read.csv("nhgis0011_csv/nhgis0011_ds262_20225_blck_grp.csv")

# Load LEHD data

lehd.2013 <- read.csv("polygon_2013_edit.csv")
lehd.2014 <- read.csv("polygon_2014_edit.csv")
lehd.2015 <- read.csv("polygon_2015_edit.csv")
lehd.2016 <- read.csv("polygon_2016_edit.csv")
lehd.2017 <- read.csv("polygon_2017_edit.csv")
lehd.2018 <- read.csv("polygon_2018_edit.csv")
lehd.2019 <- read.csv("polygon_2019_edit.csv")
lehd.2020 <- read.csv("polygon_2020_edit.csv")
lehd.2021 <- read.csv("polygon_2021_edit.csv")
lehd.2022 <- read.csv("polygon_2022_edit.csv")

# Create study datasets

denver.study.2013 <- control.treatment(blckgrp2010.denver,RTD.commuter.stations,denver2010,data.2013,lehd.2013,2013)
denver.study.2014 <- control.treatment(blckgrp2010.denver,RTD.commuter.stations,denver2010,data.2014,lehd.2014,2014)
denver.study.2015 <- control.treatment(blckgrp2010.denver,RTD.commuter.stations,denver2010,data.2015,lehd.2015,2015)
denver.study.2016 <- control.treatment(blckgrp2010.denver,RTD.commuter.stations,denver2010,data.2016,lehd.2016,2016)
denver.study.2017 <- control.treatment(blckgrp2010.denver,RTD.commuter.stations,denver2010,data.2017,lehd.2017,2017)
denver.study.2018 <- control.treatment(blckgrp2010.denver,RTD.commuter.stations,denver2010,data.2018,lehd.2018,2018)
denver.study.2019 <- control.treatment(blckgrp2010.denver,RTD.commuter.stations,denver2010,data.2019,lehd.2019,2019)
denver.study.2020 <- control.treatment(blckgrp2020.denver,RTD.commuter.stations,denver2020,data.2020,lehd.2020,2020)
denver.study.2021 <- control.treatment(blckgrp2020.denver,RTD.commuter.stations,denver2020,data.2021,lehd.2021,2021)
denver.study.2022 <- control.treatment(blckgrp2020.denver,RTD.commuter.stations,denver2020,data.2022,lehd.2022,2022)

# Expected value functions for discrete values of rent and travel times

gmr <- function(data, rent, expand = TRUE) {
  if (expand == FALSE) {
    values <- c(50,124.5,174.5,224.5,274.5,324.5,374.5,424.5,474.5,524.5,
                574.5,624.5,674.5,724.5,774.5,849.5,949.5,1124.5,1374.5,
                1749.5,2750)
    suffixes <- sprintf("%03d", 3:23)
  } else {
    values <- c(50,124.5,174.5,224.5,274.5,324.5,374.5,424.5,474.5,524.5,
                574.5,624.5,674.5,724.5,774.5,849.5,949.5,1124.5,1374.5,
                1749.5,2249.5,2749.5,3249.5,3749.5)
    suffixes <- sprintf("%03d", 3:26)
  }
  
  counts <- sapply(paste0(rent, suffixes), function(col) data[[col]])
  total <- data[[paste0(rent, "002")]]
  probs <- counts / total
  sum(values * probs, na.rm = TRUE)
}

travel <- function(data, traveltime) {
  values <- c(2.5,7,12,17,22,27,32,37,42,57,74.5,100)
  suffixes <- sprintf("%03d", 2:13)
  counts <- sapply(paste0(traveltime, suffixes), function(col) data[[col]])
  total <- data[[paste0(traveltime, "001")]]
  probs <- counts / total
  sum(values * probs, na.rm = TRUE)
}

# Create regression variables

denver.regression.vars.2013 <- denver.study.2013 %>%
  mutate(UHDE001 = as.numeric(UHDE001), UMME001 = as.numeric(UMME001)) %>% # Make sure mhi and medvalhu are numeric
  mutate(lpop = log(UEPE001), lpdensity = log(UEPE001/ALAND10), pct_transit = UFFE010/UFHE001, 
         lmhi = log(UHDE001), pct_unem = UJ8E005/UJ8E003, lhousing = log(UKNE001), 
         lmedvalhu = log(UMME001),ljobs = log(c000), 
         ljobdensity = log(c000/ALAND10), lhpjobs = log(ce03)) %>% # Create regression vars
  rowwise() %>%
  mutate(lrent = log(gmr(cur_data(), "UL8E", expand = FALSE)), 
         ltravel = log(travel(cur_data(), "UFHE"))) %>% # Run expected value functions for rent & travel time
  ungroup() %>%
  select(GISJOIN,center_dist_to_rail,nearest_station,nearest_route,ref_date,YEAR,treat,
         lpop,lpdensity,pct_transit,ltravel,lmhi,lrent,pct_unem,lhousing,lmedvalhu,ljobs,
         ljobdensity,lhpjobs,INTPTLAT10,INTPTLON10,Shape_area,
         Shape_len,geometry) # Select interest variables and fixed effects

denver.regression.vars.2014 <- denver.study.2014 %>%
  mutate(ABDPE001 = as.numeric(ABDPE001), ABITE001 = as.numeric(ABITE001)) %>%
  mutate(lpop = log(ABA1E001), lpdensity = log(ABA1E001/ALAND10), pct_transit = ABBRE011/ABBTE001, 
         lmhi = log(ABDPE001), pct_unem = ABGFE005/ABGFE003, lhousing = log(ABGWE001), 
         lmedvalhu = log(ABITE001),ljobs = log(c000), ljobdensity = log(c000/ALAND10), 
         lhpjobs = log(ce03)) %>%
  rowwise() %>%
  mutate(lrent = log(gmr(cur_data(), "ABIGE", expand = FALSE)), 
         ltravel = log(travel(cur_data(), "ABBTE"))) %>%
  ungroup() %>%
  select(GISJOIN,center_dist_to_rail,nearest_station,nearest_route,ref_date,YEAR,treat,
         lpop,lpdensity,pct_transit,ltravel,lmhi,lrent,pct_unem,lhousing,lmedvalhu,ljobs,
         ljobdensity,lhpjobs,INTPTLAT10,INTPTLON10,Shape_area,
         Shape_len,geometry) # Select interest variables and fixed effects

denver.regression.vars.2015 <- denver.study.2015 %>%
  mutate(ADNKE001 = as.numeric(ADNKE001), ADRWE001 = as.numeric(ADRWE001)) %>%
  mutate(lpop = log(ADKWE001), lpdensity = log(ADKWE001/ALAND10), pct_transit = ADLME010/ADLOE001, 
         lmhi = log(ADNKE001), pct_unem = ADPIE005/ADPIE003, lhousing = log(ADPZE001), 
         lmedvalhu = log(ADRWE001),ljobs = log(c000), ljobdensity = log(c000/ALAND10), 
         lhpjobs = log(ce03)) %>%
  rowwise() %>%
  mutate(lrent = log(gmr(cur_data(), "ADRJE")), 
         ltravel = log(travel(cur_data(), "ADLOE"))) %>%
  ungroup() %>%
  select(GISJOIN,center_dist_to_rail,nearest_station,nearest_route,ref_date,YEAR,treat,
         lpop,lpdensity,pct_transit,ltravel,lmhi,lrent,pct_unem,lhousing,lmedvalhu,ljobs,
         ljobdensity,lhpjobs,INTPTLAT10,INTPTLON10,Shape_area,
         Shape_len,geometry) # Select interest variables and fixed effects

denver.regression.vars.2016 <- denver.study.2016 %>%
  mutate(AF49E001 = as.numeric(AF49E001), AF9LE001 = as.numeric(AF9LE001)) %>%
  mutate(lpop = log(AF2LE001), lpdensity = log(AF2LE001/ALAND10), pct_transit = AF3BE010/AF3DE001, 
         lmhi = log(AF49E001), pct_unem = AF67E005/AF67E003, lhousing = log(AF7OE001), 
         lmedvalhu = log(AF9LE001),ljobs = log(c000), ljobdensity = log(c000/ALAND10), 
         lhpjobs = log(ce03)) %>%
  rowwise() %>%
  mutate(lrent = log(gmr(cur_data(), "AF88E")), 
         ltravel = log(travel(cur_data(), "AF3DE"))) %>%
  ungroup() %>%
  select(GISJOIN,center_dist_to_rail,nearest_station,nearest_route,ref_date,YEAR,treat,
         lpop,lpdensity,pct_transit,ltravel,lmhi,lrent,pct_unem,lhousing,lmedvalhu,ljobs,
         ljobdensity,lhpjobs,INTPTLAT10,INTPTLON10,Shape_area,
         Shape_len,geometry) # Select interest variables and fixed effects

denver.regression.vars.2017 <- denver.study.2017 %>%
  mutate(AH1PE001 = as.numeric(AH1PE001), AH53E001 = as.numeric(AH53E001)) %>%
  mutate(lpop = log(AHY1E001), lpdensity = log(AHY1E001/ALAND10), pct_transit = AHZRE010/AHZTE001, 
         lmhi = log(AH1PE001), pct_unem = AH3PE005/AH3PE003, lhousing = log(AH36E001), 
         lmedvalhu = log(AH53E001),ljobs = log(c000), ljobdensity = log(c000/ALAND10), 
         lhpjobs = log(ce03)) %>%
  rowwise() %>%
  mutate(lrent = log(gmr(cur_data(), "AH5QE")), 
         ltravel = log(travel(cur_data(), "AHZTE"))) %>%
  ungroup() %>%
  select(GISJOIN,center_dist_to_rail,nearest_station,nearest_route,ref_date,YEAR,treat,
         lpop,lpdensity,pct_transit,ltravel,lmhi,lrent,pct_unem,lhousing,lmedvalhu,ljobs,
         ljobdensity,lhpjobs,INTPTLAT10,INTPTLON10,Shape_area,
         Shape_len,geometry) # Select interest variables and fixed effects

denver.regression.vars.2018 <- denver.study.2018 %>%
  mutate(AJZAE001 = as.numeric(AJZAE001), AJ3QE001 = as.numeric(AJ3QE001)) %>%
  mutate(lpop = log(AJWME001), lpdensity = log(AJWME001/ALAND10), pct_transit = AJXCE010/AJXEE001, 
         lmhi = log(AJZAE001), pct_unem = AJ1CE005/AJ1CE003, lhousing = log(AJ1TE001), 
         lmedvalhu = log(AJ3QE001),ljobs = log(c000), ljobdensity = log(c000/ALAND10), 
         lhpjobs = log(ce03)) %>%
  rowwise() %>%
  mutate(lrent = log(gmr(cur_data(), "AJ3DE")), 
         ltravel = log(travel(cur_data(), "AJXEE"))) %>%
  ungroup() %>%
  select(GISJOIN,center_dist_to_rail,nearest_station,nearest_route,ref_date,YEAR,treat,
         lpop,lpdensity,pct_transit,ltravel,lmhi,lrent,pct_unem,lhousing,lmedvalhu,ljobs,
         ljobdensity,lhpjobs,INTPTLAT10,INTPTLON10,Shape_area,
         Shape_len,geometry) # Select interest variables and fixed effects

denver.regression.vars.2019 <- denver.study.2019 %>%
  mutate(ALW1E001 = as.numeric(ALW1E001), AL1HE001 = as.numeric(AL1HE001)) %>%
  mutate(lpop = log(ALUBE001), lpdensity = log(ALUBE001/ALAND10), pct_transit = ALU1E010/ALU3E001, 
         lmhi = log(ALW1E001), pct_unem = ALY3E005/ALY3E003, lhousing = log(ALZKE001), 
         lmedvalhu = log(AL1HE001),ljobs = log(c000), ljobdensity = log(c000/ALAND10), 
         lhpjobs = log(ce03)) %>%
  rowwise() %>%
  mutate(lrent = log(gmr(cur_data(), "AL04E")), 
         ltravel = log(travel(cur_data(), "ALU3E"))) %>%
  ungroup() %>%
  select(GISJOIN,center_dist_to_rail,nearest_station,nearest_route,ref_date,YEAR,treat,
         lpop,lpdensity,pct_transit,ltravel,lmhi,lrent,pct_unem,lhousing,lmedvalhu,ljobs,
         ljobdensity,lhpjobs,INTPTLAT10,INTPTLON10,Shape_area,
         Shape_len,geometry) # Select interest variables and fixed effects

denver.regression.vars.2020 <- denver.study.2020 %>%
  mutate(AMR8E001 = as.numeric(AMR8E001), AMWBE001 = as.numeric(AMWBE001)) %>%
  mutate(lpop = log(AMPVE001), lpdensity = log(AMPVE001/ALAND), pct_transit = AMQKE010/AMQME001, 
         lmhi = log(AMR8E001), pct_unem = AMT9E005/AMT9E003, lhousing = log(AMUEE001), 
         lmedvalhu = log(AMWBE001),ljobs = log(c000), ljobdensity = log(c000/ALAND), 
         lhpjobs = log(ce03)) %>%
  rowwise() %>%
  mutate(lrent = log(gmr(cur_data(), "AMVYE")), 
         ltravel = log(travel(cur_data(), "AMQME"))) %>%
  ungroup() %>%
  select(GISJOIN,center_dist_to_rail,nearest_station,nearest_route,ref_date,YEAR,treat,
         lpop,lpdensity,pct_transit,ltravel,lmhi,lrent,pct_unem,lhousing,lmedvalhu,ljobs,
         ljobdensity,lhpjobs,INTPTLAT,INTPTLON,Shape_Area,
         Shape_Leng,geometry) # Select interest variables and fixed effects

denver.regression.vars.2021 <- denver.study.2021 %>%
  mutate(AOQIE001 = as.numeric(AOQIE001), AOULE001 = as.numeric(AOULE001)) %>%
  mutate(lpop = log(AON4E001), lpdensity = log(AON4E001/ALAND), pct_transit = AOOTE010/AOOVE001, 
         lmhi = log(AOQIE001), pct_unem = AOSJE005/AOSJE003, lhousing = log(AOSOE001), 
         lmedvalhu = log(AOULE001),ljobs = log(c000), ljobdensity = log(c000/ALAND), 
         lhpjobs = log(ce03)) %>%
  rowwise() %>%
  mutate(lrent = log(gmr(cur_data(), "AOT8E")), 
         ltravel = log(travel(cur_data(), "AOOVE"))) %>%
  ungroup() %>%
  select(GISJOIN,center_dist_to_rail,nearest_station,nearest_route,ref_date,YEAR,treat,
         lpop,lpdensity,pct_transit,ltravel,lmhi,lrent,pct_unem,lhousing,lmedvalhu,ljobs,
         ljobdensity,lhpjobs,INTPTLAT,INTPTLON,Shape_Area,
         Shape_Leng,geometry) # Select interest variables and fixed effects

denver.regression.vars.2022 <- denver.study.2022 %>%
  mutate(AQP6E001 = as.numeric(AQP6E001), AQU4E001 = as.numeric(AQU4E001)) %>%
  mutate(lpop = log(AQNFE001), lpdensity = log(AQNFE001/ALAND), pct_transit = AQN5E010/AQN7E001, 
         lmhi = log(AQP6E001), pct_unem = AQR8E005/AQR8E003, lhousing = log(AQSPE001), 
         lmedvalhu = log(AQU4E001),ljobs = log(c000), ljobdensity = log(c000/ALAND), 
         lhpjobs = log(ce03)) %>%
  rowwise() %>%
  mutate(lrent = log(gmr(cur_data(), "AQURE")), 
         ltravel = log(travel(cur_data(), "AQN7E"))) %>%
  ungroup() %>%
  select(GISJOIN,center_dist_to_rail,nearest_station,nearest_route,ref_date,YEAR,treat,
         lpop,lpdensity,pct_transit,ltravel,lmhi,lrent,pct_unem,lhousing,lmedvalhu,ljobs,
         ljobdensity,lhpjobs,INTPTLAT,INTPTLON,Shape_Area,
         Shape_Leng,geometry) # Select interest variables and fixed effects

# Bind rows

denver.regression.set <- bind_rows(denver.regression.vars.2013,denver.regression.vars.2014,
                                   denver.regression.vars.2015,denver.regression.vars.2016,
                                   denver.regression.vars.2017,denver.regression.vars.2018,
                                   denver.regression.vars.2019,denver.regression.vars.2020,
                                   denver.regression.vars.2021,denver.regression.vars.2022) %>%
  mutate(YEAR = as.numeric(strtrim(YEAR, 4)) + 4, event_time = YEAR - year(ref_date))

# Regress the interest variables

population <- feols(lpop ~ i(event_time, treat, ref = -1) | GISJOIN + YEAR, data = denver.regression.set)
density <- feols(lpdensity ~ i(event_time, treat, ref = -1) | GISJOIN + YEAR, data = denver.regression.set)
transit <- feols(pct_transit ~ i(event_time, treat, ref = -1) | GISJOIN + YEAR, data = denver.regression.set)
commute_time <- feols(ltravel ~ i(event_time, treat, ref = -1) | GISJOIN + YEAR, data = denver.regression.set)
mhi <- feols(lmhi ~ i(event_time, treat, ref = -1) | GISJOIN + YEAR, data = denver.regression.set)
rent<- feols(lrent ~ i(event_time, treat, ref = -1) | GISJOIN + YEAR, data = denver.regression.set)
unemployment <- feols(pct_unem ~ i(event_time, treat, ref = -1) | GISJOIN + YEAR, data = denver.regression.set)
housing <- feols(lhousing ~ i(event_time, treat, ref = -1) | GISJOIN + YEAR, data = denver.regression.set)
medvalhu <- feols(lmedvalhu ~ i(event_time, treat, ref = -1) | GISJOIN + YEAR, data = denver.regression.set)
jobs <- feols(ljobs ~ i(event_time, treat, ref = -1) | GISJOIN + YEAR, data = denver.regression.set)
high_paying_jobs <- feols(lhpjobs ~ i(event_time, treat, ref = -1) | GISJOIN + YEAR, data = denver.regression.set)

coefplot(population)
coefplot(density)
coefplot(transit)
coefplot(commute_time)
coefplot(mhi)
coefplot(rent)
coefplot(unemployment)
coefplot(housing)
coefplot(medvalhu)
coefplot(jobs)
coefplot(high_paying_jobs)

summary(population)
summary(density)
summary(transit)
summary(commute_time)
summary(mhi)
summary(rent)
summary(unemployment)
summary(housing)
summary(medvalhu)
summary(jobs)
summary(high_paying_jobs)

models <- list(population,transit,commute_time,mhi,rent,unemployment,housing,medvalhu,
               jobs,high_paying_jobs)
names(models) <- c("Log Populataion","Percent of Commuters who Use Transit","Log Average Commute Time",
           "Log Median Household Income","Log Gross Mean Rent","Percent of Unemployment",
           "Log Total Housing Units","Log Median Value of a Housing Unit","Log Total Jobs",
           "Log Jobs over $3,333 in Monthly Earnings")
modelsummary(models, stars = TRUE)
