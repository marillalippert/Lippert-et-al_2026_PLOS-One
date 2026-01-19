library(tidyverse)
library(sf)
library(ncdf4)
library(raster)
library(lubridate)
library(GGally)
library(ggpubr)
library(nlme)

# InSitu regional comp script needs to run before this script
# Get lat long limitations from in situ colonies

insitu_coords <- read.csv("Data/Nighttime_Full_Insitu.csv", stringsAsFactors = F)
insitu_coords$Latitude <- as.numeric(insitu_coords$Latitude)
insitu_coords$Longitude <- as.numeric(insitu_coords$Longitude)

insitu_colony_coords <- insitu_coords |> filter(!is.na(Latitude)) |> 
  group_by(Region, Colony) |> 
  summarize(
    Latitude = unique(Latitude), 
    Longitude = unique(Longitude))

insitu_coords_minmax <- insitu_colony_coords |>  group_by(Region) |> 
  summarize(min_lat = min(Latitude),
            min_long = min(Longitude),
            max_lat = max(Latitude),
            max_long = max(Longitude))

# SSST Download -----------
##url: https://pae-paha.pacioos.hawaii.edu/erddap/griddap/dhw_5km.htmlTable?CRW_SST%5B(2025-04-02T12:00:00Z):1:(2025-04-02T12:00:00Z)%5D%5B(89.975):1:(-89.975)%5D%5B(-179.975):1:(179.975)%5D
server <- "https://pae-paha.pacioos.hawaii.edu/erddap/griddap/"  ## server is the website
sst_id <- "dhw_5km.nc?CRW_SST" 
##10s are the southern regions, 20s are the western region, 30s are northern
latmin <- list(min10 = format(floor(100 * as.numeric(insitu_coords_minmax[2,2])) / 100, nsmall = 2), 
               min20 = format(floor(100 * as.numeric(insitu_coords_minmax[3,2])) / 100, nsmall = 2), 
               min30 = format(floor(100 * as.numeric(insitu_coords_minmax[1,2])) / 100, nsmall = 2))
latmax <- list(max10 = format(ceiling(100 * as.numeric(insitu_coords_minmax[2, 4])) / 100, nsmall = 2), 
               max20 = format(ceiling(100 * as.numeric(insitu_coords_minmax[3,4])) / 100, nsmall = 2), 
               max30 = format(ceiling(100 * as.numeric(insitu_coords_minmax[1,4])) / 100, nsmall = 2))
longmin <- list(min10 = format(floor(100 * as.numeric(insitu_coords_minmax[2, 3])) / 100, nsmall = 2), 
                min20 = format(floor(100 * as.numeric(insitu_coords_minmax[3,3])) / 100, nsmall = 2), 
                min30 = format(floor(100 * as.numeric(insitu_coords_minmax[1,3])) / 100, nsmall = 2))
longmax <-list(max10 = format(ceiling(100 * as.numeric(insitu_coords_minmax[2, 5])) / 100, nsmall = 2), 
               max20 = format(ceiling(100 * as.numeric(insitu_coords_minmax[3,5])) / 100, nsmall = 2), 
               max30 = format(ceiling(100 * as.numeric(insitu_coords_minmax[1,5])) / 100, nsmall = 2))
start_date <- as.Date("2017-01-01")
end_date <- as.Date("2021-04-01")

site <- c("CRW_Southern_Palau_url", "CRW_Western_Palau_url", "CRW_Northern_Palau_url")

for(i in 1:length(latmin)) {
  site[i] <- paste0(server, sst_id, "%5B(", start_date, "T09:00:00Z):1:(", 
                    end_date, "T09:00:00Z)%5D%5B(", latmin[i], "):1:(", latmax[i], ")%5D%5B(", 
                    longmin[i], "):1:(", longmax[i], ")%5D")
}


if(!file.exists("Data/CRW_Southern_Palau_sst.nc")) curl::curl_download(site[1], "Data/CRW_Southern_Palau_sst.nc")
sst_southern <- raster::stack(paste0(getwd(), "/Data/CRW_Southern_Palau_sst.nc"))
#sst_southern <- terra::rast(paste0(getwd(), "/Data/CRW_Southern_Palau_sst.nc"))
#sst_southern <- raster(paste0(getwd(), "/Data/CRW_Southern_Palau_sst.nc"))
sst_southern <- as.data.frame(sst_southern, xy = TRUE, long = TRUE)
sst_southern$site <- c("Southern")

if(!file.exists("Data/CRW_Western_Palau_sst.nc")) curl::curl_download(site[2], "Data/CRW_Western_Palau_sst.nc")
sst_western <- raster::stack(paste0(getwd(), "/Data/CRW_Western_Palau_sst.nc"))
sst_western <- as.data.frame(sst_western, xy = TRUE, long = TRUE)
sst_western$site <- c("Western")

if(!file.exists("Data/CRW_Northern_Palau_sst.nc")) curl::curl_download(site[3], "Data/CRW_Northern_Palau_sst.nc")
sst_northern <- raster::stack(paste0(getwd(), "/Data/CRW_Northern_Palau_sst.nc")) 
sst_northern <- as.data.frame(sst_northern, xy = TRUE, long = TRUE)
sst_northern$site <- c("Northern")

sst <- rbind(sst_southern, sst_western, sst_northern)
write.csv(sst, "Data/Regional_CRW_SST_Full.csv", row.names = F)


