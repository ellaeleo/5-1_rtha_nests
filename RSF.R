#### INSTALL/ LOAD PACKAGES ####
#install.packages("dplyr")
library(dplyr)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("sf")
library(sf)
#install.packages("terra")
library(terra)
#install.packages("janitor")
library(janitor)
#install.packages("naturalearth")
library(rnaturalearth)
#install.packages("ggspatial")
library(ggspatial)
##install.packages("amt")
library(amt)
##install.packages("lme4")
library(lme4)
##install.packages("ape")
library(ape)
##install.packages("units")
library(units)
##install.packages("gstat")
library(gstat)

sf_use_s2(FALSE)



##### DATA TO COORDINATES #####
## Included in repo!!!
rtha_coord <-read.csv("dat_RTHA_LA.csv") |>
  ## Convert this data to include coordinates based on "Latitude" and "Longitude"
  st_as_sf(coords=c("Longitude", "Latitude"), 
           ## This number is a shortcut to indicate WGS84, the datum that we are 
           ## using to project our data
           crs = 4326)
## View the coordinate data
st_crs(rtha_coord)



##### FILTER FOR LA SHAPEFILE #####
## Included in repo!!!
la_city_shape <- read_sf(dsn = ".",
                         ## layer = is the name of the files. Shapefiles use multiple files, which
                         ## why multiple files were downloaded in the zip. R will read all of these
                         layer = "City_Boundary")
### st_intersection clips the input data (first input) based on a shape (second input)
rtha_la <- st_intersection(rtha_coord, la_city_shape)



##### PLOTTING DATA #####
ggplot(data = rtha_la) +
  geom_sf()
### To meters
rtha_la <- rtha_la |>
  st_transform(crs = st_crs(26910))
ggplot(data = rtha_la) + geom_sf() +
  coord_sf(datum = st_crs(26910))



##### PLOT RANDOM POINTS #####
la_x_y <- as_tibble(st_coordinates(rtha_la))
rtha_la <- cbind(rtha_la, la_x_y)
## Make track
rtha_tk <- rtha_la |> 
  mutate(ID = seq.int(nrow(rtha_la))) |> 
  make_track(.x = X,
             .y = Y,
             .t = Year,
             id = ID,
             crs = st_crs(26910))
## Random points
r_pts <- rtha_tk |> 
  random_points(n = 71)

plot(r_pts)



##### LOAD RASTER #####
### Included in repo!!!
la_landuse <- rast("~/R_projects/5-1_rtha_nests/la_landuse_3.tif")

la_landuse_proj <- project(la_landuse, "+proj=utm +zone=10")
plot(la_landuse_proj)



#### CREATE BUFFER WITH TRUE AND REPLICATED POINTS ####
rsf1_geom <- st_as_sf(
  r_pts, coords = c("x_", "y_"), crs = 26910)

buffer_rtha_r <- st_buffer(rsf1_geom, dist = 2.41)



#### EXTRACT COVARIATES FROM LAND USE ####
buff_extract <- extract(la_landuse_proj, buffer_rtha_r, weights = F) |>
  ## Label true and replicated points
  mutate(case_ = ifelse(ID <= 71, 0, 1)) |>
  ## Create a categorical variable based off of the bands
  mutate(Land_class = ifelse(Band_1 >= 11 
                             & Band_1 < 12, 
                             "Water",
                             ifelse(Band_1 >= 12 
                                    & Band_1 < 21,
                                    "Snow",
                                    ifelse(Band_1 >= 21 
                                           & Band_1 < 22,
                                           "D_OS",
                                           ifelse(Band_1 >= 22 
                                                  & Band_1 < 23,
                                                  "D_LS",
                                                  ifelse(
                                                    Band_1 >= 23 
                                                    & Band_1 < 24,
                                                    "D_MS",
                                                    ifelse(Band_1 >= 24 
                                                           & Band_1 < 31,
                                                           "D_HS",
                                                           ifelse(Band_1 >= 31 
                                                                  & Band_1 < 81
                                                                  | Band_1 >= 90,
                                                                  "Natural",
                                                                  ifelse(Band_1 >= 81 
                                                                         & Band_1 < 90,
                                                                         "Ag", NA)))))))))




#### FIT RESOURCE SELECTION FUNCTION ####
buff_extract |> 
  fit_rsf(formula = case_ ~ Land_class) |> 
  summary()



#### TEST FOR SPATIAL AUTOCORRELATION ####

rtha_sa <- r_pts |> 
  extract_covariates(la_landuse_proj)

rtha_sa_geom <- st_as_sf(
  rtha_sa, coords = c("x_", "y_"), crs = 26910)

variogram(Band_1 ~1, data = rtha_sa_geom) |> plot()

mod <- lm(Band_1 ~ case_,
          data = rtha_sa_geom)

variogram(residuals(mod) ~ 1, rtha_sa_geom) |> plot()

inverse <- function(x) {
  1/x
}
distmat <- rsf1_geom |> 
  st_distance() |> 
  drop_units()
distmat[distmat == 0] <- 1
invdist <- distmat |> inverse()

rtha_sa_geom <- rtha_sa_geom |> 
  mutate(resid = mod |> 
           residuals())

Moran.I(rtha_sa_geom$resid, invdist)

