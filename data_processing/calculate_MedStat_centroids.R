################################################################################
# Centroid of population density per MedStat region
# - load the data population densities and shapefile of Medstat regions
# - cut out the pop densities per shapefile
# - calculate the centroid per shapefile




### DATA ####
rm(list=ls())

# packages
library(raster); library(sf); library(ggplot2); library(knitr); library(RColorBrewer); library(exactextractr); library(dplyr)


# load population densities of Switzerland
pop <- raster("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data-raw/population/data/gpw_v4_population_density_rev11_2010_30_sec_3.asc")

# load spahefile data
mesh <- st_read("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data-raw/MedStat/MEDSTAT_AREAS_2019.shp", quiet = TRUE)

#####


### CROP TO SWITZERLAND ####

# longitudes 5 to 12 and latitudes 45 to 48
crop_extent <- extent(5, 12, 45, 48)

pop_sui <- crop(pop, crop_extent)

## plot
# max min values
zlim_values <- c(0, 5000)
# number of colors
num_colors <- 5
# number of legend intervals
custom_breaks <- seq(0, 5000, by = 1000)
# color pallete
colors <- colorRampPalette(brewer.pal(5, "Greys"))(num_colors)


plot(pop_sui,
     col = colors,
     zlim = zlim_values, ylim = c(45,48), xlim = c(5,12),
     main = "Population density of Switzerland",
     breaks = custom_breaks)

#####


### CUT OUT SF FROM RASTER ####

# Check CRS of both raster and shapefile
crs_raster <- crs(pop_sui)
crs_shapefile <- st_crs(mesh$geometry)

# Reproject the shapefile to match the raster's CRS if necessary
if (crs_raster != crs_shapefile) {
  mesh$geometry <- st_transform(mesh$geometry, crs = crs_raster)
}


# Function to calculate weighted centroid
calculate_weighted_centroid <- function(extract_data) {
  # Extract the population densities and coordinates (centroids of raster cells)
  densities <- extract_data$value
  x_coords <- extract_data$x
  y_coords <- extract_data$y

  # Calculate the weighted x and y centroids
  weighted_x <- sum(x_coords * densities, na.rm = TRUE) / sum(densities, na.rm = TRUE)
  weighted_y <- sum(y_coords * densities, na.rm = TRUE) / sum(densities, na.rm = TRUE)

  return(c(weighted_x, weighted_y))
}

# Extract both raster values and coordinates for each polygon
centroids <- exact_extract(pop_sui, mesh$geometry, fun = calculate_weighted_centroid, include_xy = TRUE, summarize_df = TRUE)
centroids <- t(as.data.frame(centroids))

# coordinates to mesh file
mesh$centroid_x <- centroids[,1]
mesh$centroid_y <- centroids[,2]


# Plot the raster
plot(pop_sui,
     col = colors,
     zlim = zlim_values, ylim = c(45.5,48), xlim = c(5,12),
     main = "Population density of Switzerland",
     breaks = custom_breaks)

# Add the shapefile polygons to the plot
plot(st_geometry(mesh$geometry), add = TRUE, border = "red", lwd = .2)

# Plot the weighted centroids
points(mesh$centroid_x, mesh$centroid_y, col = "blue", pch = 19, cex = .3)
#####

## WRITE CSV ####

st_write(mesh, "C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data-raw/MedStat/MedStat_Areas_2019_w_centroids_popdensity.shp")
######
