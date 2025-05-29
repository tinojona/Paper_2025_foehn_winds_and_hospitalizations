################################################################################
# Centroid of population density per MedStat region
################################################################################
# - load the data population densities and shapefile of Medstat regions
# - cut out the pop densities per shapefile
# - calculate the centroid per shapefile



### DATA
#-------

rm(list=ls())

# packages
library(raster); library(sf); library(ggplot2); library(knitr); library(RColorBrewer); library(exactextractr); library(dplyr);library(terra)


# load population densities of Switzerland
pop <- raster("../data-raw/SUI_population_densities/gpw_v4_population_density_rev11_2010_30_sec_3.asc")

# load spahefile data
mesh <- st_read("data_nonsensitive/MedStat_shapefiles/raw/MEDSTAT_AREAS_2019.shp", quiet = TRUE)

#--------


### CROP TO SWITZERLAND

#--------
pop_sui <- crop(pop, extent(5, 12, 45, 48))

## plot switzerlands population density
# plot(pop_sui,
#      col = colorRampPalette(brewer.pal(5, "Greys"))(5),
#      zlim = c(0, 400),
#      breaks = c(0, 50, 100, 200, 300, 400),
#      asp = 1,
#      main = "Population density of Switzerland")

#-------


### CUT OUT SF FROM RASTER
#-------

# Reproject coordinate reference system of shapefiles
mesh$geometry <- st_transform(mesh$geometry, crs = crs(pop_sui))



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
centroids <- as.data.frame(t(centroids))
colnames(centroids) = c("x","y")

# # Plot the raster
# plot(pop_sui,
#      col = colorRampPalette(brewer.pal(5, "Greys"))(5),
#      zlim =  c(0, 400),
#      asp = 1,
#      main = "Population density of Switzerland",
#      breaks =  c(0, 50, 100, 200, 300, 400))
#
# # Add the shapefile polygons to the plot
# plot(st_geometry(mesh$geometry), add = TRUE, border = "red", lwd = .4)
#
# # Plot the weighted centroids
# points(centroids$x, centroids$y, col = "blue", pch = 19, cex = .4)

#-----



### CONVERT COORDINATES TO SHAPEFILE POINTS ####
#-----

# add columns
centroids$MDSTID <-  mesh$MDSTID
centroids$MDST04 <-  mesh$MDST04

# Convert the data frame to an sf object (geometry points)
sf_centroids <- st_as_sf(centroids,
                         coords = c("x", "y"),
                         crs = crs(pop_sui))

# convert mesh and centroids crs to Meteoschweiz
mesh$geometry <- st_transform(mesh$geometry, crs= 2056)
sf_centroids$geometry <- st_transform(sf_centroids$geometry, crs = 2056)

#-----


### WRITE CSV
#----

# new mesh (crs adjusted)
st_write(mesh, "data_nonsensitive/Medstat_shapefiles/processed/MedStat_csr_adjusted.shp")

# centroids
st_write(sf_centroids, "data_nonsensitive/Medstat_shapefiles/centroids/MedStat_centroids_popdensity.shp")

#----
