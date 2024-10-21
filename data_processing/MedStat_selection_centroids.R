################################################################################
# Determine the MedStat regions for different buffers that overlap with the centroids
# - load the MedStat file with centroids
# - load the station locations
# - cut out the pop densities per shapefile
# - calculate the centroid per shapefile


library(sf)
library(raster)
library(ggplot2)
library(gridExtra)

### DATA #####
rm(list=ls())


# load shapefile of regions
mesh <- st_read("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MetStatRegions/centroids/shapefiles/MedStat_csr_adjusted.shp", quiet = TRUE)

# load shapefile of centroids
centroids <- st_read("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MetStatRegions/centroids/shapefiles/MedStat_centroids_popdensity.shp", quiet = TRUE)

# coordinates for the 8 stations
df_stations <- data.frame(station = c("Davos", "Chur", "Altdorf", "Montana", "Visp", "Magadino", "Lugano", "Poschiavo"),
                          x = c("2783519", "2759489", "2690181", "2601709", "2631151", "2715480", "2717874", "2801994"),
                          y = c("1187459", "1193182", "1193564", "1127489", "1128024", "1113161", "1095883", "1136249"),
                          MDSTID = c("GR06200", "GR06001", "UR03002", "VS09802", "VS09605", "TI08001", "TI08204", "GR06804"))


######



### DETERMINE LIST OF MEDSTAT REGIONS #####

# buffer_size_list = seq(4000, 15000, by = 1000)
buffer_size_list = 8000
buffer_list = list()

for(i in buffer_size_list){
  # buffer size 10k
  buffer_size <- i

  # list of Meshdat regions that are set to FALSE -> TRUE when they intersect with a buffer
  df_regions <- data.frame(MDST04= mesh$MDST04, Index = rep(0, length(mesh$MDSTID)))

  # calculate the buffer around the 8 stations
  for (i in 1:8) {

    # change format of station coordinates into sf as a point
    point_sf <- st_as_sf(df_stations[i,2:3], coords = c("x", "y"), crs = 2056)

    # Create the buffer
    buffer <- st_buffer(point_sf, dist = buffer_size)
    # save buffer
    buffer_list[[i]] = buffer

    # Find regions that intersect with the buffer
    intersects_with_buffer <- st_intersects(centroids$geometry, buffer, sparse = FALSE)

    # Add selected regions to regions list
    df_regions$Index <- df_regions$Index +  intersects_with_buffer

  }





  # Extract all regions that intersect in a list
  df_regions$ID <- rep(NA, nrow(mesh))

  # for every area that at least intersects once
  for(i in 1:nrow(df_regions)){
    if(df_regions$Index[i] > 0){

      # save the area code once again, all others NA
      df_regions$ID[i] <- df_regions$MDST04[i]

    }
  }


  # create a list of all medstat regions that fall within a buffer
  list_regions <- na.omit(df_regions$ID)
  # write.csv(list_regions, file = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MetStatRegions/centroids/MDSTID_MetRegions_", buffer_size, ".csv"))



  # plot the selected regions
  index_stations = match(df_stations$MDSTID, mesh$MDSTID)
  index_regions = match(list_regions, mesh$MDST04)


  map_plot <- ggplot() +
    geom_sf(data = mesh$geometry, fill = "grey", alpha = 0.5) +  # Adjust color and transparency
    geom_sf(data = mesh$geometry[index_regions], fill = "skyblue1", alpha = 0.7) +   # Adjust color and transparency
    geom_sf(data = centroids$geometry[index_regions], alpha = 0.7, cex = 0.5) +
    # geom_sf(data = mesh$geometry[index_stations], fill = "brown1", alpha = 0.5) +
    geom_sf(data = buffer_list[[1]], fill = "brown1", alpha = 0.4) +
    geom_sf(data = buffer_list[[2]], fill = "brown1", alpha = 0.4) +
    geom_sf(data = buffer_list[[3]], fill = "brown1", alpha = 0.4) +
    geom_sf(data = buffer_list[[4]], fill = "brown1", alpha = 0.4) +
    geom_sf(data = buffer_list[[5]], fill = "brown1", alpha = 0.4) +
    geom_sf(data = buffer_list[[6]], fill = "brown1", alpha = 0.4) +
    geom_sf(data = buffer_list[[7]], fill = "brown1", alpha = 0.4) +
    geom_sf(data = buffer_list[[8]], fill = "brown1", alpha = 0.4) +
    theme_minimal() +
    labs(title = paste0("selected regions with buffer: ", buffer_size))

  print(map_plot)
}
#####


### CONCLUSIONS #####
## for buffer size....:

# 2000, number of MedStat regions = 4
# 3000, n = 8, but none for Visp, Montana, Davos
# 4000, n = 13, but none for Visp and Davos, from here on I will sort regions to stations
# 5000, n = 18, but none for Davos
# 6000, n = 24 !all represented!
# 7000, n = 26
# 8000, 32
# 9000, 35
# 10k, 38
# 11k, 41
# 12k, 46
# 13k, 46
# 14k, 50
# 15k, 53 this does not really make sense anymore, I will stop sorting




#####





# ### FOR TESTING MEDSTAT #####
#
# map_plot2 <- ggplot() +
#   geom_sf(data = mesh$geometry, fill = "grey", alpha = 0.5) +
#   geom_sf(data = mesh$geometry[index_regions], fill = "skyblue1", alpha = 0.7) +
#   geom_sf(data = centroids$geometry[index_regions], alpha = 0.7) +
#   geom_sf(data = mesh$geometry[index_stations], fill = "brown1", alpha = 0.5) +
#
#   # geom_sf(data = mesh$geometry[mesh$MDST04 == "SZ05" ], fill = "green", alpha = 0.7) +
#   # geom_sf(data = mesh$geometry[mesh$MDST04 == "SG24" ], fill = "blue", alpha = 0.7) +
#   # geom_sf(data = mesh$geometry[mesh$MDST04 == "TI01" ], fill = "orange", alpha = 0.7) +
#
#   geom_sf(data = buffer_list[[1]], fill = "brown1", alpha = 0.4) +
#   geom_sf(data = buffer_list[[2]], fill = "brown1", alpha = 0.4) +
#   geom_sf(data = buffer_list[[3]], fill = "brown1", alpha = 0.4) +
#   geom_sf(data = buffer_list[[4]], fill = "brown1", alpha = 0.4) +
#   geom_sf(data = buffer_list[[5]], fill = "brown1", alpha = 0.4) +
#   geom_sf(data = buffer_list[[6]], fill = "brown1", alpha = 0.4) +
#   geom_sf(data = buffer_list[[7]], fill = "brown1", alpha = 0.4) +
#   geom_sf(data = buffer_list[[8]], fill = "brown1", alpha = 0.4) +
#   theme_minimal() +
#   labs(title = paste0("selected regions with buffer: ", buffer_size))
# print(map_plot2)

map_plot2 <- ggplot() +
  geom_sf(data = mesh$geometry, fill = "grey", alpha = 0.5) +
  geom_sf(data = mesh$geometry, fill = "skyblue1", alpha = 0.7) +
  geom_sf(data = centroids$geometry, alpha = 0.7, cex = 0.5) +
  geom_sf(data = mesh$geometry[index_stations], fill = "brown1", alpha = 0.5) +

  # geom_sf(data = mesh$geometry[mesh$MDST04 == "SZ05" ], fill = "green", alpha = 0.7) +
  # geom_sf(data = mesh$geometry[mesh$MDST04 == "SG24" ], fill = "blue", alpha = 0.7) +
  # geom_sf(data = mesh$geometry[mesh$MDST04 == "TI01" ], fill = "orange", alpha = 0.7) +

  # geom_sf(data = buffer_list[[1]], fill = "brown1", alpha = 0.4) +
  # geom_sf(data = buffer_list[[2]], fill = "brown1", alpha = 0.4) +
  # geom_sf(data = buffer_list[[3]], fill = "brown1", alpha = 0.4) +
  # geom_sf(data = buffer_list[[4]], fill = "brown1", alpha = 0.4) +
  # geom_sf(data = buffer_list[[5]], fill = "brown1", alpha = 0.4) +
  # geom_sf(data = buffer_list[[6]], fill = "brown1", alpha = 0.4) +
  # geom_sf(data = buffer_list[[7]], fill = "brown1", alpha = 0.4) +
  # geom_sf(data = buffer_list[[8]], fill = "brown1", alpha = 0.4) +
  theme_minimal() +
  labs(title = paste0("selected regions with buffer: ", buffer_size))
print(map_plot2)

grid.arrange(map_plot2, map_plot, ncol = 2)
#
# #####
