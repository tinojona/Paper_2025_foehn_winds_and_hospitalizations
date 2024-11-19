################################################################################
# Determine the MedStat regions for different buffers that overlap with the centroids
# - load the MedStat file with centroids
# - load the station locations
# - cut out the pop densities per shapefile
# - calculate the centroid per shapefile


library(sf); library(raster); library(ggplot2); library(gridExtra); library(viridis); library(terra)

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



# PLOT FOR THE PAPER


# colors
colors <- viridis(8, option = "viridis")

# station locations
station_loc <- st_as_sf(df_stations[,2:3], coords = c("x", "y"), crs = 2056)

# raster data from Swisstopo https://www.swisstopo.admin.ch/en/height-model-dhm25-200m#Additional-information
raster_data <- rast("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data-raw/elevation/DHM200.asc")
crs(raster_data) <- "EPSG:21781" # assign projection based on LV03 LN02


# reproject shapefiles based on raster
crs_raster <- crs(raster_data)

# print(crs_raster)

mesh$geometry <- st_transform(mesh$geometry, crs = crs_raster)
station_loc$geometry <- st_transform(station_loc$geometry, crs = crs_raster)

# print(crs(mesh$geometry))
# print(crs(station_loc$geometry))

# crop the raster to the extend of the shapefile
mesh_spatvector <- terra::vect(mesh$geometry)
raster_cropped <- terra::mask(raster_data, mesh_spatvector)

# convert to data frame for plotting
raster_df <- as.data.frame(raster_cropped, xy = TRUE, na.rm = TRUE); colnames(raster_df) <- c("x", "y", "elevation")



png("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/plots/paper/map_plot.png", width = 1800, height = 1400, res = 300)


ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = elevation)) +
  # scale_fill_viridis_c() +
  scale_fill_gradient(
    low = "white", high = "black", # Set the color range for the elevation
    name = "Elevation [m]",           # Title for the color legend
    limits = c(0, 4500), # Set limits for the legend
    na.value = "transparent"      # Set how NA values are represented
  ) +
  coord_fixed() +
  theme_minimal() +
  xlab("") +
  ylab("") +

  geom_sf(data = mesh$geometry, fill = NA, color = "white") +
  geom_sf(data = mesh$geometry[index_regions[1:3]], fill = colors[1], alpha = 0.6, color = "black") +
  geom_sf(data = mesh$geometry[index_regions[4:6]], fill = colors[6], alpha = 0.6, color = "black") +
  geom_sf(data = mesh$geometry[index_regions[7]], fill = colors[4], alpha = 0.6, color = "black") +
  geom_sf(data = mesh$geometry[index_regions[8]], fill = colors[2], alpha = 0.6, color = "black") +
  geom_sf(data = mesh$geometry[index_regions[9:13]], fill = colors[8], alpha = 0.6, color = "black") +
  geom_sf(data = mesh$geometry[index_regions[14:25]], fill = colors[5], alpha = 0.6, color = "black") +
  geom_sf(data = mesh$geometry[index_regions[26:27]], fill = colors[3], alpha = 0.6, color = "black") +
  geom_sf(data = mesh$geometry[index_regions[28:32]], fill = colors[7], alpha = 0.6, color = "black") +
  geom_sf(data = station_loc, alpha = 1, fill = "black", stroke = .2) +

  # # Set the color scale for the elevation data (fill based on elevation)

  #
  theme(legend.position = "top",    # Position the legend on the right side
        legend.title = element_text(size = 10),  # Font size for the legend title
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.margin = unit(c(0,0,0,0), "cm")) +

# Customize legend appearance (optional)
  guides(fill = guide_colorbar(barwidth = 10, barheight = .5, title.position = "top", title.hjust = 0.5)) #+



dev.off()





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





# # ### FOR TESTING MEDSTAT #####
# #
# # map_plot2 <- ggplot() +
# #   geom_sf(data = mesh$geometry, fill = "grey", alpha = 0.5) +
# #   geom_sf(data = mesh$geometry[index_regions], fill = "skyblue1", alpha = 0.7) +
# #   geom_sf(data = centroids$geometry[index_regions], alpha = 0.7) +
# #   geom_sf(data = mesh$geometry[index_stations], fill = "brown1", alpha = 0.5) +
# #
# #   # geom_sf(data = mesh$geometry[mesh$MDST04 == "SZ05" ], fill = "green", alpha = 0.7) +
# #   # geom_sf(data = mesh$geometry[mesh$MDST04 == "SG24" ], fill = "blue", alpha = 0.7) +
# #   # geom_sf(data = mesh$geometry[mesh$MDST04 == "TI01" ], fill = "orange", alpha = 0.7) +
# #
# #   geom_sf(data = buffer_list[[1]], fill = "brown1", alpha = 0.4) +
# #   geom_sf(data = buffer_list[[2]], fill = "brown1", alpha = 0.4) +
# #   geom_sf(data = buffer_list[[3]], fill = "brown1", alpha = 0.4) +
# #   geom_sf(data = buffer_list[[4]], fill = "brown1", alpha = 0.4) +
# #   geom_sf(data = buffer_list[[5]], fill = "brown1", alpha = 0.4) +
# #   geom_sf(data = buffer_list[[6]], fill = "brown1", alpha = 0.4) +
# #   geom_sf(data = buffer_list[[7]], fill = "brown1", alpha = 0.4) +
# #   geom_sf(data = buffer_list[[8]], fill = "brown1", alpha = 0.4) +
# #   theme_minimal() +
# #   labs(title = paste0("selected regions with buffer: ", buffer_size))
# # print(map_plot2)
#
# map_plot2 <- ggplot() +
#   geom_sf(data = mesh$geometry, fill = "grey", alpha = 0.5) +
#   geom_sf(data = mesh$geometry[index_regions], fill = "skyblue1", alpha = 0.7) +
#   geom_sf(data = centroids$geometry, alpha = 0.7, cex = 0.5) +
#   geom_sf(data = mesh$geometry[index_stations], fill = "brown1", alpha = 0.5) +
#
#   # geom_sf(data = mesh$geometry[mesh$MDST04 == "SZ05" ], fill = "green", alpha = 0.7) +
#   # geom_sf(data = mesh$geometry[mesh$MDST04 == "SG24" ], fill = "blue", alpha = 0.7) +
#   # geom_sf(data = mesh$geometry[mesh$MDST04 == "TI01" ], fill = "orange", alpha = 0.7) +
#
#   # geom_sf(data = buffer_list[[1]], fill = "brown1", alpha = 0.4) +
#   # geom_sf(data = buffer_list[[2]], fill = "brown1", alpha = 0.4) +
#   # geom_sf(data = buffer_list[[3]], fill = "brown1", alpha = 0.4) +
#   # geom_sf(data = buffer_list[[4]], fill = "brown1", alpha = 0.4) +
#   # geom_sf(data = buffer_list[[5]], fill = "brown1", alpha = 0.4) +
#   # geom_sf(data = buffer_list[[6]], fill = "brown1", alpha = 0.4) +
#   # geom_sf(data = buffer_list[[7]], fill = "brown1", alpha = 0.4) +
#   # geom_sf(data = buffer_list[[8]], fill = "brown1", alpha = 0.4) +
#   theme_minimal() +
#   # labs(title = paste0("selected regions with buffer: ", buffer_size))
# print(map_plot2)
#
# grid.arrange(map_plot2, map_plot, ncol = 2)
#
# #####


