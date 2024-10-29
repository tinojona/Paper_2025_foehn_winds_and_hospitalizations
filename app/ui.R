
#
# ### LIBRARIES ----
# library(sf)
# library(raster)
# library(ggplot2)
# library(gridExtra)
#
# #####
#
#
# ### DATA ----
# # coordinates for the 8 stations
# df_stations <- data.frame(station = c("Davos", "Chur", "Altdorf", "Montana", "Visp", "Magadino", "Lugano", "Poschiavo"),
#                           x = c("2783519", "2759489", "2690181", "2601709", "2631151", "2715480", "2717874", "2801994"),
#                           y = c("1187459", "1193182", "1193564", "1127489", "1128024", "1113161", "1095883", "1136249"),
#                           MDSTID = c("GR06200", "GR06001", "UR03002", "VS09802", "VS09605", "TI08001", "TI08204", "GR06804"))
#
# # load shapefile of regions
# mesh <- st_read("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MetStatRegions/centroids/shapefiles/MedStat_csr_adjusted.shp", quiet = TRUE)
#
# # load shapefile of centroids
# centroids <- st_read("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MetStatRegions/centroids/shapefiles/MedStat_centroids_popdensity.shp", quiet = TRUE)
#
# #####
#




ui <- page_sidebar(

  # App title
  title = "Foehn as a modifier of the temperature hospitalization association",

  # Buffer radius selection (in meters)
  sidebar = sidebar(
    helpText("The buffer radius determines the MedStat regions that will be included in the further analysis"),

    selectInput("buffer_radius",
                label = "Choose a buffer radius [m]",
                choices = list(2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000),
                selected = 8000),

    helpText("From experience: 8000m buffer radius shows the largest effect of foehn")
  ),

  # Enclosed area for plot
  card(
    card_header("Figure of the selected MedStat regions")
  ),

  # save output of buffer radius for server file
  textOutput("selected_radius")
)
