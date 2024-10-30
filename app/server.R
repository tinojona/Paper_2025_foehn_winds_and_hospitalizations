# server.R

# Required libraries
library(dplyr)
library(sf)
library(raster)
library(ggplot2)
library(gridExtra)
library(zoo)
library(tidyr)
library(plotly)


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
# Define server logic
server <- function(input, output) {

  # Reactive expression to load and preprocess data based on buffer radius input
  data_means <- reactive({
    buffer <- as.numeric(input$buffer_radius)  # Get buffer radius from input

    # Load the data dynamically based on buffer radius
    data <- read.csv(paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MedStat_aggregated/centroid_aggregated/hosp_buffer_", buffer, ".csv"))

    # Preprocess data
    data$date <- as.Date(data$date)
    data$station <- as.factor(data$station)

    # Index to include only stratum that have hosp counts
    data$stratum_dow <- as.factor(data$stratum_dow)
    data$stratum <- as.factor(data$stratum)
    ind_dow <- tapply(data$all, data$stratum_dow, sum)
    ind <- tapply(data$all, data$stratum, sum)

    data <- data %>%
      mutate(y64 = a014y + a1564y) %>%
      mutate(o64 = a6574y + a7584y + a85plusy)

    # daily means
    data_daily_mean = data %>%
      mutate(daymonth = format(date, "%m-%d")) %>%
      group_by(daymonth) %>%
      summarise(
        across(c(all, mal, fem, y64, o64, cvd, resp, temp), mean),
        mean_f_id = mean(f_id),
        p90_f_id = quantile(f_id, 0.9) # Calculate 90th percentile
      ) %>%
      mutate(daymonth = as.Date(paste0("2000-", daymonth)))

    return(data_daily_mean)  # Return the processed data
  })

  # Plot the buffer radius plot
  output$selected_radius_plot <- renderPlot({
    buffer_size <- as.numeric(input$buffer_radius)  # Get the buffer size directly from input
    buffer_list <- list()

    # Initialize region dataframe to keep track of regions intersecting the buffer
    df_regions <- data.frame(MDST04 = mesh$MDST04, Index = rep(0, length(mesh$MDSTID)))

    # Calculate buffers around the 8 stations
    for (i in 1:8) {
      point_sf <- st_as_sf(df_stations[i, 2:3], coords = c("x", "y"), crs = 2056)

      # Create and save the buffer
      buffer <- st_buffer(point_sf, dist = buffer_size)
      buffer_list[[i]] <- buffer

      # Find regions intersecting with the buffer
      intersects_with_buffer <- st_intersects(centroids$geometry, buffer, sparse = FALSE)
      df_regions$Index <- df_regions$Index + intersects_with_buffer
    }

    # Extract regions that intersect
    df_regions$ID <- rep(NA, nrow(mesh))
    for (j in 1:nrow(df_regions)) {
      if (df_regions$Index[j] > 0) {
        df_regions$ID[j] <- df_regions$MDST04[j]
      }
    }

    # Create a list of regions within the buffer
    list_regions <- na.omit(df_regions$ID)

    # Match station and region indices for plotting
    index_stations <- match(df_stations$MDSTID, mesh$MDSTID)
    index_regions <- match(list_regions, mesh$MDST04)

    # Plotting the map
    map_plot <- ggplot() +
      geom_sf(data = mesh$geometry, fill = "grey", alpha = 0.5) +
      geom_sf(data = mesh$geometry[index_regions], fill = "skyblue1", alpha = 0.7) +
      geom_sf(data = centroids$geometry[index_regions], alpha = 0.7, cex = 0.5) +
      theme_minimal() +
      labs(title = paste0("Selected regions with buffer size: ", buffer_size))

    # Add each buffer as a layer to the map
    for (buf in buffer_list) {
      map_plot <- map_plot + geom_sf(data = buf, fill = "brown1", alpha = 0.4)
    }

    print(map_plot)
  })

  output$histograms <- renderPlot({
    # req(input$generate_stats)  # Ensure the checkbox is checked
    # data_daily_mean_use <- data_means()  # Call the reactive expression to get the latest data
    # req(data_daily_mean_use)  # Ensure data_daily_mean is available

    # Generate histograms for multiple variables
    p1 <- ggplot(data[which(data$f_id != 0),], aes(x = f_id[which(data$f_id != 0)])) +
      geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
      theme_minimal() + labs(title = "Histogram of f_id", x = "f_id", y = "Count")

    p2 <- ggplot(data, aes(x = all)) +
      geom_histogram(binwidth = 1, fill = "orange", color = "black") +
      theme_minimal() + labs(title = "Histogram of all", x = "all", y = "Count")

    p3 <- ggplot(data, aes(x = temp)) +
      geom_histogram(binwidth = 1, fill = "brown1", color = "black") +
      theme_minimal() + labs(title = "Histogram of temp", x = "temp", y = "Count")

    gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
  })


}


