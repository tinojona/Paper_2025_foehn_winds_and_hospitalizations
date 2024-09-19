


# aggregate by station
# while sorting areas to stations
# for different buffers

# packages
library(dplyr)
rm(list = ls())


# read hospitalization data
data = read.csv("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data-raw/MedStat/med_stat_data/hosp_daily.csv", header =T)
data$station = NA

# get regions sorted to station files
files_regions = list.files("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/MetStatRegions")[c(2,4,7,10)]

# recreate all the buffer sies that are in files_regions, in THE SAME ORDER
buffer_sizes = c(1000,10000,2500,5000)


# start loop for 4 buffers
for(i in 1:4){

  # get buffer size and sorted regions
  buffer = buffer_sizes[i]
  regions = read.csv(paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/MetStatRegions/", files_regions[i]))

  # get colnames of the stations for later
  station_names = colnames(regions)

  # start loop for different station names
  for (j in 1:8) {

    # current station name
    station_current = station_names[j]

    # get all region names of that station
    region_names = regions[regions[,j] != "",j ]

    # assign the station name to data$station when one of the region names is present
    # first create index
    index <- sapply(data$ID_WOHNREGION, function(row_text) {
      any(sapply( region_names,function(word) grepl(word, row_text)))
    })

    data$station[index] = station_current


    # summarize when the station name is present in data$station by date
    aggregated_by_station <- data %>%
      filter(station == station_current) %>%  # Filter out rows where station is NA
      group_by(DT_EINTRITTSDAT) %>%           # Group by the date column
      summarise(across(all:uri, sum),         # across al variables we sum
                station = station_current)    # but we keep the station name

    # add temperature and foehn data to the aggregated_by_station


    # save the data
    file_name = paste0(station_current, "_buffer_", buffer, ".csv")
    write.csv(aggregated_by_station, file = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/MedStat_aggregated/by_station/", file_name))

  }

  # save the data of one buffer in one data.frame
  # aggregated_by_buffer = rbind(aggregated_by_buffer ,aggregated_by_station)

}


