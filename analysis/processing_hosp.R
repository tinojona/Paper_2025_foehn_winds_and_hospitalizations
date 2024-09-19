


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

# get foehn data for every station
files_foehn = list.files("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/foehn_processed")

# get temperature data for every station
files_temp = list.files("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/temp_processed")

# start loop for 4 buffers
for(i in 1:4){

  # get buffer size and sorted regions
  buffer = buffer_sizes[i]
  regions = read.csv(paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/MetStatRegions/", files_regions[i]))

  # get colnames of the stations for later
  station_names = colnames(regions)

  # create empty aggregated_by_buffer data.frame for later here
  aggregated_by_buffer = data.frame(DT_EINTRTTSDAT = character(),
                                    all = numeric(),
                                    a014 = numeric(),
                                    a1564y = numeric(),
                                    a6574y = numeric(),
                                    a7584y = numeric(),
                                    a85plusy = numeric(),
                                    mal = numeric(),
                                    fem = numeric(),
                                    inf = numeric(),
                                    ment = numeric(),
                                    cvd = numeric(),
                                    resp = numeric(),
                                    uri = numeric(),
                                    station = character())

  # start loop for different station names
  for (j in 1:8) {

    # to clear the index for later
    data$station = NA

    # current station name
    station_current = station_names[j]
    # print(station_current) # for code checking

    # get all region names of that station, exclude empty entries
    region_names = regions[regions[,j] != "",j ]
    # print(region_names) # for code checking

    # assign the station name to data$station when one of the region names is present
    # first create index
    index <- sapply(data$ID_WOHNREGION, function(row_text) {
      any(sapply( region_names,function(word) grepl(word, row_text, ignore.case = TRUE)))
    })

    data$station[index] = station_current


    # summarize when the station name is present in data$station by date
    aggregated_by_station <- data %>%
      filter(station == station_current) %>%  # Filter out rows where station is NA
      group_by(DT_EINTRITTSDAT) %>%           # Group by the date column
      summarise(across(all:uri, sum),         # across al variables we sum
                station = station_current)    # but we keep the station name

    # rename date column so that it matches with foehn and temp later
    aggregated_by_station <- aggregated_by_station %>%
      rename("date" = DT_EINTRITTSDAT)

    # insert dates that got excluded by the aggregation process




    # cbind the foehn data to the appropriate station
    foehn_data = read.csv(paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/foehn_processed/", files_foehn[j]))
    foehn_data <- foehn_data %>%
      rename("date" = time_conv)

    aggregated_by_station <- aggregated_by_station %>%
      left_join(foehn_data[,2:3], by = "date")


    # cbind the temp data
    temp_data = read.csv(paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/temp_processed/", files_temp[j]))
    temp_data$date <- as.character(as.Date(strptime(as.character(temp_data$time), format = "%Y%m%d")))

    aggregated_by_station <- aggregated_by_station %>%
      left_join(temp_data[,3:4], by = "date")


    # save the data per station and buffer
    file_name = paste0(station_current, "_buffer_", buffer, ".csv")
    write.csv(aggregated_by_station, file = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/MedStat_aggregated/by_station/", file_name))


    # append the data to the buffer aggregated dataset
    aggregated_by_buffer = rbind(aggregated_by_buffer, aggregated_by_station)

  }

  # save the combined buffer set
  write.csv(aggregated_by_buffer, file = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/MedStat_aggregated/hosp_buffer_", buffer, ".csv"))

}


