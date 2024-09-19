################################################################################
### In this file I:
### - aggregated the hospitalization data of the regions (depending on the buffer)
###   to their corresponding meteorological station
### - as some days were missing, I added the dates and inserted 0 hospitalizations!
###   and added a column with the day of the week
### - to this data, I saved the corresponding foehn and temp data, maybe inducing some NAs
###   for dates that no temp / foehn data was apparent
### - I saved both the data per station, and per buffer size combined of all stations




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
                                    station = character(),
                                    dow = character())

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

    ## assign the station name to data$station when one of the region names is present
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

    ## insert dates that got excluded by the aggregation process
    start_date = min(aggregated_by_station$date)
    end_date = max(aggregated_by_station$date)
    end_date_alter = "2019-12-31"

    if(end_date > end_date_alter){
      end_date <- end_date_alter
    }

    # create empty df with continous dates
    df_with_dates = data.frame(date = as.character(seq.Date(from = as.Date(start_date), to = as.Date(end_date), by = "day")))

    # merge df with hosp data
    df_full <- df_with_dates %>%
      left_join(aggregated_by_station, by = "date")

    # fill NA in station with station name
    df_full$station = station_current

    # fill NA in values with 0
    df_full[is.na(df_full)] <- 0

    # reasign proper df name
    aggregated_by_station <- df_full

    print(nrow(aggregated_by_station)) # TODO delete later



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

    print(nrow(aggregated_by_station)) # TODO delete later


    # create day of the week column
    aggregated_by_station$date <- as.Date(aggregated_by_station$date)
    aggregated_by_station$dow <- weekdays(aggregated_by_station$date)


    # save the data per station and buffer
    file_name = paste0(station_current, "_buffer_", buffer, ".csv")
    write.csv(aggregated_by_station, file = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/MedStat_aggregated/by_station/", file_name))


    # append the data to the buffer aggregated dataset
    aggregated_by_buffer = rbind(aggregated_by_buffer, aggregated_by_station)

  }

  # save the combined buffer set
  write.csv(aggregated_by_buffer, file = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/MedStat_aggregated/hosp_buffer_", buffer, ".csv"))

}


