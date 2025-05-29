################################################################################
### HOSPITALIZATION PROCESSING
################################################################################
### In this file I:
### - aggregated the hospitalization data of the regions (depending on the buffer)
###   to their corresponding meteorological station
### - as some days were missing, I added the dates and inserted 0 hospitalizations!
###   and added a column with the day of the week
### - to this data, I saved the corresponding foehn and temp data, maybe inducing some NAs
###   for dates that no temp / foehn data was apparent
### - I saved the data per buffer size combined of all stations



### DATA
#----

rm(list=ls())

# libraries
library(tidyverse); library(lubridate)

# read hospitalization data
data = read.csv("../data-raw/MedStat/hosp_daily.csv", header =T) |>
  mutate(station = NA)

# get foehn data for every station
files_foehn = list.files("data_nonsensitive/foehn_wind_station_aggregated/")

# CAREFUL! Sensitivity analysis with only full foehn in aggregation
# files_foehn = files_foehn[grepl("sensitivity", files_foehn)]

# get temperature data for every station
files_temp = list.files("data_nonsensitive/temperature_station_processed")

# get regions sorted to station files
files_regions = list.files("data_nonsensitive/Medstat_regions_per_buffer_size")

#----


### FOR LOOP over all buffer sizes and consequently Medstat combinations
#----

# define the buffer size you are interested in
buffer_sizes = c(5000,6000,7000,seq(9000,15000,1000))

# start loop for 4 buffers
for(i in buffer_sizes){

  # save buffer for file writing in the end
  buffer = i

  # get the file for the selected buffer size
  medstat_file = files_regions[grepl( paste0("_",as.character(i)), files_regions )]

  # load the regions file
  regions = read.csv(paste0("data_nonsensitive/Medstat_regions_per_buffer_size/", medstat_file))

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
  for (j in 1:ncol(regions)) {

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
    aggregated_by_station <- data |>

      # Filter out rows where station is NA
      filter(station == station_current) |>

      # group by date
      group_by(DT_EINTRITTSDAT) |>

      # sum across hospitalization subpopulations
      summarise(across(all:uri, sum),
                station = station_current) |>

      # rename date column so that it matches with foehn and temp later
      rename("date" = DT_EINTRITTSDAT)


    ## insert dates that got excluded by the aggregation process
    start_date = min(aggregated_by_station$date, na.rm = TRUE)
    end_date = max(aggregated_by_station$date, na.rm = TRUE)
    end_date_alter = "2019-12-31"

    if(as.Date(end_date) >= as.Date(end_date_alter)){
      end_date <- end_date_alter
    }



    # create empty data frame with continuous dates
    df_with_dates = data.frame(date = as.character(seq.Date(from = as.Date(start_date), to = as.Date(end_date), by = "day")))


    # merge empty data frame with time series with hospitalization data
    df_full <- df_with_dates %>%
      left_join(aggregated_by_station, by = "date")


    # fill NA in station with station name
    df_full$station = station_current

    # fill NA in values with 0
    df_full[is.na(df_full)] <- 0

    # reassign proper data frame name
    aggregated_by_station <- df_full




    # cbind the foehn data to the appropriate station
    # find the correct file
    station_current_up = substr(toupper(station_current), 1,3)

    matching_strings <- grep(station_current_up, files_foehn, value = TRUE)

    foehn_data = read.csv(paste0("data_nonsensitive/foehn_wind_station_aggregated/", matching_strings)) |>
      rename("date" = time_conv)

    aggregated_by_station <- aggregated_by_station %>%
      left_join(foehn_data[,2:3], by = "date")




    # cbind the temp data
    # find the correct file
    matching_strings <- grep(station_current_up, files_temp, value = TRUE)
    temp_data = read.csv(paste0("data_nonsensitive/temperature_station_processed/", matching_strings)) |>
      mutate(date = time_conv)

    aggregated_by_station <- aggregated_by_station |>
      left_join(temp_data[,3:4], by = "date")




    # create day of the week column
    aggregated_by_station <- aggregated_by_station |>
      mutate(date = as.Date(date),

             # create base variables for the stratum
             dow = as.factor(weekdays(date)),
             year = as.factor(format(date, "%Y")),
             month = as.factor(format(date, "%m")),

             # create stratum
             stratum_dow = paste(station, year, month, dow, sep="-"),
             stratum = paste(station, year, month, sep="-")
             )


    # append the data to the buffer aggregated dataset
    aggregated_by_buffer = rbind(aggregated_by_buffer, aggregated_by_station)

  }

  # only allow complete cases
  aggregated_by_buffer <- aggregated_by_buffer[complete.cases(aggregated_by_buffer), ]

  # save the combined buffer set
  write.csv(aggregated_by_buffer, file = paste0("../data/Medstat_hospitalizations_aggregated/hosp_buffer_", buffer, ".csv"))

}


