################################################################################
### FOEHN WIND PROCESSING
################################################################################
### - load the foehn data from meteoschweiz (deleted the first three empty rows)
### - converted UTC to CET
### - aggregated the 10 minute foehn index to daily counts
### - define a binary variable using threshold of 6 hours of full foehn wind
### - saved this data in a folder below the data folder of this project

# CAREFUL!
# in the output files:
# change the accronyms: MVE -> MON & ROB -> POS

### Load data
#----

rm(list = ls())

# libraries
library(lubridate); library(tidyverse)

# List all files in the foehn wind folder
file_names <- list.files(path = "data/environmental_RAW/foehn_wind_station_raw/data")

#----




### FOR LOOP OVER ALL STATIONS
#----
for(i in 1:length(file_names)){

  # load complete path to file
  file = paste0("data/environmental_RAW/foehn_wind_station_raw/data/", file_names[i])

  # load data
  data = read.table(file, header = TRUE) |>
    mutate(
      # convert time column from seconds since 1970, change timezone and drop subdaily information
      time_conv = as.POSIXct(strptime(time, format = "%Y%m%d%H%M"), tz = "UTC"),
      time_conv = as.Date(with_tz(time_conv, tzone = Sys.timezone()), format = "%Y-%m-%d %H:%M:%s"),

      # convert foehn data to numeric
      f_id = as.numeric(wcc006s0),

      # CAREFUL ! SENSITIVITY ANALYSIS WITH DIFFERENT FOEHN WINDS AGGREGATION
      f_id_sens = ifelse(f_id == 2, 2, 0)
           )


  # aggregate by day, sum
  data_agg_daily_sum <- aggregate(f_id ~ time_conv, data = data, FUN = function(x) sum(x, na.rm = TRUE))
  # data_agg_daily_sum <- aggregate(f_id_sens ~ time_conv, data = data, FUN = function(x) sum(x, na.rm = TRUE))


  data_agg_daily_sum <- data_agg_daily_sum |>
    mutate(f_id_bin = ifelse(f_id > 72, 0, 1),
           f_id_bin_rev = ifelse(f_id_bin == 0, 1, 0))

  # data_agg_daily_sum <- data_agg_daily_sum |>
  #   mutate(f_id_bin = ifelse(f_id_sens > 72, 0, 1),
  #          f_id_bin_rev = ifelse(f_id_bin == 0, 1, 0))

  # extract station abbreviation
  station_abbr = substring(file_names[i], 14,16)

  # path for file
  path_for_file = paste0("data/environmental_processed/foehn_wind_station_aggregated/",station_abbr,"_daily_aggregated.csv")
  # path_for_file = paste0("data_nonsensitive/foehn_wind_station_aggregated/",station_abbr,"_daily_aggregated_sensitivity_onlyfullfoehnaggregation.csv")

  # save aggregated file
  write.csv(data_agg_daily_sum, file = path_for_file)

}

#----
