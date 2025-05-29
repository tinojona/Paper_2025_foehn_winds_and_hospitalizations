################################################################################
### TEMPERATURE PROCESSING
################################################################################
### In this file I:
### - extracted the temperature data from meteoschweiz (deleted the first three empty rows)
### - saved the daily mean temperature an the date in an extra file

# CAREFUL!
# in the output files:
# change the accronyms: MVE -> MON & ROB -> POS

### DATA
#----

rm(list=ls())

# libraries
library(lubridate); library(tidyverse)

# List all files in the folder
file_names <- list.files(path = "data_nonsensitive/temperature_station_raw/data/")

#----



### LOOP OVER ALL FILES
#----

for(i in file_names){

  # load complete path to file
  file = paste0("data_nonsensitive/temperature_station_raw/data/", i)

  # load data
  data = read.table(file, header = TRUE) |>
    mutate(
      time = as.character(time),

      # Convert to date-time object using strptime
      time_conv = as.Date(time, format = "%Y%m%d"),

      # rename temperature column
      temp = tre200d0
    ) |>
    select(temp, time_conv)


  # extract station abbreviation
  station_abbr = substring(i, 14,16)

  # path for file
  path_for_file = paste0("data_nonsensitive/temperature_station_processed/temp_",station_abbr,".csv")

  # save aggregated file
  write.csv(data, file = path_for_file)


}


