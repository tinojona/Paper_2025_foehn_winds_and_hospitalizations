################################################################################
### In this file I:
### - extracted the foehn data from meteoschweiz (deleted the first three empty rows)
### - converted UTC to CET
### - aggregated the 10 minute foehn index to daily counts
###                                       - alternative, take threshold of eg 6 hours and then assign foehn or no foehn
### - saved this data in a folder below the data folder of this project

library(lubridate) # for UTC to CET
rm(list = ls())

# Specify the folder path
folder_path <- "C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data-raw/foehn/data"

# List all files in the folder
file_names <- list.files(path = folder_path)

# create empty NA record
foehn_NA = data.frame(stations = rep(NA,8),
                      raw = rep(NA,8),
                      aggregated = rep(NA,8))

for(i in 1:8){

  # load complete path to file
  file = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data-raw/foehn/data/", file_names[i])

  # load data
  data = read.table(file, header = TRUE)

  # Convert to date-time object using strptime
  data$time_conv <- as.POSIXct(strptime(data$time, format = "%Y%m%d%H%M"), tz = "UTC")

  # convert time_conv from UTC to CET
  data$time_conv = as.Date(with_tz(data$time_conv, tzone = Sys.timezone()), format = "%Y-%m-%d %H:%M:%s")


  # convert foehn data to numeric
  data$f_id <- as.numeric(data$wcc006s0)

  # aggregate by day, sum
  data_agg_daily_sum <- aggregate(f_id ~ time_conv, data = data, FUN = function(x) sum(x, na.rm = TRUE))

  # extract station abbreviation
  station_abbr = substring(file_names[i], 14,16)

  # path for file
  path_for_file = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/foehn_processed/",station_abbr,"_daily_aggregated.csv")

  # save aggregated file
  write.csv(data_agg_daily_sum, file = path_for_file)

  # extract and save perc of NA in f_id before and after aggegation
  foehn_NA$stations[i] = station_abbr
  foehn_NA$raw[i] = sum(is.na(data)) / nrow(data)
  foehn_NA$aggregated[i] = sum(is.na(data_agg_daily_sum)) / nrow(data_agg_daily_sum)



}

# write NA file
path_for_NA = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/NA/NA_record_foehn.csv")
write.csv(foehn_NA, file = path_for_NA)
