################################################################################
### In this file I:
### - check how many dates are missing for every station and save that in the NA folder

### REASON WHY I HAVE SO MANY NAS IN MY FINAL DATA:
### after merging hosp and foehn data, I get a lot of NAs in foehn because,
### foehn records only begin in 2006 and after for Magadino and Lugano but hosp in 1998.

### I have no missing dates in my raw data but some NAs
### These NAs lead to missing dates that are then filled with NAs when merging with hosp data

rm(list=ls())


# load the raw data
files_list = list.files(path = "C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data-raw/foehn/data")

# create empty dates percentage record
foehn_dates = data.frame(stations = rep(NA,8),
                     missing = rep(NA,8))


for(i in 1:8){

  # load complete path to file
  file = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data-raw/foehn/data/", files_list[i])

  # load data
  data = read.table(file, header = TRUE)

  # Convert to date-time object using strptime
  data$time_conv <- as.POSIXct(strptime(data$time, format = "%Y%m%d%H%M"), tz = "UTC")

  # convert time_conv from UTC to CET
  data$time_conv = as.Date(with_tz(data$time_conv, tzone = Sys.timezone()), format = "%Y-%m-%d %H:%M:%s")

  # get number of individual days in og dataset
  n_dates_og = length(unique(data$time_conv))

  # make day series of dates
  start_date = min(data$time_conv)
  end_date = max(data$time_conv)
  series = seq.Date(from = as.Date(start_date), to = as.Date(end_date), by = "day")

  # calculate ratio
  ratio_og = n_dates_og / length(series)


  # get station name
  station_abbr = substring(files_list[i], 14,16)

  # save number of presented dates in data frame
  foehn_dates$stations[i] = station_abbr
  foehn_dates$missing[i] = ratio_og


}

# write NA file
path_for_NA = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/NA/NA_dates_raw.csv")
write.csv(foehn_dates, file = path_for_NA)


### one more time for the aggregated data that has been processed
rm(list=ls())
# load the raw data
files_list = list.files(path = "C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/foehn_processed")

# create empty dates percentage record
foehn_dates = data.frame(stations = rep(NA,8),
                         missing = rep(NA,8))


for(i in 1:8){

  # load complete path to file
  file = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/foehn_processed/", files_list[i])

  # load data
  data = read.csv(file, header = TRUE)

  #  Convert to date object
  data$time_conv <- as.Date(data$time_conv)

  # get number of individual days in og dataset
  n_dates_og = length(unique(data$time_conv))

  # make day series of dates
  start_date = min(data$time_conv)
  end_date = max(data$time_conv)
  series = seq.Date(from = as.Date(start_date), to = as.Date(end_date), by = "day")

  # calculate ratio
  ratio_og = n_dates_og / length(series)


  # get station name
  station_abbr = substring(files_list[i], 1,3)

  # save number of presented dates in data frame
  foehn_dates$stations[i] = station_abbr
  foehn_dates$missing[i] = ratio_og
  # # aggregate the data
  # # convert foehn data to numeric
  # data$f_id <- as.numeric(data$wcc006s0)
  #
  # # aggregate by day, sum
  # data_agg_daily_sum <- aggregate(f_id ~ time_conv, data = data, FUN = function(x) sum(x, na.rm = TRUE))
  #
  # # extract station abbreviation
  # station_abbr = substring(file_names[i], 14,16)

}

# write NA file
path_for_NA = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/NA/NA_dates_aggr.csv")
write.csv(foehn_dates, file = path_for_NA)




## load NAs percentage in raw data
## load missing dates percentage in aggregated data
rm(list=ls())

data_raw = read.csv("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/NA/NA_record_foehn.csv")

data_aggr = read.csv("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/NA/NA_dates_aggr.csv")


summer = data_raw$raw + data_aggr$missing




## lets look into the final product
rm(list=ls())

data = read.csv("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/MedStat_aggregated/CLM_buffersize_5000.csv")

1 - (sum(complete.cases(data)) / nrow(data))
# -> 16% missing data





