


# Specify the folder path
folder_path <- "C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data-raw/foehn/data"

# List all files in the folder
file_names <- list.files(path = folder_path)

for(i in file_names){

  # load complete path to file
  file = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data-raw/foehn/data/", i)

  # load data
  data = read.table(file, header = TRUE)

  # Convert to date-time object using strptime
  data$time_conv <- as.Date(strptime(data$time, format = "%Y%m%d%H%M"), format = "%Y%m%d%H%M")

  # convert foehn data to numeric
  data$f_id <- as.numeric(data$wcc006s0)

  # aggregate by day, sum
  data_agg_daily_sum <- aggregate(f_id ~ time_conv, data = data, sum)

  # extract station abbreviation
  station_abbr = substring(i, 14,16)

  # path for file
  path_for_file = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/foehn_processed/",station_abbr,"_daily_aggregated.csv")

  # save aggregated file
  write.csv(data_agg_daily_sum, file = path_for_file)


}


