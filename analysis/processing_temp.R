


# Specify the folder path
folder_path <- "C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data-raw/temp/data"

# List all files in the folder
file_names <- list.files(path = folder_path)

for(i in file_names){

  # load complete path to file
  file = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data-raw/temp/data/", i)

  # load data
  data = read.table(file, header = TRUE)

  data$time <- as.character(data$time)

  # Convert to date-time object using strptime
  data$time_conv <- as.Date(data$time, format = "%Y%m%d")

  # rename temperature column
  names(data)[names(data) == "tre200d0"] <- "temp"


  # extract station abbreviation
  station_abbr = substring(i, 14,16)

  # path for file
  path_for_file = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/temp_processed/temp_",station_abbr,".csv")

  # save aggregated file
  write.csv(data[,2:3], file = path_for_file)


}


