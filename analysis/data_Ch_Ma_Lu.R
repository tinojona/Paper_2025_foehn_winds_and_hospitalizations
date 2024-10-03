################################################################################
### In this file I:
### - extract the aggregated station region data for multiple buffers
### - delete all rows are not station Chur, Magadino or Lugano
### - and save this data in the same folder again

rm(list=ls())

# get all the files in that folder
files_buffer = list.files("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MedStat_aggregated")

# filter for ones that have hosp_buffer in them
files_buffer = grep("hosp_buffer", files_buffer, value = TRUE)

buffers = c(1000,10000,2500,5000)

# start loop through these files
for(i in 1:4){
  file = files_buffer[i]
  buffer_size = buffers[i]

  data = read.csv(paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MedStat_aggregated/", file))

  subdata = data[data$station == "Chur" | data$station == "Lugano" | data$station == "Magadino",]

  # save the new data set
  write.csv(subdata, file = paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MedStat_aggregated/CLM_buffersize_", buffer_size, ".csv"))

}
