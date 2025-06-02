################################################################################
### DESCRIPTIVE STATISTICS tables
################################################################################
# - contribution and distribution of hospitalizations, foehn wind and temperature
# - description of the Medstat regions



### Data
#----

rm(list=ls())

# libraries
library(dplyr); library(tidyr); library(plotly);library(zoo); library(viridis);
library(knitr);library(kableExtra);library(webshot2); library(magick)

# colors
colors <- viridis(3, option = "viridis")

# data
data = read.csv("../data/Medstat_hospitalizations_aggregated/hosp_buffer_8000.csv") |>
  mutate(date = as.Date(date),
         stratum_dow = as.factor(stratum_dow),
         y64 = a014y + a1564y,
         o64 = a6574y + a7584y + a85plusy,
         f_id_binary = ifelse(f_id >= 72, 1, 0),
         station_year = paste0(station, year))


# index to include only stratum that have hosp counts
ind_dow = tapply(data$all, data$stratum_dow, sum)

# MMT function
source("functions/findmin.R")

# define the maximum lag distance we account for
maxlago <- 3

# time series of averages across all stations
data_time_serie = data |>
  group_by(date) |>
  summarise(

    # averaging over stations
    mean_f_id = mean(f_id, na.rm = TRUE),
    mean_temp = mean(temp, na.rm = TRUE),
    mean_all = mean(all, na.rm = TRUE)) |>

  mutate(

    # moving averages over 30 days
    MA_f_id = rollmean(mean_f_id, k = 30, fill = NA, align = "center"),
    MA_all = rollmean(mean_all, k = 30, fill = NA, align = "center"),
    MA_temp = rollmean(mean_temp, k = 30, fill = NA, align = "center")
  )


# start dates of station records
start_dates = data.frame(station = as.vector(unique(data$station)),
                         date = rep("no",8))

for (i  in 1:nrow(start_dates)) {
  subs = data[data$station == as.character(unique(data$station)[i]),]
  start_dates$date[i] <- as.character(min(subs$date))
}

start_dates$date <- as.Date(start_dates$date)


# daily means
data_daily_mean = data |>
  mutate(daymonth = format(date, "%m-%d")) |>
  group_by(daymonth) |>
  summarise(

    # averages
    across(c(mal, fem, y64, o64, cvd, resp), mean),
    mean_f_id = mean(f_id),
    mean_temp = mean(temp),
    mean_all = mean(all),

    # percentiles
    p50_f_id = quantile(f_id, 0.5),
    p95_f_id = quantile(f_id, 0.95),
    p75_f_id = quantile(f_id, 0.8),
    p90_f_id = quantile(f_id, 0.9),
    p25_temp = quantile(temp, 0.25),
    p50_temp = quantile(temp, 0.5),
    p75_temp = quantile(temp, 0.75),
    p95_temp = quantile(temp, .95),
    p99_temp = quantile(temp, .95),
    p25_all = quantile(all, 0.25),
    p50_all = quantile(all, 0.5),
    p75_all = quantile(all, 0.75),
  ) |>
  mutate(
    daymonth = as.Date(paste0("2000-", daymonth)),

    # moving averages of percentiles
    MA50 = rollmean(p50_f_id, k = 30, fill = NA, align = "center"),
    MA95 = rollmean(p95_f_id, k = 30, fill = NA, align = "center"),
    MA75 = rollmean(p75_f_id, k = 30, fill = NA, align = "center"),
    MA90 = rollmean(p90_f_id, k = 30, fill = NA, align = "center"),

    MA25_temp = rollmean(p25_temp, k = 30, fill = NA, align = "center"),
    MA50_temp = rollmean(p50_temp, k = 30, fill = NA, align = "center"),
    MA75_temp = rollmean(p75_temp, k = 30, fill = NA, align = "center"),
    MA95_temp = rollmean(p95_temp, k = 30, fill = NA, align = "center"),
    MA99_temp = rollmean(p99_temp, k = 30, fill = NA, align = "center"),
    MA25_all = rollmean(p25_all, k = 30, fill = NA, align = "center"),
    MA50_all = rollmean(p50_all, k = 30, fill = NA, align = "center"),
    MA75_all = rollmean(p75_all, k = 30, fill = NA, align = "center")
         )

# ticks annual cycle x axis
monthly_ticks <- data_daily_mean$daymonth[!duplicated(format(data_daily_mean$daymonth, "%Y-%m"))]


# for Medstat table
station_MedStat = read.csv("data_nonsensitive/Medstat_regions_per_buffer_size/MDSTID_MetRegions_8000_per_station.csv", header = TRUE)

station_MedStat_full <- rbind(c("ALT","CHU","DAV","LUG","MAG","MVE","ROB","VIS"), #Abbreviation
                              c(46.887042,46.870622,46.812956,46.003833,46.160019,46.298772,46.347194,46.302875),  #x coord
                              c(8.621806,9.530814,9.843492,8.960117,8.933608,7.460761,10.062964,7.842903),  #y coord
                              station_MedStat)
rownames(station_MedStat_full) <- c("Abbreviation", "X Coordinate", "Y Coordinate",
                                    "1 MedStat regions", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")


# for station tables
stations = as.vector(unique(data$station))

#----


### Table 1
#----

# Descriptive statistics for the whole data set and for every station. Temperature
# statistics are based on annual means in °C. Foehn statistics are based on annual
# sums and the percentage indicates the ratio of the sum of foehn per station
# compared to the sum of foehn from all stations. The binary foehn percentage shows the
# percentage of days with foehn with a binary threshold of 72 which corresponds to 6 h of
# full foehn. Hospitalization characteristics are based on annual sums. The percentage of
# hospitalizations indicates the contribution of hospitalizations per station to the whole data set.
# NA corresponds to the missing days. Standard deviations are given in brackets.

# empty data frame to be filled
gen_stats = data.frame(period = rep(NA,9),
                       years = rep(NA,9),
                       temp = rep(NA, 9),
                       foehn = rep(NA,9),
                       per_tot_foehn = rep(NA,9),
                       foehn_bin_perc = rep(NA,9),
                       hosp = rep(NA,9),
                       per_tot_hosp = rep(NA,9),
                       row.names = c("all stations",stations))

tot_foehn = sum(data$f_id)
tot_hosp = sum(data$all)

# yearly mean data weather
subs_y_all_weather = data |>
  group_by(station_year) |>
  summarize(mittel_temp = mean(temp, na.rm = TRUE),
            summe_foehn = sum(f_id, na.rm = TRUE))

# yearly mean data hosp
subs_y_all_hosp = data |>
  group_by(station_year) |>
  summarize(summe = sum(all, na.rm = TRUE))

# for all stations
gen_stats$period[1] = as.character(min(data$date))
gen_stats$years[1] = length(unique(data$station_year))
gen_stats$temp[1] = paste0( sprintf("%.1f",round(mean(subs_y_all_weather$mittel_temp),digits=1)), " [", sprintf("%.1f",round(sd(subs_y_all_weather$mittel_temp), digits= 1)), "]")
gen_stats$foehn[1] = paste0( round(mean(subs_y_all_weather$summe_foehn)/12,digits=0), " [", round(sd(subs_y_all_weather$summe_foehn)/12, digits= 0), "]")
gen_stats$per_tot_foehn[1] = "100 %"
gen_stats$foehn_bin_perc[1] = paste0(sprintf("%.1f",round(length(which(data$f_id_binary>0))/nrow(data)*100, digits = 3)), " %")
gen_stats$hosp[1] = paste0( round(mean(subs_y_all_hosp$summe),digits=0), " [", round(sd(subs_y_all_hosp$summe), digits= 0), "]")
gen_stats$per_tot_hosp[1] = "100 %"

# loop through stations
for (i in 1:length(stations)) {
  sta = stations[i]

  subs = data[data$station == sta,]

  # yearly mean data weather
  subs_y_all_weather = subs |>
    group_by(station_year) |>
    summarize(mittel_temp = mean(temp, na.rm = TRUE),
              summe_foehn = sum(f_id, na.rm = TRUE))

  # yearly mean data hosp
  subs_y_all_hosp = subs |>
    group_by(station_year) |>
    summarize(summe = sum(all, na.rm = TRUE))

  # populate the table
  gen_stats$period[i+1] = as.character(min(subs$date))
  gen_stats$years[i+1] = length(unique(subs$station_year))
  gen_stats$temp[i+1] = paste0( sprintf("%.1f",round(mean(subs_y_all_weather$mittel_temp),digits=1)), " [", sprintf("%.1f",round(sd(subs_y_all_weather$mittel_temp), digits= 1)), "]")
  gen_stats$foehn[i+1] = paste0( round(mean(subs_y_all_weather$summe_foehn)/12,digits=0), " [", round(sd(subs_y_all_weather$summe_foehn)/12, digits= 0), "]")
  gen_stats$per_tot_foehn[i+1] =  paste0( sprintf("%.1f",round( sum(subs$f_id) / tot_foehn*100, digits=1)), " %")
  gen_stats$foehn_bin_perc[i+1] = paste0( sprintf("%.1f",round(length(which(subs$f_id_binary>0))/nrow(subs)*100, digits = 1)) , " %")
  gen_stats$hosp[i+1] = paste0( round(mean(subs_y_all_hosp$summe),digits=0), " [", round(sd(subs_y_all_hosp$summe), digits= 0), "]")
  gen_stats$per_tot_hosp[i+1] = paste0(sprintf("%.1f",round(sum(subs$all) / tot_hosp *100, digits=1)) , " %")
}

# save the table
write.csv(gen_stats, file = "data_nonsensitive/descriptive_statistics_tables/Table1.csv")

#----


### sTable 1
#----

# Abbreviation and location of meteorological measurement stations and their
# assigned MedStat regions.

# save the table
write.csv(station_MedStat_full, file = "data_nonsensitive/descriptive_statistics_tables/sTable1.csv")

#----


### sTable 3
#----

# Contribution in percent of different subpopulations to all-cause hospitalizations
# Area counts across all stations and for every station.

# empty data frame
subgroups = data.frame(total = rep(NA,9),
                       mal = rep(NA, 9),
                       row.names = c("all stations", stations))

# across all stations
subgroups$total[1] = sum(data$all)
subgroups$mal[1] = paste( sprintf("%.1f",round( sum(data$mal) / sum(data$all)*100, digits = 1)) , " %")
subgroups$fem[1] = paste( sprintf("%.1f",round( sum(data$fem) / sum(data$all)*100, digits = 1)), " %")

subgroups$y64[1] = paste( sprintf("%.1f",round( sum(data$y64) / sum(data$all)*100, digits = 1)), " %")
subgroups$o64[1] = paste( sprintf("%.1f",round( sum(data$o64) / sum(data$all)*100, digits = 1)), " %")

subgroups$cvd[1] = paste( sprintf("%.1f",round( sum(data$cvd) / sum(data$all)*100, digits = 1) ), " %")
subgroups$resp[1] = paste( sprintf("%.1f",round( sum(data$resp) / sum(data$all)*100, digits =1)), " %")
subgroups$inf[1] = paste( sprintf("%.1f",round( sum(data$inf) / sum(data$all)*100, digits = 1) ), " %")
subgroups$uri[1] = paste( sprintf("%.1f",round( sum(data$uri) / sum(data$all)*100, digits = 1) ), " %")
subgroups$ment[1] = paste( sprintf("%.1f",round( sum(data$ment) / sum(data$all)*100, digits = 1)), " %")


# for every individual station
for (i in 1:length(unique(data$station))) {
  subs = data[data$station == unique(data$station)[i],]

  subgroups$total[i+1] = sum(subs$all)

  subgroups$mal[i+1] = paste( sprintf("%.1f",round( sum(subs$mal) / sum(subs$all)*100, digits =1)) , " %")
  subgroups$fem[i+1] = paste( sprintf("%.1f",round( sum(subs$fem) / sum(subs$all)*100, digits = 1)) , " %")

  subgroups$y64[i+1] = paste( sprintf("%.1f",round( sum(subs$y64) / sum(subs$all)*100, digits = 1)) , " %")
  subgroups$o64[i+1] = paste( sprintf("%.1f",round( sum(subs$o64) / sum(subs$all)*100, digits = 1)) , " %")

  subgroups$cvd[i+1] = paste( sprintf("%.1f",round( sum(subs$cvd) / sum(subs$all)*100, digits = 1)) , " %")
  subgroups$resp[i+1] = paste( sprintf("%.1f",round( sum(subs$resp) / sum(subs$all)*100, digits =1)), " %")
  subgroups$inf[i+1] = paste( sprintf("%.1f",round( sum(subs$inf) / sum(subs$all)*100, digits = 1) ), " %")
  subgroups$uri[i+1] = paste( sprintf("%.1f",round( sum(subs$uri) / sum(subs$all)*100, digits = 1) ), " %")
  subgroups$ment[i+1] = paste( sprintf("%.1f",round( sum(subs$ment) / sum(subs$all)*100, digits = 1)), " %")


}

# assign full column names
colnames(subgroups) <- c("all", "male", "female", "<65 years", ">64 years", "cardiovascular", "respiratory", "infectious", "genitourinary", "mental")

# save the table
write.csv(subgroups, file = "data_nonsensitive/descriptive_statistics_tables/sTable3.csv")

#----

