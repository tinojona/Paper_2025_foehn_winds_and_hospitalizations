################################################################################
### DESCRIPTIVE STATISTICS
################################################################################
# - plots of distribution, seasonality, time seires for temperature, foehn wind and all-cause hospitalizations



### Data
#----

rm(list=ls())

# libraries
library(dplyr); library(tidyr); library(plotly);library(zoo); library(viridis);
library(knitr);library(kableExtra);library(webshot2); library(magick)

# colors
colors <- viridis(3, option = "viridis")

# data
data = read.csv("/Volumes/FS/_ISPM/CCH/Tino/master_thesis/data/Medstat_hospitalizations_aggregated/hosp_buffer_8000.csv") |>
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

### sFigure 2
#----

# (a) daily foehn score distribution excluding 0-foehn days, (b) daily mean foehn wind
# score with 30-day moving averages of the 50th, 75th, 90th percentile, (c) daily mean
# temperature distribution, (d) daily averages of daily mean temperature with 30-day moving
# averages of the 25th, 50th, 75th percentile, (e) daily all-cause hospitalization distribution,
# (f) daily mean all-cause hospitalizations with 30-day moving averages of the 25th, 50th, 75th
# percentile.


png("output/figures/sFigure2.png", width = 1800, height = 2300, res = 300)

par(
  # mfcol =c(2,3),
    mfrow = c(3,2),
    mar = c(4,4,.5,1),
    mgp = c(2, .5, 0))


hist(data$f_id[data$f_id!=0],
     breaks = 40,
     col = colors[2],
     xlim = c(0,300),
     xlab = "foehn wind",
     main = "",
     # cex.axis = 0.6,
     ylab = "frequency")
text(285, 2250, labels = "(a)", pos = 2)

plot(data_daily_mean$daymonth,
     data_daily_mean$mean_f_id,
     xaxt = "n", col = colors[2],
     type = "p",
     xlab = "month",
     ylab = "foehn wind",
     main = "",
     ylim = c(0,140),
     pch = 16,
     cex = .6,
     # cex.axis = 0.6,
     # lwd = 2,
     bty = "n")

axis(1,
     at = monthly_ticks + 15,
     labels = substr(format(monthly_ticks, "%b"),1,1)#,cex.axis = 0.6
)

lines(data_daily_mean$daymonth, data_daily_mean$MA50, col = 1, lwd = 2, lty = 1)
lines(data_daily_mean$daymonth, data_daily_mean$MA75, col = 1, lwd = 2, lty = 2)
lines(data_daily_mean$daymonth, data_daily_mean$MA90, col = 1, lwd = 2, lty = 3)
# lines(data_daily_mean$daymonth, data_daily_mean$MA95, col = 5, lwd = 2, lty = 1)

legend(monthly_ticks[6]+15, 120,
       legend = c("p50", "p75","p90"), bty = "n",
       col = c(1,1,1) ,lwd = 1, lty = c(1,2,3), cex = 0.8,
       ncol = 2)

text(monthly_ticks[12], 135,
     labels = "(b)", pos = 2)

hist(data$temp,
     breaks = 30,
     col = colors[1],
     # xlim = c(0,300),
     xlab = "temperature [\u00B0C]",
     main = "",
     # cex.axis = 0.6,
     ylab = "frequency")
text(28, 3950, labels = "(c)", pos = 2)

plot(data_daily_mean$daymonth,
     data_daily_mean$mean_temp,
     xaxt = "n", col = colors[1],
     type = "p",
     xlab = "month",
     ylab = "temperature [\u00B0C]",
     main = "",
     ylim = c(-5,25),
     pch = 16,
     cex = .6,
     # cex.axis = 0.6,
     bty = "n")

axis(1,
     at = monthly_ticks + 15,
     labels = substr(format(monthly_ticks, "%b"),1,1)#, cex.axis = 0.6
)

lines(data_daily_mean$daymonth, data_daily_mean$MA25_temp, col = 1, lwd = 2, lty = 3)
lines(data_daily_mean$daymonth, data_daily_mean$MA50_temp, col = 1, lwd = 2, lty = 1)
lines(data_daily_mean$daymonth, data_daily_mean$MA75_temp, col = 1, lwd = 2, lty = 2)
# lines(data_daily_mean$daymonth, data_daily_mean$MA95_temp, col = 1, lwd = 2, lty = 2)
# lines(data_daily_mean$daymonth, data_daily_mean$MA99_temp, col = 2, lwd = 2, lty = 1)

legend(monthly_ticks[2]+15, -2,
       legend = c("p25", "p50","p75"), bty = "n",
       col = c(1,1,1),lwd = 1, lty = c(3,1,2), cex = 0.8,
       ncol = 3)

text(monthly_ticks[12], 23, labels = "(d)", pos = 2)


hist(data$all,
     breaks = 50,
     col = colors[3],
     xlim = c(0,20),
     xlab = "hospitalizations",
     main = "",
     # cex.axis = 0.6,
     ylab = "frequency")
text(19, 21000, labels = "(e)", pos = 2)

plot(data_daily_mean$daymonth,
     data_daily_mean$mean_all,
     xaxt = "n", col = colors[3],
     type = "p",
     xlab = "month",
     ylab = "hospitalizations",
     main = "",
     ylim = c(0,6),
     pch = 16,
     cex = .6,
     # cex.axis = 0.6,
     bty = "n")

axis(1,
     at = monthly_ticks + 15,
     labels = substr(format(monthly_ticks, "%b"),1,1)#, cex.axis = 0.6
)

lines(data_daily_mean$daymonth, data_daily_mean$MA25_all, col = 1, lwd = 2, lty = 3)
lines(data_daily_mean$daymonth, data_daily_mean$MA50_all, col = 1, lwd = 2, lty = 1)
lines(data_daily_mean$daymonth, data_daily_mean$MA75_all, col = 1, lwd = 2, lty = 2)

legend(monthly_ticks[3]+5, 5.75,
       legend = c("p25", "p50","p75"), bty = "n",
       col = c(1,1,1),lwd = 1, lty = c(3,1,2), cex = 0.8,
       ncol = 3)

text(monthly_ticks[12], 6, labels = "(f)", pos = 2)

dev.off()




#----

### sFigure 3
#----

# (a) 30-day moving average time series of mean daily foehn, (b) 30-day moving average
# time series of average daily mean temperature, (c) 30-day moving average of mean
# daily all-cause hospitalization. (1) start of data from Altdorf, Chur, Davos, Montana,
# (2) start of data from Magadino, (3) start of data from Lugano, (4) start of data from Visp,
# (5) start of data from Poschiavo. (x) redefinition of the MedStat regions by the Swiss
# Federal Office for Statistics.

png("output/figures/sFigure3.png", width = 1600, height = 1400, res = 300)

# margins
par(mfrow = c(3, 1),
    mar = c(0, 4, 0, 0.5),  # Top two plots have zero bottom margin
    oma = c(4, 0, 0, 0),
    mgp = c(2, .5, 0))      # Outer bottom margin for x-axis label

# foehn plot
plot(data_time_serie$date, data_time_serie$MA_f_id,
     type = "n",
     xlab = "",
     ylab = "foehn wind",
     xaxt = "n"
)

for(i in start_dates$date){  abline(v = i, col = 1, lty = 2)}

abline(v = as.Date("2008-01-01"), col = "brown1", lty = 1)

lines(data_time_serie$date, data_time_serie$MA_f_id,
      lty = 1,
      lwd = 2,
      col = colors[2],)

text(data_time_serie$date[nrow(data_time_serie)], 78, labels = "(a)", pos = 2)

# temp plot
plot(data_time_serie$date, data_time_serie$MA_temp,
     type = "n",
     xlab = "",
     ylab = "temperature [\u00B0C]",
     xaxt = "n"
)

for(i in start_dates$date){  abline(v = i, col = 1, lty = 2)}

abline(v = as.Date("2008-01-01"), col = "brown1", lty = 1)

lines(data_time_serie$date, data_time_serie$MA_temp,
      lty = 1,
      lwd = 2,
      col = colors[1],)

text(data_time_serie$date[nrow(data_time_serie)], -2, labels = "(b)", pos = 2)

# hosp plot
plot(data_time_serie$date, data_time_serie$MA_all,
     type = "n",
     xlab = "",
     ylab = "hospitalizations",
)

mtext("time", side = 1, line = 2, outer = TRUE, cex = 0.8)

for(i in start_dates$date){  abline(v = i, col = 1, lty = 2)}

abline(v = as.Date("2008-01-01"), col = "brown1", lty = 1)

lines(data_time_serie$date, data_time_serie$MA_all,
      lty = 1,
      lwd = 2,
      col = colors[3])

text(data_time_serie$date[nrow(data_time_serie)], 1.8, labels = "(c)", pos = 2)

text(start_dates$date[1]+85, 4.8, labels = "(1)", pos = 2, cex = .8)
text(start_dates$date[4]+125, 4.8, labels = "(3)", pos = 2, cex = .8)
text(start_dates$date[5]+85, 4.8, labels = "(2)", pos = 2, cex = .8)
text(start_dates$date[7]-70, 4.8, labels = "(5)", pos = 4, cex = .8)
text(start_dates$date[8]+85, 4.8, labels = "(4)", pos = 2, cex = .8)
text(as.Date("2008-01-01")+70, 5.2, labels = "(x)", pos = 2, col = "brown1", cex = .8)

dev.off()

#----


