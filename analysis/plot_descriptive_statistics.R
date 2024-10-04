################################################################################
### In this file I:
### - generate plots of descriptive statistics for temperature, foehn and hospitalizations
### - I focus on trends, seasonality and distribution (histogram, boxplots)
### - mean, standard deviation, median
### - these investigations will be done combined, because doing them station wise doesnt add anything as they will be analysed together
### - color scheme: temp: rgb(1,.1,.1,0.6), foehn: rgb(0.1,0.1,1,0.7) , MA : rgb(0.2,0.2,0.2,0.8)


# table notes for study site section
# data availability time frames for the three stations, for temp, foehn, hosp...
# canton, elevation, south/north of the alps, rows with NAs


## Preamble ####
rm(list=ls())

# packages
library(dplyr); library(tidyr); library(plotly);library(zoo)

# load data
# here for buffer of 5 km for the three stations
data <- read.csv("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MedStat_aggregated/CLM_buffersize_5000.csv")

# get date format
data$date = as.Date(data$date)

# remove rows with NAs
data <- data[complete.cases(data),]

# create two age groups below and above 64
data$a64_below = data$a014y + data$a1564y
data$a64_above = data$a6574y + data$a7584y + data$a85plusy

# daily means
data_daily_mean = data %>%
  mutate(daymonth = format(date, "%m-%d")) %>%
  group_by(daymonth) %>%
  summarise(
    across(c(all, mal, fem, a64_below, a64_above, cvd, resp, temp), mean),
    mean_f_id = mean(f_id),
    p90_f_id = quantile(f_id, 0.9) # Calculate 90th percentile
  ) %>%
  mutate(daymonth = as.Date(paste0("2000-", daymonth)))


# per stations sets
chur = data[data$station=="Chur",]
maga = data[data$station=="Magadino",]
luga = data[data$station=="Lugano",]

######


### TEMPERATURE #####
# table notes:
# per station: trend?, mean, sd, p01, p99
## some stats about the stations for a table
print(paste0(chur$station[1], " start: ", chur$date[1], ", end: ", chur$date[nrow(chur)],
             ", mean temperature: ", round(mean(chur$temp),2), ", standard deviation: ", round(sd(chur$temp),2) ))
print(paste0(maga$station[1], " start: ", maga$date[1], ", end: ", maga$date[nrow(maga)],
             ", mean temperature: ", round(mean(maga$temp),2), ", standard deviation: ", round(sd(maga$temp),2) ))
print(paste0(luga$station[1], " start: ", luga$date[1], ", end: ", luga$date[nrow(luga)],
             ", mean temperature: ", round(mean(luga$temp),2), ", standard deviation: ", round(sd(luga$temp),2) ))


## seasonality
# ticks for ploting
monthly_ticks <- data_daily_mean$daymonth[!duplicated(format(data_daily_mean$daymonth, "%Y-%m"))]
par(mfrow=c(1,2))
# plot
plot(data_daily_mean$daymonth, data_daily_mean$temp, xaxt = "n", col = rgb(1,.1,.1,0.6), ylim = c(0,25),
     xlab = "month", ylab = "daily mean temperature [°C]", pch = 20, cex.axis = 0.6, cex.label = 0.7, bty = "n")
axis(1, at = monthly_ticks + 15, labels = substr(format(monthly_ticks, "%b"),1,1), cex.axis = 0.6)


## distribution
hist(data$temp, breaks = 20, col = rgb(1,.1,.1,0.6),
     xlim = c(-15,33), # ylim = c(0,1000),
     xlab = "daily mean temperature [°C]", main = "", cex.axis = 0.6,
     ylab = "frequency")

#####


### FOEHN ####
# table notes
# per station: trend?, mean yearly, sd, mean daily, sd, 0s,
## some stats about the stations for a table
print(paste0(chur$station[1], ", mean daily foehn: ", round(mean(chur$f_id),2), ", standard deviation: ", round(sd(chur$f_id),2) ))
print(paste0(maga$station[1], ", mean daily foehn: ", round(mean(maga$f_id),2), ", standard deviation: ", round(sd(maga$f_id),2) ))
print(paste0(luga$station[1], ", mean daily foehn: ", round(mean(luga$f_id),2), ", standard deviation: ", round(sd(luga$f_id),2) ))


# seasonality
par(mfrow=c(1,2))
plot(data_daily_mean$daymonth, data_daily_mean$mean_f_id, xaxt = "n", col = rgb(0.1,0.1,1,0.7),
     xlab = "month", ylab = "daily foehn score", ylim = c(0,130), pch = 20, cex.axis = 0.6, bty = "n")
axis(1, at = monthly_ticks + 15, labels = substr(format(monthly_ticks, "%b"),1,1), cex.axis = 0.6)

# data_daily_mean$MA365 = rollmean(data_daily_mean$mean_f_id, k = 30, fill = NA, align = "center")
# lines(data_daily_mean$daymonth, data_daily_mean$MA365, col = rgb(0.2,0.2,0.2,0.8), lwd = 3.5)

data_daily_mean$MA365 = rollmean(data_daily_mean$p90_f_id, k = 30, fill = NA, align = "center")
lines(data_daily_mean$daymonth, data_daily_mean$MA365, col = rgb(0.2,0.2,.2,0.8), lwd = 3, lty = 1)
# legend("topleft",  legend = "30-day MA mean", bty = "n",
#        col = rgb(0.2,0.2,0.2,0.8),lwd = 3, , lty = 1)
legend(monthly_ticks[6], 120,  legend = "30-day MA p90", bty = "n",
       col = rgb(0.2,0.2,.2,0.7),lwd = 3, lty = 1, cex = 0.7)






## distribution
hist(data$f_id[data$f_id!=0], breaks = 20, col = rgb(0.1,0.1,1,0.7),
     xlim = c(0,300),  #ylim = c(0,900),
     xlab = "daily foehn score", main = "", ylab = "frequency", cex.axis = 0.6)

print(paste0("Percentage of 0 foehn score of all days: ", round(length(data$f_id[data$f_id==0]) / nrow(data),4) *100, "%"      ))






#####


### HOSPITALIZTAIONS ####
# table notes per station
# trend?, total, 0, mean yearly all, sd, mean daily, sd,
# percentages for of total: mal fem age64_below age64_above, and 2 disease groups


# TODO
# h
## seasonality
par(mfrow=c(1,2))

# all
plot(data_daily_mean$daymonth, data_daily_mean$all, xaxt = "n", col = rgb(.1,.8,.1,.6), xlab = "month", pch=20,
     ylab = "daily hospitalizations", main = "", ylim = c(0,12), cex.axis = 0.6, bty = "n")
axis(1, at = monthly_ticks + 15, labels = substr(format(monthly_ticks, "%b"),1,1), cex.axis = 0.6)
# data_daily_mean$MA365 = rollmean(data_daily_mean$all, k = 30, fill = NA, align = "center")
# lines(data_daily_mean$daymonth, data_daily_mean$MA365, col = rgb(0.2,0.2,0.7,0.7), lwd = 3)

# disease groups ! you need to mention that these are not all groups! (csd+resp!=100%)
# cvd
points(data_daily_mean$daymonth, data_daily_mean$cvd, col = rgb(.1,.4,.2,.6), type = "p", pch = 20)
# data_daily_mean$MA365 = rollmean(data_daily_mean$cvd, k = 30, fill = NA, align = "center")
# lines(data_daily_mean$daymonth, data_daily_mean$MA365, col = rgb(0.2,0.2,0.7,0.7), lwd = 3)

points(data_daily_mean$daymonth, data_daily_mean$resp, col = rgb(.9,.9,.1,.5),  type = "p", pch = 20)
# data_daily_mean$MA365 = rollmean(data_daily_mean$resp, k = 30, fill = NA, align = "center")
# lines(data_daily_mean$daymonth, data_daily_mean$MA365, col = rgb(0.2,0.7,0.2,0.7), lwd = 3)

legend(monthly_ticks[1], 12,  legend = "all", bty = "n",
       col = rgb(.1,.8,.1,.6) , cex = 0.7, pch = 20, pt.cex = 1.5)
legend(monthly_ticks[4], 12,  legend = "cvd", bty = "n",
       col = rgb(.1,.4,.2,.6) , cex = 0.7, pch = 20, pt.cex = 1.5)
legend(monthly_ticks[7], 12,  legend = "resp", bty = "n",
       col = rgb(.9,.9,.1,.5) , cex = 0.7, pch = 20, pt.cex = 1.5)


## distribution
hist(data$all, breaks = 40, col = rgb(.1,.8,.1,.6),
     xlim = c(0,25),  ylim = c(0,4000),
     xlab = "daily hospitalizations", main = "", ylab = "frequency",
     cex.axis = 0.6)
#######














# ## TEMPERATURE station wise (non useful anymore) ####
#
# # trend with yearly moving average
# plot(chur$date, chur$temp,
#      cex = 0.5,
#      col = 2,
#      ylab = "Daily Mean Temperature [°C]",
#      xlab = "Date",
#      ylim = c(-15,35),
#      main = paste0(chur$station[1], " Temperature Series"))
# chur$MA365 = rollmean(chur$temp, k = 365, fill = NA, align = "center")
# lines(chur$date, chur$MA365, col = rgb(0.2,0.2,0.2,0.7), lwd = 5)
# grid()
# legend("bottomright",  legend = "365-day MA", bty = "n", col = rgb(0.2,0.2,0.2,0.7), lwd = 5)
#
# plot(luga$date, luga$temp,
#      cex = 0.5,
#      col = 2,
#      ylab = "Daily Mean Temperature [°C]",
#      xlab = "Date",
#      ylim = c(-5,30),
#      main = paste0(luga$station[1], " Temperature Series"))
# luga$MA365 = rollmean(luga$temp, k = 365, fill = NA, align = "center")
# lines(luga$date, luga$MA365, col = rgb(0.2,0.2,0.2,0.7), lwd = 5)
# grid()
# legend("bottomright",  legend = "365-day MA", bty = "n", col = rgb(0.2,0.2,0.2,0.7), lwd = 5)
#
# plot(maga$date, maga$temp,
#      cex = 0.5,
#      col = 2,
#      ylab = "Daily Mean Temperature [°C]",
#      xlab = "Date",
#      ylim = c(-10,30),
#      main = paste0(maga$station[1], " Temperature Series"))
# maga$MA365 = rollmean(maga$temp, k = 365, fill = NA, align = "center")
# lines(maga$date, maga$MA365, col = rgb(0.2,0.2,0.2,0.7), lwd = 5)
# grid()
# legend("bottomright",  legend = "365-day MA", bty = "n", col = rgb(0.2,0.2,0.2,0.7), lwd = 5)
#
#
#
# # seasonality
# chur_d_mean <- chur %>%
#   mutate(daymonth = format(date, "%m-%d")) %>%
#   group_by(daymonth) %>%
#   summarise(mean_temp = mean(temp, na.rm = TRUE)) %>%
#   mutate(daymonth = as.Date(paste0("2000-", daymonth)))
#
# luga_d_mean <- luga %>%
#   mutate(daymonth = format(date, "%m-%d")) %>%
#   group_by(daymonth) %>%
#   summarise(mean_temp = mean(temp, na.rm = TRUE)) %>%
#   mutate(daymonth = as.Date(paste0("2000-", daymonth)))
#
# maga_d_mean <- maga %>%
#   mutate(daymonth = format(date, "%m-%d")) %>%
#   group_by(daymonth) %>%
#   summarise(mean_temp = mean(temp, na.rm = TRUE)) %>%
#   mutate(daymonth = as.Date(paste0("2000-", daymonth)))
#
# # monthly ticks
# monthly_ticks <- chur_d_mean$daymonth[!duplicated(format(chur_d_mean$daymonth, "%Y-%m"))]
#
# plot(chur_d_mean$daymonth, chur_d_mean$mean_temp, xaxt = "n", col = 2,
#      xlab = "Time", ylab = "Daily Mean Temperature [°C]", main = "Chur")
# axis(1, at = monthly_ticks, labels = substr(format(monthly_ticks, "%b"),1,1))
#
# plot(chur_d_mean$daymonth, luga_d_mean$mean_temp, xaxt = "n", col = 2,
#      xlab = "Time", ylab = "Daily Mean Temperature [°C]", main = "Lugano")
# axis(1, at = monthly_ticks, labels = substr(format(monthly_ticks, "%b"),1,1))
#
# plot(chur_d_mean$daymonth, maga_d_mean$mean_temp, xaxt = "n", col = 2,
#      xlab = "Time", ylab = "Daily Mean Temperature [°C]", main = "Magadino")
# axis(1, at = monthly_ticks, labels = substr(format(monthly_ticks, "%b"),1,1))
#
# # distribution
# hist(chur$temp, breaks = 40, freq = FALSE, col = 2,
#      xlim = c(-15,35), ylim = c(0,0.05),
#      xlab = "Daily Mean Temperature [°C]", main = "Chur")
# chur_mean = round(mean(chur$temp),2)
# chur_median = round(median(chur$temp),2)
# text(-12,0.045,labels = paste0("mean: ", chur_mean), pos = 4)
# text(-12,0.04,labels = paste0("median: ", chur_median), pos = 4)
#
# hist(luga$temp, breaks = 40, freq = FALSE, col = 2,
#      xlim = c(-10,30), ylim = c(0,0.05),
#      xlab = "Daily Mean Temperature [°C]", main = "Lugano")
# luga_mean = round(mean(luga$temp),2)
# luga_median = round(median(luga$temp),2)
# text(-9,0.04,labels = paste0("mean: ", luga_mean), pos = 4)
# text(-9,0.035,labels = paste0("median: ", luga_median), pos = 4)
#
# hist(maga$temp, breaks = 40, freq = FALSE, col = 2,
#      xlim = c(-10,30), ylim = c(0,0.05),
#      xlab = "Daily Mean Temperature [°C]", main = "Magadino")
# maga_mean = round(mean(maga$temp),2)
# maga_median = round(median(maga$temp),2)
# text(-9,0.045,labels = paste0("mean: ", maga_mean), pos = 4)
# text(-9,0.04,labels = paste0("median: ", maga_median), pos = 4)
# #####
