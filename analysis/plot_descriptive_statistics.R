################################################################################
### In this file I:
### - generate plots of descriptive statistics for temperature, foehn and hospitalizations
### - I focus on trends, seasonality and distribution (histogram, boxplots)
### - mean, standard deviation, median
### - these investigations will be done per station (Chur, Magadino, Lugano)


## Preamble ####
rm(list=ls())

# packages
library(dplyr); library(tidyr); library(plotly)

# load data
# here for buffer of 5 km for the three stations
data <- read.csv("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/MedStat_aggregated/CLM_buffersize_5000.csv")

# get date format
data$date = as.Date(data$date)

# remove rows with NAs
# data <- data[complete.cases(data)]
# TODO

chur = data[data$station=="Chur",]
maga = data[data$station=="Magadino",]
luga = data[data$station=="Lugano",]


# get a list of the three station names for later looping
station_names = unique(data$station)
######


## TEMPERATURE ####

# trend
plot(chur$date, chur$temp,
     cex = 0.5,
     col = 2,
     ylab = "Daily Mean Temperature [°C]",
     xlab = "Date",
     ylim = c(-15,35),
     main = paste0(chur$station[1], " Temperature Series"))

plot(luga$date, luga$temp,
     cex = 0.5,
     col = 2,
     ylab = "Daily Mean Temperature [°C]",
     xlab = "Date",
     ylim = c(-5,30),
     main = paste0(luga$station[1], " Temperature Series"))

plot(maga$date, maga$temp,
     cex = 0.5,
     col = 2,
     ylab = "Daily Mean Temperature [°C]",
     xlab = "Date",
     ylim = c(-10,30),
     main = paste0(maga$station[1], " Temperature Series"))

print(paste0(chur$station[1], " Start"))

#####
