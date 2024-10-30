################################################################################
# GNM Temperature


rm(list=ls())


### PACKAGES ####
library(dlnm);library(splines);library(ggplot2);library(viridis);library(gnm)
library(dplyr)
# mmt function
source("functions/findmin.R")
######

### DATA ####

# buffer size for data file read
buffer = 8000

# data
data = read.csv(paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MedStat_aggregated/centroid_aggregated/hosp_buffer_", buffer, ".csv"))

# reasign original types
data$date = as.Date(data$date)
data$station <- as.factor(data$station)


# index to include only stratum that have hosp counts
data$stratum_dow = as.factor(data$stratum_dow); data$stratum = as.factor(data$stratum)
ind_dow = tapply(data$all, data$stratum_dow, sum); ind = tapply(data$all, data$stratum, sum)




# age categories
data <- data %>%
  mutate(y64 = a014y + a1564y) %>%
  mutate(o64 = a6574y + a7584y + a85plusy)

# define the maximum lag distance we account for
maxlago <- 3

#####

### CROSSBASIS TEMPERATURE ####

cb.temp <- crossbasis(data$temp,
                      lag=21,
                      argvar=list(fun="ns", knots = quantile(data$temp, c(.5,.9), na.rm=TRUE)),
                      arglag=list(fun="ns", knots = logknots(21,3)),
                      group = data$station)

#####


# model
mod <- gnm(all ~ cb.temp, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

# prediction preliminary
pred <- crosspred(cb.temp, mod, cumul=FALSE, cen = 0)

# find min value for plotting
min1 <- findmin(cb.temp,pred,from=quantile(data$temp, .1),to=quantile(data$temp, .9))

# prediction final
prednew <- crosspred(cb.temp, mod, cumul=FALSE, cen = min1)



plot(prednew,              ## cumulative exposure
     "overall",
     col = 2,
     ci.arg = list(density = 20, col = 2 ,angle = 45),
     xlab = "Temperature",
     ylab = "Cumulative Response",
     lwd = 2,
     main = "")
