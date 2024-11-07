################################################################################
# GNM FOEHN as a modifier of temperature

rm(list=ls())

### PACKAGES ####
library(dlnm);library(splines);library(ggplot2);library(viridis);library(gnm);library(mgcv)
library(dplyr)
######


### DATA ####

# buffer size for data file read
buffer = 8000

data = read.csv(paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MedStat_aggregated/centroid_aggregated/hosp_buffer_", buffer, ".csv"))

data$date = as.Date(data$date)
data$station <- as.factor(data$station)

# index to include only stratum that have hosp counts
data$stratum_dow = as.factor(data$stratum_dow); data$stratum = as.factor(data$stratum)
ind_dow = tapply(data$all, data$stratum_dow, sum); ind = tapply(data$all, data$stratum, sum)


# create larger age groups
data <- data %>%
  mutate(y64 = a014y + a1564y) %>%
  mutate(o64 = a6574y + a7584y + a85plusy)

# define the maximum lag distance we account for
maxlago <- 3

# ONLY 2008 ONWARDS
# data <- data[data$date >= "2008-01-01", ]

# mmt function
source("functions/findmin.R")

#####

### CROSSBASIS TEMPERATURE ######
# crossbasis temp
cb.temp <- crossbasis(data$temp,
                      lag=21,
                      argvar=list(fun="ns", knots = quantile(data$temp, c(.5,.9), na.rm=TRUE)),
                      arglag=list(fun="ns", knots = logknots(21,3)),
                      group = data$station)
#####



### VISUALIZATION  ####

# different thresholds

par(mfrow=c(2,2))
for( i in c(10,35,80,140)){ #c(0,5,60,120

# binary foehn
foehn_bin <- ifelse(data$f_id > i, 0, 1)
foehn_bin_rev <- ifelse(foehn_bin == 1, 0, 1)


modif <- cb.temp * foehn_bin
modif_rev <- cb.temp *foehn_bin_rev

# model with and without foehn
mod_modif <- gnm(inf ~  cb.temp + modif, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
mod_modif_rev <- gnm(inf ~  cb.temp + modif_rev, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

# prediction with and without foehn
pred_modif  <- crosspred(cb.temp, mod_modif, cen = 20, cumul=FALSE)
pred_modif2 <- crosspred(cb.temp, mod_modif_rev, cen = 20, cumul=FALSE)

# get min value of both predictions and use for centering
min1 <- findmin(cb.temp,pred_modif,from=quantile(data$temp, .1),to=quantile(data$temp, .9))
min2 <- findmin(cb.temp,pred_modif2,from=quantile(data$temp, .1),to=quantile(data$temp, .9))

pred_modif_new     <- crosspred(cb.temp, mod_modif, cen = min1, cumul=FALSE)
pred_modif_rev_new <- crosspred(cb.temp, mod_modif_rev, cen = min2, cumul=FALSE)


plot(pred_modif_new,              ## cumulative exposure
     "overall",
     col = 2,
     ci.arg = list(density = 20, col = 2 ,angle = -45),
     lwd = 2,
     main = paste0("Overall, binary thr.=",i, ", ", as.character(mod_modif$formula[2])),
     ylim = c(0.7,3))

lines(pred_modif_rev_new,           ## cumulative exposure
     "overall",
     col = 4,
     ci = "area",
     ci.arg = list(density = 20, col = 4 ,angle = 45),
     lwd = 2)

legend("topright", legend = c("foehn days", "non foehn d."), col = c(2,4), lwd = 2)
}



#####

dev.off()

