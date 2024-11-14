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

#### different thresholds ----
par(mfrow=c(2,2))
for( i in c(10,35,72,140)){ #c(0,5,60,120

# binary foehn
foehn_bin <- ifelse(data$f_id > i, 0, 1)
foehn_bin_rev <- ifelse(foehn_bin == 1, 0, 1)


modif <- cb.temp * foehn_bin
modif_rev <- cb.temp *foehn_bin_rev

# model with and without foehn
mod_modif <- gnm(fem ~  cb.temp + modif, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
mod_modif_rev <- gnm(fem ~  cb.temp + modif_rev, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

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


#### different model types -----
# gnm without eliminate stratum
par(mfrow=c(2,2))
for( i in c(10,35,72,140)){ #c(0,5,60,120

  # binary foehn
  foehn_bin <- ifelse(data$f_id > i, 0, 1)
  foehn_bin_rev <- ifelse(foehn_bin == 1, 0, 1)

  modif <- cb.temp * foehn_bin
  modif_rev <- cb.temp *foehn_bin_rev

  # model with and without foehn
  mod_modif <- gnm(inf ~  cb.temp + modif, data = data,  family=quasipoisson())
  mod_modif_rev <- gnm(inf ~  cb.temp + modif_rev, data = data,  family=quasipoisson())

  pred_modif_new     <- crosspred(cb.temp, mod_modif, cen = 17, cumul=FALSE) # min1
  pred_modif_rev_new <- crosspred(cb.temp, mod_modif_rev, cen = 17, cumul=FALSE) # min2


  plot(pred_modif_new,              ## cumulative exposure
       "overall",
       col = 2,
       ci.arg = list(density = 20, col = 2 ,angle = -45),
       lwd = 2,
       main = paste0("Overall, binary thr.=",i, ", ", as.character(mod_modif$formula[2])),
       ylim = c(0.3,3))

  lines(pred_modif_rev_new,           ## cumulative exposure
        "overall",
        col = 4,
        ci = "area",
        ci.arg = list(density = 20, col = 4 ,angle = 45),
        lwd = 2)

  legend("topleft", legend = c("foehn days", "non foehn d."), col = c(2,4), lwd = 2)
}


# glm without eliminate
for( i in c(10,35,80,140)){ #c(0,5,60,120

  # binary foehn
  foehn_bin <- ifelse(data$f_id > i, 0, 1)
  foehn_bin_rev <- ifelse(foehn_bin == 1, 0, 1)


  modif <- cb.temp * foehn_bin
  modif_rev <- cb.temp *foehn_bin_rev

  # model with and without foehn
  mod_modif <- glm(inf ~  cb.temp + modif, data = data,  family=quasipoisson())
  mod_modif_rev <- glm(inf ~  cb.temp + modif_rev, data = data,  family=quasipoisson())

  pred_modif_new     <- crosspred(cb.temp, mod_modif, cen = 17, cumul=FALSE) # min1
  pred_modif_rev_new <- crosspred(cb.temp, mod_modif_rev, cen = 17, cumul=FALSE) # min2


  plot(pred_modif_new,              ## cumulative exposure
       "overall",
       col = 2,
       ci.arg = list(density = 20, col = 2 ,angle = -45),
       lwd = 2,
       main = paste0("Overall, binary thr.=",i, ", ", as.character(mod_modif$formula[2])),
       ylim = c(0.3,3))

  lines(pred_modif_rev_new,           ## cumulative exposure
        "overall",
        col = 4,
        ci = "area",
        ci.arg = list(density = 20, col = 4 ,angle = 45),
        lwd = 2)

  legend("topleft", legend = c("foehn days", "non foehn d."), col = c(2,4), lwd = 2)
}



#####


#### summer analysis (months 6:8), nothing significant -----
# summer data
data_summer = data[data$month %in% c(6:8),]
# grouping index
data_summer$station_year = paste0(data_summer$station, data_summer$year)
# crossbasis temp
cb.temp <- crossbasis(data_summer$temp,
                      lag=21,
                      argvar=list(fun="ns", knots = quantile(data_summer$temp, c(.5,.9), na.rm=TRUE)),
                      arglag=list(fun="ns", knots = logknots(21,3)),
                      group = data_summer$station_year)
par(mfrow=c(2,2))
for( i in c(10,35,80,140)){ #c(0,5,60,120
  foehn_bin <- ifelse(data_summer$f_id > i, 0, 1) # binary foehn index
  foehn_bin_rev <- ifelse(foehn_bin == 1, 0, 1)
  modif <- cb.temp * foehn_bin # modified cb.temp
  modif_rev <- cb.temp *foehn_bin_rev
  mod_modif <- gnm(all ~  cb.temp + modif, data = data_summer,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
  mod_modif_rev <- gnm(all ~  cb.temp + modif_rev, data = data_summer,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
  pred_modif  <- crosspred(cb.temp, mod_modif, cen = 20, cumul=FALSE)
  pred_modif2 <- crosspred(cb.temp, mod_modif_rev, cen = 20, cumul=FALSE)
  min1 <- findmin(cb.temp,pred_modif,from=quantile(data_summer$temp, .1),to=quantile(data_summer$temp, .9))
  min2 <- findmin(cb.temp,pred_modif2,from=quantile(data_summer$temp, .1),to=quantile(data_summer$temp, .9))

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


#### winter analysis (months 12:2), nothing significant -----
# winter data
data_winter = data[data$month %in% c(12,1,2),]
# grouping index
data_winter$station_year = paste(data_winter$station, ifelse(data_winter$month %in% 1:11, data_winter$year, data_winter$year+1))


# Filter out groups with fewer than 22 observations
# one group has only 13 observations
group_counts <- table(data_winter$station_year)
valid_groups <- names(group_counts[group_counts >= 21])
data_winter <- subset(data_winter, station_year %in% valid_groups)

# crossbasis temp
cb.temp <- crossbasis(data_winter$temp,
                      lag=21,
                      argvar=list(fun="ns", knots = quantile(data_winter$temp, c(.5,.9), na.rm=TRUE)),
                      arglag=list(fun="ns", knots = logknots(21,3)),
                      group = data_winter$station_year)

for( i in c(10,35,80,140)){ #c(0,5,60,120
  foehn_bin <- ifelse(data_winter$f_id > i, 0, 1) # binary foehn index
  foehn_bin_rev <- ifelse(foehn_bin == 1, 0, 1)
  modif <- cb.temp * foehn_bin # modified cb.temp
  modif_rev <- cb.temp *foehn_bin_rev
  mod_modif <- gnm(all ~  cb.temp + modif, data = data_winter,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
  mod_modif_rev <- gnm(all ~  cb.temp + modif_rev, data = data_winter,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
  pred_modif  <- crosspred(cb.temp, mod_modif, cen = 20, cumul=FALSE)
  pred_modif2 <- crosspred(cb.temp, mod_modif_rev, cen = 20, cumul=FALSE)
  min1 <- findmin(cb.temp,pred_modif,from=quantile(data_winter$temp, .1),to=quantile(data_winter$temp, .9))
  min2 <- findmin(cb.temp,pred_modif2,from=quantile(data_winter$temp, .1),to=quantile(data_winter$temp, .9))

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





# crossbasis temp
cb.temp <- crossbasis(data$temp,
                      lag=21,
                      argvar=list(fun="ns", knots = quantile(data$temp, c(.5,.9), na.rm=TRUE)),
                      arglag=list(fun="ns", knots = logknots(21,3)),
                      group = data$station)

# binary foehn
foehn_bin <- ifelse(data$f_id >= 72, 0, 1)
foehn_bin_rev <- ifelse(foehn_bin == 1, 0, 1)

# modifier functions
modif <- cb.temp * foehn_bin
modif_rev <- cb.temp * foehn_bin_rev

# groups
groups_id = colnames(data)[c(3,9,10,24,25, 13, 14, 11, 15, 12)]


par(mfrow=c(4,3))


  colvar  =  groups_id[i]

  formula  <- as.formula(paste0(colvar, "~ cb.foehn + modif"))
  formula2 <- as.formula(paste0(colvar, "~ cb.foehn + modif_rev"))

  # model with and without foehn
  mod_modif     <-gnm(formula, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
  mod_modif_rev <- gnm(formula2, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

  # prediction with and without foehn
  pred_modif   <- crosspred(cb.temp, mod_modif, cen = 20, cumul=FALSE)
  pred_modif2  <- crosspred(cb.temp, mod_modif_rev, cen = 20, cumul=FALSE)


  # get min value of both predictions and use for centering
  min1 <- findmin(cb.temp,pred_modif,from=quantile(data$temp, .1),to=quantile(data$temp, .9))
  min2 <- findmin(cb.temp,pred_modif2,from=quantile(data$temp, .1),to=quantile(data$temp, .9))

  # predict with new min values
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


  legend("topright", legend = c("temp + foehn", "temp"), col = c(2,4), lwd = 2)
