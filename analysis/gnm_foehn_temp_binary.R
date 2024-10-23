################################################################################
# GNM FOEHN
# - try out many different options to model foehn when combined with temperature
# - determine the best binary code for foehn
# - because we are in binary, there are limited options for the argvar and arglag
# - argvar = integer, this is fixed due to our foehn only being 0 and 1
# - arglag = integer, because we only have three lags

# notes:


### PACKAGES ####
library(dlnm);library(splines);library(ggplot2);library(viridis);library(gnm)

######

rm(list=ls())

# buffer size for data file read
buffer = 8000

### DATA ####
data = read.csv(paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MedStat_aggregated/centroid_aggregated/hosp_buffer_", buffer, ".csv"))

data$date = as.Date(data$date)

# index to include only stratum that have hosp counts
data$stratum_dow = as.factor(data$stratum_dow); data$stratum = as.factor(data$stratum)
ind_dow = tapply(data$all, data$stratum_dow, sum); ind = tapply(data$all, data$stratum, sum)

data$station <- as.factor(data$station)
#####

### CROSSBASIS TEMPERATURE ####

cb.temp <- crossbasis(data$temp,
                      lag=21,
                      argvar=list(fun="ns", knots = quantile(data$temp, c(.1,.75,.9), na.rm=TRUE)),
                      arglag=list(fun="ns", knots = logknots(21,3)),
                      group = data$station)

#####

### FUNCTION qAIC ####
# q-AIC computation
QAIC <- function(model) {
  phi <- summary(model)$dispersion
  loglik <- sum(dpois(model$y, model$fitted.values, log=TRUE))
  return(-2*loglik + 2*summary(model)$df[3]*phi)
}

######

### ARGVAR ARGLAG DEFINITION ####

#  define var function
argvar = list(fun="integer")

#  define lag function
arglag = list(fun="integer")

#####


bin_thresholds = c(20,40,80,120,160, 180 ,200,240,270)

### ALGORITHM ####

# define the maximum lag distance we account for
maxlago <- 3

# create an empty matrix to store the qAIC
qaic_tab <- matrix(NA,
                   nrow = length(bin_thresholds),
                   ncol=1,
                   dimnames = list(c(bin_thresholds)))

## Run the model for each combination
for (i in 1:length(bin_thresholds)){

  # define binary data based on bin_threshold
  foehn_bin <- ifelse(data$f_id > i, 1, 0)

  # crossbasis
  cb.f_id <- crossbasis(foehn_bin,
                        lag=maxlago,
                        argvar=argvar,
                        arglag=arglag,
                        group=data$station)

  # model
  mod <- gnm(all ~ cb.f_id + cb.temp,
             data=data,
             eliminate=stratum,
             subset=ind_dow>0,
             family=quasipoisson())

  # save qAIC in qaic_tab
  qaic_tab[i,] <- QAIC(mod)

}


# Check model with lowest Q-AIC score
min_qaic = min(qaic_tab, na.rm = TRUE)

# extract location of minimum value
min_position <- which(qaic_tab == min_qaic, arr.ind = TRUE)

# extract name of col and row and save them for plotting (the functions)
opt_var <- rownames(qaic_tab)[min_position[1]]

# print results
print(paste0("Minimum value: ", round(min_qaic, 1), digits = 1))
cat("Binary Threshold:", opt_var, "\n")

View(qaic_tab)
#####

### VISUALIZATION of the best performing combination ####

# binary foehn
foehn_bin <- ifelse(data$f_id > 180, 1, 0) #  as.numeric(opt_var) c(20,40,80,120,160,200,240,270)


# crossbasis
cb.foehn <- crossbasis(foehn_bin,
                       lag = maxlago,
                       argvar = argvar,
                       arglag = arglag,
                       group = data$station)

# model with and without foehn
mod_nm <- gnm(all ~  cb.temp, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
mod_nm_f <- gnm(all ~ cb.foehn + cb.temp, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

# prediction with and without foehn
pred_nm <- crosspred(cb.temp, mod_nm, cen = 28, cumul=FALSE)
pred_nm_f <- crosspred(cb.temp, mod_nm_f, cen = 28, cumul=FALSE)




plot(pred_nm,              ## cumulative exposure
     "overall",
     col = 2,
     ci.arg = list(density = 10, col = 2 ,angle = -45),
     lwd = 5,
     main = "Overall cumulative exposure-response")

lines(pred_nm_f,           ## cumulative exposure
     "overall",
     col = 4,
     ci = "area",
     ci.arg = list(density = 10, col = 4 ,angle = 45),
     lwd = 2)

legend("topright", legend = c("temp", "temp+binary foehn"), col = c(2,4), lwd = 2)

#
#
# plot(pred_nm_f,              ## exposure at specific lag
#      "slices",
#      lag  = 0,
#      ci = "area",
#      col = 3,
#      ci.arg = list(density = 20, col = 3 ,angle = -45),
#      xlab = "Exposure (Foehn)",
#      ylab = "Relative Risk (RR)",
#      main = "Exposure-Response at Lag of 0",
#      lwd = 2
# )
#
# plot(pred_nm,              ## exposure at specific lag
#      "slices",
#      lag  = 1,
#      ci = "area",
#      col = 4,
#      ci.arg = list(density = 20, col = 4 ,angle = -45),
#      xlab = "Exposure (Foehn)",
#      ylab = "Relative Risk (RR)",
#      main = "Exposure-Response at Lag of 1",
#      lwd = 2
# )
#
#
# plot(pred_nm,              ## exposure at foehn increase
#      "slices",
#      lag  = 3,
#      ci = "area",
#      col = 5,
#      ci.arg = list(density = 20, col = 5 ,angle = -45),
#      xlab = "Exposure (Foehn)",
#      ylab = "Relative Risk (RR)",
#      main = "Exposure-Response at Lag of 3",
#      lwd = 2
# )
#
# for(i in c(50, 120, 200, 280)){
# plot(pred_nm,              ## exposure at foehn increase
#      "slices",
#      var  = i,
#      ci = "area",
#      col = 5,
#      ci.arg = list(density = 20, col = 5 ,angle = -45),
#      xlab = "Exposure (Foehn)",
#      ylab = "Relative Risk (RR)",
#      main = paste0("Exposure-Response at Exposure of ", i),
#      lwd = 2
# )
# }

par(mfrow=c(1,1))
plot(pred_nm,            ## 3D Plot
     "3d",
     xlab = "Exposure (Foehn)",
     ylab = "Lag (Days)",
     zlab = "Relative Risk",
     theta = 210,                    # The azimuthal angle (horizontal rotation)
     phi = 30,                      # The colatitude (vertical rotation).
     # col = "skyblue2",              # color of the 3d surface
     border = "black", # color of the borders
     #ticktype = "simple",           # type of grid on surface, alt.: "detailed"
     #nticks = 10                    # number of ticks
     r = 0.02
)

#####

### OTHER MODELS ####

# binary foehn
foehn_bin <- ifelse(data$f_id > 80, 1, 0) #  as.numeric(opt_var) c(20,40,80,120,160,200,240,270)


# crossbasis
cb.foehn <- crossbasis(foehn_bin,
                       lag = maxlago,
                       argvar = argvar,
                       arglag = arglag,
                       group = data$station)

# model with and without foehn
mod_nm <- glm(all ~  cb.temp + ns(X, 7*20), data = data,  family=quasipoisson())

mod_nm_f <- glm(all ~ cb.foehn + cb.temp + ns(X, 7*20), data = data,  family=quasipoisson())

# prediction with and without foehn
pred_nm <- crosspred(cb.temp, mod_nm, cen = 18, cumul=FALSE)

pred_nm_f <- crosspred(cb.temp, mod_nm_f, cen = 18, cumul=FALSE)




plot(pred_nm,              ## cumulative exposure
     "overall",
     col = 2,
     ci.arg = list(density = 20, col = 2 ,angle = -45),
     lwd = 2,
     main = "Overall cumulative exposure-response")

lines(pred_nm_f,           ## cumulative exposure
      "overall",
      col = 4,
      ci = "area",
      ci.arg = list(density = 20, col = 4 ,angle = 45),
      lwd = 2)

#####

dev.off()

