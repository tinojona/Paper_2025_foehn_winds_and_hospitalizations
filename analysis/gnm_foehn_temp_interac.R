################################################################################
# GNM FOEHN as a modifier of temperature
# - the lower the binary threshold (0 is best) the better the model performance

# notes:


### PACKAGES ####
library(dlnm);library(splines);library(ggplot2);library(viridis);library(gnm);library(mgcv)

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

### CROSSBASIS TEMPERATURE ######
# crossbasis temp
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


bin_thresholds = c(0,5)#,10,20,40,80,120,160, 180 ,200,240,270)

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
  # careful: the main term represents the risk when the specific interaction term is 0
  foehn_bin <- ifelse(data$f_id > i, 0, 1)

  # interaction term
  modif <- cb.temp * foehn_bin

  # model
  mod <- gnm(all ~ cb.temp + modif,
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
par(mfrow=c(2,2))
for( i in c(10,35,80,140)){ #c(0,5,60,120

# binary foehn
foehn_bin <- ifelse(data$f_id > i, 0, 1) #  as.numeric(opt_var) c(20,40,80,120,160,200,240,270)
foehn_bin_rev <- ifelse(foehn_bin == 1, 0, 1)


modif <- cb.temp * foehn_bin
modif_rev <- cb.temp *foehn_bin_rev

# model with and without foehn
mod_modif <- gnm(all ~  cb.temp + modif, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
mod_modif_rev <- gnm(all ~  cb.temp + modif_rev, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

# return significance of modif term in model
# coef_summary <- summary(mod_modif)$coefficients
# p_value_mod_modif <- coef_summary["mod_modif", "Pr(>|z|)"]
# p_value_mod_modif


# prediction with and without foehn
pred_modif <- crosspred(cb.temp, mod_modif, cen = 20, cumul=FALSE)
pred_modif_rev <- crosspred(cb.temp, mod_modif_rev, cen = 20, cumul=FALSE)



plot(pred_modif,              ## cumulative exposure
     "overall",
     col = 2,
     ci.arg = list(density = 10, col = 2 ,angle = -45),
     lwd = 2,
     main = paste0("Overall cum exp-resp: modifier, binary thr=",i))

lines(pred_modif_rev,           ## cumulative exposure
     "overall",
     col = 4,
     ci = "area",
     ci.arg = list(density = 10, col = 4 ,angle = 45),
     lwd = 2)

legend("topright", legend = c("temp+modif", "temp-modif"), col = c(2,4), lwd = 2)
}
#####

dev.off()

