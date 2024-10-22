################################################################################
# GNM FOEHN
# - try out many different options to model foehn when combined with temperature
# - determine the best binary code for foehn

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

### BINARY CODES ####
# TODO
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
# two lists of argvar and arglag arguments
v_var <- list(list(fun = "integer"),
              list(fun="strata")
)

v_lag <- list(list(fun="integer"),
              list(fun="strata", breaks = 1)
              )

#####


# TODO
# restructure the algorithm so that the matrix at the end has a third dimension:
# the different binary definitions

### ALGORITHM ####

# define the maximum lag distance we account for
maxlago <- 3

# create an empty matrix to store the qAIC
qaic_tab <- matrix(NA,
                   nrow = length(v_var),
                   ncol=length(v_lag),
                   dimnames = list(c(v_var), c(v_lag)))

## Run the model for each combination
for (i in 1:length(v_var)){

  # extract variable function
  argvar = v_var[[i]]

  for (j in 1:length(v_lag)) {

    #  extract lag function
    arglag = v_lag[[j]]

    # crossbasis
    cb.f_id <- crossbasis(data$f_id,
                          lag=maxlago,
                          argvar=argvar,
                          arglag=arglag,
                          group=data$station)

    # model
    mod <- gnm(all ~ cb.f_id + cb.temp,
               data=data,
               eliminate=stratum,
               subset=ind>0,
               family=quasipoisson())

    # save qAIC in qaic_tab
    qaic_tab[i,j] <- QAIC(mod)
  }
}


# Check model with lowest Q-AIC score
min_qaic = min(qaic_tab, na.rm = TRUE)

# extract location of minimum value
min_position <- which(qaic_tab == min_qaic, arr.ind = TRUE)

# extract name of col and row and save them for plotting (the functions)
opt_var <- rownames(qaic_tab)[min_position[1]]
opt_lag <- colnames(qaic_tab)[min_position[2]]

# print results
print(paste0("Minimum value: ", round(min_qaic, 1), digits = 1))
cat("Var function:", opt_var, "\n")
cat("Lag function:", opt_lag, "\n")

#####

### VISUALIZATION of the best performing combination ####

# crossbasis
cb.foehn <- crossbasis(data$f_id,lag = 3,
                       argvar = eval(parse(text = opt_var)), # list(fun = "thr",  thr.value = 40), #
                       arglag = eval(parse(text = opt_lag)), # list(fun="integer"), #
                       group = data$station)
# model
mod_nm <- gnm(all ~ cb.foehn + cb.temp, data = data,  family=quasipoisson(), eliminate=stratum, subset=ind>0)

# prediction
pred_nm <- crosspred(cb.foehn, mod_nm, at=0:288, cumul=FALSE, cen = 0)


par(mfrow=c(2,2))

plot(pred_nm,              ## cumulative exposure
     "overall",
     col = 2,
     ci.arg = list(density = 20, col = 2 ,angle = -45),
     xlab = "Exposure (Foehn)",
     ylab = "Cumulative Response",
     lwd = 2,
     main = "Overall cumulative exposure-response")

plot(pred_nm,              ## exposure at specific lag
     "slices",
     lag  = 0,
     ci = "area",
     col = 3,
     ci.arg = list(density = 20, col = 3 ,angle = -45),
     xlab = "Exposure (Foehn)",
     ylab = "Relative Risk (RR)",
     main = "Exposure-Response at Lag of 0",
     lwd = 2
)

plot(pred_nm,              ## exposure at specific lag
     "slices",
     lag  = 1,
     ci = "area",
     col = 4,
     ci.arg = list(density = 20, col = 4 ,angle = -45),
     xlab = "Exposure (Foehn)",
     ylab = "Relative Risk (RR)",
     main = "Exposure-Response at Lag of 1",
     lwd = 2
)


plot(pred_nm,              ## exposure at foehn increase
     "slices",
     lag  = 3,
     ci = "area",
     col = 5,
     ci.arg = list(density = 20, col = 5 ,angle = -45),
     xlab = "Exposure (Foehn)",
     ylab = "Relative Risk (RR)",
     main = "Exposure-Response at Lag of 3",
     lwd = 2
)

for(i in c(50, 120, 200, 280)){
plot(pred_nm,              ## exposure at foehn increase
     "slices",
     var  = i,
     ci = "area",
     col = 5,
     ci.arg = list(density = 20, col = 5 ,angle = -45),
     xlab = "Exposure (Foehn)",
     ylab = "Relative Risk (RR)",
     main = paste0("Exposure-Response at Exposure of ", i),
     lwd = 2
)
}

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



dev.off()

