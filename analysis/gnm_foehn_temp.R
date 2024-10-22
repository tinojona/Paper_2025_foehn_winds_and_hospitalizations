################################################################################
# GNM FOEHN
# - try out many different options to model foehn when combined with temperature
# - determine the best

# - foehn (0:288) is not clear enough, CI too large, to many variations depending on the model type
# - 8k buffer has the strongest association -> continue analysis with this one, other ones for sensitivity analysis
# - linear for variable modeling seems to perform best because too little data (20%) actually has foehn in it, all others are 0
#   thats why we change our approach to treshold and binary definition of the foehn impact







# notes:
# -
# - 4K BUFFER: 113763.41 (worse than without temp), var: lin, lag: strata (integer +1.5, but same performance as one ns)
# - 6K BUFFER: 147566.91 (""), vars: "", (integer +3, equal to ns), very small response
# - 8K BUFFER:161159.61 (decent compared to foehn only), vars: "", (Int +3, equal to ns), larger response
# - 10k BUFFER: 179451.21, vars, "" (Int +3, baad), similar response as 8k
# - 12K BUFFER: 194425.11, vars: "", (Int +2, ""), smaller response
# - 14k BUFFER: 198525.11 (worse than without temp)
# conclusions:
# - with accounting for temperature we have higher qAICs
# - we also loose all the <1 RR associations and only have >1 RR associations
# - the stronger the foehn the stronger the association
# - visually, at 8k buffer, we have the strongest association

# seasonal tests at 8k
# months 6:9
# - small >1 RR total, <1 RR at lag 0 decreasing with foehn score
# - >1 RR driven bei positive lags
# - when lag: int, we see a similar picture but with a more lagged effect
# months 11:3
# - var: ns equalknors (4), lag: strata.
# - small foehn-> >>1 RR, medium foehn-> <1 RR, large foehn-> >>1 RR,
# - with "lin" and "integer, classic cold response to temperature with high overall response

### PACKAGES ####
library(dlnm);library(splines);library(ggplot2);library(viridis);library(gnm)

######

rm(list=ls())

# buffer size for data file read
buffer = 8000

### DATA ####
data = read.csv(paste0("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MedStat_aggregated/centroid_aggregated/hosp_buffer_", buffer, ".csv"))

data$date = as.Date(data$date)

# SUBSET FOR SEASONAL ANALYSIS
# data <- subset(data, month %in% c(6:9))
# data <- subset(data, month %in% c(11,12,1,2,3))

# index to include only stratum that have hosp counts
data$stratum_dow = as.factor(data$stratum_dow); data$stratum = as.factor(data$stratum)
ind_dow = tapply(data$all, data$stratum_dow, sum); ind = tapply(data$all, data$stratum, sum)

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
v_var <- list(list(fun="ns", knots = quantile(data$f_id, c(.8, .9), na.rm=TRUE),Boundary=range(data$f_id)),
              list(fun="ns", knots = quantile(data$f_id, c(.8, .9, .95), na.rm=TRUE),Boundary=range(data$f_id)),
              list(fun="ns", knots = quantile(data$f_id, c(.9, .95, .99), na.rm=TRUE),Boundary=range(data$f_id)),
              list(fun="strata", breaks = equalknots(data$f_id, nk = 3)),
              list(fun="strata", breaks = equalknots(data$f_id, nk = 4)),
              list(fun="strata", breaks = equalknots(data$f_id, nk = 5)),
              list(fun="ns", knots = equalknots(data$f_id, nk=2) ,Boundary=range(data$f_id)),
              list(fun="ns", knots = equalknots(data$f_id, nk=3) ,Boundary=range(data$f_id)),
              list(fun="ns", knots = equalknots(data$f_id, nk=4) ,Boundary=range(data$f_id)),
              list(fun="lin")
)

v_lag <- list(list(fun="integer"),
              list(fun="strata", breaks = 1),
              list(fun="ns", knots = 1),
              list(fun="ns", knots = c(1,2))
)

#####

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

# # read different buffer data for sensitivity analysis
# data = read.csv("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MedStat_aggregated/centroid_aggregated/hosp_buffer_10000.csv")
# # index to include only stratum that have hosp counts
# data$stratum_dow = as.factor(data$stratum_dow); data$stratum = as.factor(data$stratum)
# ind_dow = tapply(data$all, data$stratum_dow, sum); ind = tapply(data$all, data$stratum, sum)

# crossbasis
cb.foehn <- crossbasis(data$f_id,lag = 3,
                       argvar = eval(parse(text = opt_var)), # list(fun="lin"), #
                       arglag = eval(parse(text = opt_lag)), # list(fun="integer") #
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

# par(mfrow=c(1,2))
# pacf(data$all, ylim=c(-0.1,1), main="Original response variable")
# pacf(residuals(mod_nm), ylim=c(-0.1,.1), na.action=na.pass,
#      main="Residuals from the regression model")

#####



dev.off()

