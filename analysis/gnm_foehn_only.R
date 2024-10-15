################################################################################
# GNM FOEHN
# - try out many different options to model foehn
# - determine the best

# notes:
# - foehn has a maximum LAG of 3 -> fun = "integer" advised
# - 4K BUFFER: 104166.81 (significantly less rows), var: lin, lag: strata (integer +4, worse than both ns)
# - 6K BUFFER: 134934.11 (all stations represented) var: lin, lag: strata (integer +3, but worse than ns with 1 knot), stratum with dow is way better than without, very little response
# - 8K BUFFER: 148163.11 var: list(fun = "strata", breaks = c(57.6, 115.2, 172.8, 230.4)) (lin +0.6, lin better for lag = integer), lag: strata (integer +4), a lot higher response
# - 10k BUFFER: 166022.31 var: lin, lag: strata (integer +3), stratum with dow is way better than without
# - 12K BUFFER: 180278.71 var: lin, lag: strata (integer +3.3)
# - 14k BUFFER: 184111.31 var: lin, lag: strata (integer +3), stratum with dow is way better than without
# conclusions:
# - smallest buffer radius possible -> 6k
# - 6k has most accurate response, 8k has a lot more response
# - argvar definetly "lin", arglag "stratum" (break at 1)
# - stratum with dow is far more accurate

# seasonal tests
# months 6:9
# - 8k, var: lin, lag: strata: overall slightly >1 RR
# - HOWEVER <1 RR for lag0, but >1 RR for lag1 etc
# months 11:3
# - 8k, var: Var function: list(fun = "ns", knots = c(`80%` = 4, `90%` = 72.29), Boundary = c(0, 288)), lag: strata
# - overall <1 RR between 80-200 foehn, clear lagged impact of heat during cold season
# - very complex pattern in winter

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
    mod <- gnm(all ~ cb.f_id,
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
                       argvar = eval(parse(text = opt_var)), # list(fun="lin") #
                       arglag = eval(parse(text = opt_lag)) # list(fun="integer") #
                       )
# model
mod_nm <- gnm(all ~ cb.foehn, data = data,  family=quasipoisson(), eliminate=stratum, subset=ind>0)
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
     lag  = 2,
     ci = "area",
     col = 5,
     ci.arg = list(density = 20, col = 5 ,angle = -45),
     xlab = "Exposure (Foehn)",
     ylab = "Relative Risk (RR)",
     main = "Exposure-Response at Lag of 2",
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
plot(pred_nm)

# par(mfrow=c(1,2))
# pacf(data$all, ylim=c(-0.1,1), main="Original response variable")
# pacf(residuals(mod_nm), ylim=c(-0.1,.1), na.action=na.pass,
#      main="Residuals from the regression model")

#####



dev.off()
