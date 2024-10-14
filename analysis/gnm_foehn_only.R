################################################################################
# GNM FOEHN
# - try out many different options to model foehn
# - determine the best

# - stratum is used to work with time trends and seasonality and station differences
# - foehn has a maximum LAG of 3 -> fun = "integer",
# - strata for lag works better than integer
# - lin for var ergibt sinn weil warum sollte der effect von foehn nicht lin sein


### PACKAGES ####
library(dlnm);library(splines);library(ggplot2);library(viridis);library(gnm)

######



### DATA ####
rm(list=ls())
data = read.csv("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MedStat_aggregated/centroid_aggregated/hosp_buffer_10000.csv")

data$date = as.Date(data$date)

# index to include only stratum that have hosp counts
data$stratum_dow = as.factor(data$stratum_dow); data$stratum = as.factor(data$stratum)
ind_dow = tapply(data$all, data$stratum_dow, sum); ind = tapply(data$all, data$stratum, sum)



# data$stratum = with(data, factor(paste(station, year, month, sep="-"))) -> inferior to above
ind = tapply(data$all, data$stratum_dow, sum)

#####



### FUNCTIONS ####
# q-AIC computation
QAIC <- function(model) {
  phi <- summary(model)$dispersion
  loglik <- sum(dpois(model$y, model$fitted.values, log=TRUE))
  return(-2*loglik + 2*summary(model)$df[3]*phi)
}

######


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


### ALGORITHM ####

# distribution of foehn (for argvar)
# quantile(data$f_id, seq(.75,1,0.01))
# mean(data$f_id[data$f_id != 0])
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
               eliminate=stratum_dow,
               subset=ind>0,
               family=quasipoisson())

    # save qAIC in qaic_tab
    qaic_tab[i,j] <- QAIC(mod)
  }
}


# Check model with lowest Q-AIC score
min_qaic = min(qaic_tab)

# extract location of minimum value
min_position <- which(qaic_tab == min_qaic, arr.ind = TRUE)

# extract name of col and row and save them for plotting (the functions)
opt_var <- rownames(qaic_tab)[min_position[1]]
opt_lag <- colnames(qaic_tab)[min_position[2]]

# print results
print(paste0("Minimum value:", round(min_qaic, 1), digits = 1))
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
                       argvar = eval(parse(text = opt_var)),
                       arglag = eval(parse(text = opt_lag)) # list(fun="integer")
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
     var  = 150,
     ci = "area",
     col = 5,
     ci.arg = list(density = 20, col = 5 ,angle = -45),
     xlab = "Exposure (Foehn)",
     ylab = "Relative Risk (RR)",
     main = "Exposure-Response at Exposure of 150",
     lwd = 2
)


#####
