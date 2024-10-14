################################################################################
# GNM FOEHN
# - try out many different options to model foehn
# - determine the best

# - stratum is used to work with time trends and seasonality and station differences
# - foehn has a maximum LAG of 3 -> fun = "integer",
# - strata for lag doesnt make sense because temperature exposure doesnt behave like this (?)
# - lin for var ergibt sinn weil warum sollte der effect von foehn nichtlin sein

#TODO adjust the algorithm so that..
#...the functions are easier defined outside the loop
#...the results are saved in a matrix
#...this matrix has row and colnames (the funcitons)
#...the functions of the minimum qAIC are directly extracted and used for plotting

### PACKAGES ####
library(dlnm);library(splines);library(ggplot2);library(viridis);library(gnm)

######



### DATA ####
rm(list=ls())
data = read.csv("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MedStat_aggregated/centroid_aggregated/hosp_buffer_10000.csv")

# multiple time variables for stratum
data$date  = as.Date(data$date); data$dow   = as.factor(data$dow)
data$year  = as.factor(format(data$date, "%Y")); data$month = as.factor(format(data$date, "%m"))

# create stratum and index for deealing with trends/seasonality
data$stratum = with(data, factor(paste(station, year, month, dow, sep="-")))
# data$stratum = with(data, factor(paste(station, year, month, sep="-"))) -> inferior to above
ind = tapply(data$all, data$stratum, sum)

#####



### FUNCTIONS ####
# q-AIC computation
QAIC <- function(model) {
  phi <- summary(model)$dispersion
  loglik <- sum(dpois(model$y, model$fitted.values, log=TRUE))
  return(-2*loglik + 2*summary(model)$df[3]*phi)
}

######



### ALGORITHM ####

# distribution of foehn (for argvar)
quantile(data$f_id, seq(.75,1,0.01))


# two vectors of argvar and arglag arguments

v_var <-

v_lag <-





# CREATE TWO MATRIX WITH THE DIFFERENT COMB OF FUNCTIONS
comb <- expand.grid(1:length(v_var),1:length(v_lag)) # the number of specifications for argvar and arglag below

# CREATE AN EMPTY MATRIX BY F_ID MATRIX TO STORE THE QAIC
qaic_tab <- matrix(NA, nrow = 1, ncol=nrow(comb))



# RUN THE MODEL FOR EACH COMBINATION OF FUNCTIONS
for (j in 1:nrow(comb)){

  qaic <- rep(NA, 1) # what is the purpose of this?

  # define maximum lags !for foehn = 3!
  maxlago <- 3

  # different functions for argvar: exposure response association
  argvar1o <- list(fun="ns", knots = quantile(data$f_id, c(.8, .9), na.rm=TRUE),Boundary=range(data$f_id))
  # argvar2o <- list(fun="ns", knots = quantile(data$f_id, c(.8, .9, .95), na.rm=TRUE),Boundary=range(data$f_id))
  # argvar3o <- list(fun="ns", knots = quantile(data$f_id, c(.78, .86, .94), na.rm=TRUE),Boundary=range(data$f_id))
  # are both worse than number 1
  argvar2o <- list(fun="strata", breaks = 1)
  argvar3o <- list(fun="strata", breaks = 2)
  argvar4o <- list(fun="strata", breaks = equalknots(data$f_id, nk=3))
  argvar5o <- list(fun="ns", knots = equalknots(data$f_id, nk=2) ,Boundary=range(data$f_id))
  argvar6o <- list(fun="ns", knots = equalknots(data$f_id, nk=3) ,Boundary=range(data$f_id))
  # argvar6o <- list(fun="bs", degree=2,knots=quantile(data$f_id, c(.8, .9), na.rm=T))
  # argvar7o <- list(fun="bs", degree=2,knots=quantile(data$f_id, c(.8, .9, .95), na.rm=TRUE))
  # argvar6o <- list(fun="bs", degree=2,knots=quantile(data$f_id, c(.78, .86, .94), na.rm=TRUE))
  # perform both worse than the below, bs is inferior to lin
  argvar7o <- list(fun="lin")

  # function for arglag: exposure lag association
  # maxlago = 3 -> integer only needed
  arglag1o<-list(fun="integer")
  arglag2o<-list(fun="ns", knots = 1)
  arglag3o<-list(fun="ns", knots = c(1,2))
  arglag4o<-list(fun="strata", breaks = 1)


  # crossbasis
  cb.f_id <- crossbasis(data$f_id, lag=maxlago,
                        argvar=get(paste("argvar",comb[j,1],"o", sep="")),
                        arglag=get(paste("arglag",comb[j,2],"o", sep="")),
                        group=data$station)

  # MODEL
  mod <- gnm(all ~ cb.f_id,
             eliminate=stratum,
             family=quasipoisson(), data=data, na.action="na.exclude",
             subset=ind>0)

  qaic_tab[,j] <- QAIC(mod)

}

# Check model with lowest Q-AIC score
# which.min(qaic_tab) # qaic_tab


# create matrix of results
results <- matrix(qaic_tab, nrow = length(v_var), ncol = length(v_lag), byrow = FALSE)


# convert to data frame


# assign col and row names


# extract location of minimum value


# extract name of col and row and save them for plotting (the functions)
opt_var <-

opt_lag <-



# which(results == min(results), arr.ind = TRUE)
# row = argvar, col = arglag



#####




### VISUALIZATION of the best performing combination ####

# read different buffer data for sensitivity analysis
data = read.csv("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MedStat_aggregated/centroid_aggregated/hosp_buffer_10000.csv")
data$date  = as.Date(data$date); data$dow   = as.factor(data$dow); data$year  = as.factor(format(data$date, "%Y")); data$month = as.factor(format(data$date, "%m"))

data$stratum = with(data, factor(paste(station, year, month, dow, sep="-"))); ind = tapply(data$all, data$stratum, sum)

# crossbasis
cb.foehn <- crossbasis(data$f_id,lag=3,
                       argvar=opt_var, #
                       arglag=opt_lag
                       )

mod_nm <- gnm(all ~ cb.foehn, data = data,  family=quasipoisson(), eliminate=stratum, subset=ind>0)

pred_nm <- crosspred(cb.foehn, mod_nm, at=0:288, cumul=FALSE, cen = 0)


par(mfrow=c(2,2))
plot(pred_nm,              ## cumulative exposure
     "overall",
     col = 2,
     ci.arg = list(density = 20, col = 2 ,angle = -45),
     xlab = "Exposure",
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
