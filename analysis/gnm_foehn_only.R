################################################################################
# GNM FOEHN
# - try out many different options to model foehn
# - determine the best



### PACKAGES ####
library(dlnm);library(splines);library(ggplot2);library(viridis);library(gnm)

######



### DATA ####
rm(list=ls())
data = read.csv("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MedStat_aggregated/centroid_aggregated/hosp_buffer_10000.csv")

# multiple time variables for stratum
data$date  = as.Date(data$date)
data$year  = as.factor(format(data$date, "%Y"))
data$month = as.factor(format(data$date, "%m"))
data$dow   = as.factor(data$dow)

# create stratum and index
data$stratum = with(data, factor(paste(station, year, month, dow, sep="-")))
# data$stratum = with(data, factor(paste(station, year, month, sep="-")))
ind = tapply(data$all, data$stratum, sum)

#####



### ALGORITHM ####

# CREATE TWO MATRIX WITH THE DIFFERENT COMB OF FUNCTIONS
comb <- expand.grid(1:7,1:3) # the number of specifications for argvar and arglag below

# CREATE AN EMPTY MATRIX BY TEMPERATURE MATRIX TO STORE THE QAIC
qaic_tab <- matrix(NA, nrow = 1, ncol=nrow(comb))

# FUNCTION FOR COMPUTING THE Q-AIC
QAIC <- function(model) {
  phi <- summary(model)$dispersion
  loglik <- sum(dpois(model$y, model$fitted.values, log=TRUE))
  return(-2*loglik + 2*summary(model)$df[3]*phi)
}

# RUN THE MODEL FOR EACH COMBINATION OF FUNCTIONS
for (j in 1:nrow(comb)){

  qaic <- rep(NA, 1) # what is the purpose of this?

  # EXTRACT THE dat_sum_rt <- ?what is that? just the maximum lags?
  maxlago <- 3

  # SET DIFFERENT VAR FUNCTIONS FOR TMEAN
  argvar1o <- list(fun="ns", knots = quantile(data$f_id, c(0.50, 0.90), na.rm=TRUE),Boundary=range(data$f_id))
  argvar2o <- list(fun="ns", knots=quantile(data$f_id, c(0.75, 0.90), na.rm=TRUE),Boundary=range(data$f_id))
  argvar3o <- list(fun="ns", knots=quantile(data$f_id, c(0.50, 0.75, 0.90), na.rm=TRUE),Boundary=range(data$f_id))
  argvar4o <- list(fun="bs", degree=2,knots=quantile(data$f_id, c(0.50, 0.90), na.rm=T))
  argvar5o <- list(fun="bs", degree=2,knots=quantile(data$f_id, c(0.75, 0.90), na.rm=TRUE))
  argvar6o <- list(fun="bs", degree=2,knots=quantile(data$f_id, c(0.50, 0.75, 0.90), na.rm=TRUE))
  argvar7o <- list(fun="lin")

  # DIFFERENT FUNCTIONS OF THE LAG DIMENSION
  arglag1o<-list(fun="ns", knots = 1)
  arglag2o<-list(fun="strata", breaks=1)
  arglag3o<-list(fun="integer")


  #DEFINE CROSSBASIS
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

which.min(qaic_tab) # Check model with lowest Q-AIC score


#####
