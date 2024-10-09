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
  
  qaic <- rep(NA, 1)
  
  # EXTRACT THE dat_sum_rt
  maxlago <- 3
  
  # SET DIFFERENT VAR FUNCTIONS FOR TMEAN
  argvar1o <- list(fun="ns",knots=quantile(data$tmean, c(0.50, 0.90), na.rm=TRUE),Boundary=range(data$tmean))
  argvar2o <- list(fun="ns",knots=quantile(data$tmean, c(0.75, 0.90), na.rm=TRUE),Boundary=range(data$tmean))
  argvar3o <- list(fun="ns",knots=quantile(data$tmean, c(0.50, 0.75, 0.90), na.rm=TRUE),Boundary=range(data$tmean))
  argvar4o <- list(fun="bs",degree=2,knots=quantile(data$tmean, c(0.50, 0.90), na.rm=T))
  argvar5o <- list(fun="bs",degree=2,knots=quantile(data$tmean, c(0.75, 0.90), na.rm=TRUE))
  argvar6o <- list(fun="bs",degree=2,knots=quantile(data$tmean, c(0.50, 0.75, 0.90), na.rm=TRUE))
  argvar7o <- list(fun="lin")
  
  # DIFFERENT FUNCTIONS OF THE LAG DIMENSION
  arglag1o<-list(fun="ns", knots = 1)
  arglag2o<-list(fun="strata", breaks=1)
  arglag3o<-list(fun="integer")
  

  # DEFINE THE STRATA
  data$stratum <- with(data, factor(paste(station, year, month, sep="-")))
  
  # ELIMINATE EMPTY STRATA (FOR CORRECT COMPUTATION OF CI IN gnm)
  ind <- tapply(data[,idcol], data$stratum, sum)[data$stratum]
  
  #DEFINE CROSSBASIS TEMPERATURE
  cbtmean <- crossbasis(data$tmean, lag=maxlago, 
                        argvar=get(paste("argvar",comb[j,1],"o", sep="")), 
                        arglag=get(paste("arglag",comb[j,2],"o", sep="")),
                        group=group)
  
  mod <- gnm(y ~ cbtmean, 
             eliminate=stratum,
             family=quasipoisson(), data=data, na.action="na.exclude",
             subset=ind>0)
  
  qaic_tab[,j] <- QAIC(mod)
  
}

which.min(qaic_tab) # Check model with lowest Q-AIC score