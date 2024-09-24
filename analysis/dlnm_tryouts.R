rm(list=ls())

# packages
library(dlnm)
library(splines)

# load data (for one station to begin with)
data = read.csv("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/MedStat_aggregated/hosp_buffer_5000.csv")

# HANDLE NAs HERE SET TO 0
data$f_id[is.na(data$f_id)] <- 0
data_short = data[1:5114,]



# crossbasis functions for foehn and temp
cb2.foehn <- crossbasis(data_short$f_id, lag=5, argvar=list(fun="lin"),
                        arglag=list(fun="poly",degree=4))
cb2.temp <- crossbasis(data_short$temp, lag=3, argvar=list(df=5),
                       arglag=list(fun="strata",breaks=1))

# for multiple stations the group function cuts so that there is no continuous lag from
# one station to the other; the same for temp
# cb2.foehn <- crossbasis(data$foehn, lag=5,
#                      argvar=list(fun="lin"), arglag=list(fun="integer"),
#                      group= data$station)


# model for impact on all hospitalisations
model2 <- glm(all ~ cb2.foehn + cb1.temp + ns(X, 7*14) + dow,
              family=quasipoisson(), data_short)

# predict the impact of foehn on hosp
pred1.pm <- crosspred(cb2.foehn, model2, at=0:288, bylag=0.2, cumul=FALSE)

par(mfrow=c(2,2))
# Different plots for visualization
plot(pred1.pm, "overall", main = "Overall cumulative exposure-response")         # Overall cumulative exposure-response
plot(pred1.pm, "slices", var = 30, main = "Lag-response for specific exposure") # Lag-response for specific exposure
plot(pred1.pm, "3d", main ="3D Plot")              # 3D plot
plot(pred1.pm, "contour", main = "contour plot")         # Contour plot

dev.off()










# plot the association, var defines the 1 day delta of foehn
par(mfrow=c(1,2))
plot(pred1.pm, "slices", var=50, col=3, ylab="RR", ci.arg=list(density=15,lwd=2),
     main="Association with a 50-unit increase in foehn")
plot(pred1.pm, "slices", var=50, col=2, cumul=FALSE, ylab="Cumulative RR",
     main="Cumulative association with a 50-unit increase in foehn")



# cumulative plot
plot(pred1.pm,"overall",xlab="foehn",ci="l", col=3, ylim=c(0.7,1.3), lwd=2,
     ci.arg=list(col=1,lty=3), main="Overall cumulative association for 5 lags")

# 3D plot
plot(pred1.pm, xlab="foehn", zlab="RR", theta=200, phi=40, lphi=30,
     main="3D graph of foehn effect")


par(mfrow=c(1,1))
plot(pred1.pm,"contour",xlab="foehn",main="Overall cumulative association for 5 lags")







