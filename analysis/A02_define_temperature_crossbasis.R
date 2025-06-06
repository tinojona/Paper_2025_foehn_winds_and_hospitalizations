################################################################################
### Temperature Crossbasis
################################################################################
# - visualize and justify the choice of the temperature crossbasis


### Data
#----

rm(list=ls())

# libraries
library(dlnm);library(splines);library(ggplot2);library(viridis);library(gnm);library(dplyr)

# colors
colors <- viridis(3, option = "viridis")

# data
data = read.csv("../data/Medstat_hospitalizations_aggregated/hosp_buffer_8000.csv") |>
  mutate(date = as.Date(date),
         stratum_dow = as.factor(stratum_dow),
         y64 = a014y + a1564y,
         o64 = a6574y + a7584y + a85plusy)


# index to include only stratum that have hosp counts
ind_dow = tapply(data$all, data$stratum_dow, sum)

# q-AIC computation
source("functions/qAIC_determination.R")

# MMT function
source("functions/findmin.R")

# define the maximum lag distance we account for
maxlago <- 3

#----


### CROSSBASIS TEMPERATURE
#----

cb.temp <- crossbasis(data$temp,
                      lag=21,
                      argvar=list(fun="ns", knots = quantile(data$temp, c(.5,.9), na.rm=TRUE)),
                      arglag=list(fun="ns", knots = logknots(21,3)),
                      group = data$station)

# cb.temp <- crossbasis(data$temp,
#                       lag=21,
#                       argvar=list(fun="ns", knots = quantile(data$temp, c(.1,.75,.9), na.rm=TRUE)),
#                       arglag=list(fun="ns", knots = logknots(21,3)),
#                       group = data$station)

#----


### Visualize the chosen model
#----

# model
mod <- gnm(all ~ cb.temp, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
# print(QAIC(mod)) 135721.5

# prediction preliminary
pred <- crosspred(cb.temp, mod, cumul=FALSE, cen = 0)

# find min value for plotting
min1 <- findmin(cb.temp,pred,from=quantile(data$temp, .1),to=quantile(data$temp, .9))

# prediction final
prednew <- crosspred(cb.temp, mod, cumul=FALSE, cen = min1)


png("figures/Model3.png", width = 1000, height =1000, res = 300)
par(mfrow=c(1,1),
    mar = c(3,3,.5,.5),
    mgp = c(1.8, .5, 0))


plot(prednew,              ## cumulative exposure
     "overall",
     col = colors[1],
     ci = "area",
     ci.arg = list(col = alpha(colour = colors[1], 0.25)), # list(col = alpha(colour = foehn_col, .15)),
     xlab = "temperature [\u00B0C]",
     ylab = "relative risk",
     lwd = 2,
     main = "",
     cex.axis = 1,
     cex.lab = 1)

dev.off()

#----




### Comparison of our crossbasis to Gasparrinis
#----

# original crossbasis
cb.tempog <- crossbasis(data$temp,
                      lag=21,
                      argvar=list(fun="ns", knots = quantile(data$temp, c(.1,.75,.9), na.rm=TRUE)),
                      arglag=list(fun="ns", knots = logknots(21,3)),
                      group = data$station)


# model
modog <- gnm(all ~ cb.tempog, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
# print(QAIC(modog)) 135724.9


# prediction final
prednewog <- crosspred(cb.tempog, modog, cumul=FALSE, cen = min1)


png("figures/Model3_comparison_with_original_crossbasis.png", width = 1000, height =1000, res = 300)
par(mfrow=c(1,1),
    mar = c(3,3,.5,.5),
    mgp = c(1.8, .5, 0))

plot(prednew,
     "overall",
     col = colors[1],
     ci = "area",
     ci.arg = list(col = alpha(colour = colors[1], 0.25)),
     xlab = "temperature [\u00B0C]",
     ylab = "relative risk",
     lwd = 2,
     main = "",
     cex.axis = 0.7,
     cex.lab = 0.7)

lines(prednewog,
     "overall",
     col = colors[3],
     ci = "area",
     ci.arg = list(col = alpha(colour = colors[3], 0.25)),
     lwd = 2,
     main = "")

legend("topright", legend = c("adjusted", "not adjusted"), col = c(colors[1], colors[3]), lwd = 1.5, bty = "n", cex = 0.8)

dev.off()

#----



