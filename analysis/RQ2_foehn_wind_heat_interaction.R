################################################################################
### RQ2: the foehn wind heat interaction
################################################################################



### Data
#----

rm(list=ls())

# libraries
library(dplyr); library(tidyr); library(plotly);library(zoo); library(viridis);
library(knitr);library(kableExtra);library(webshot2); library(magick)

# colors
colors <- viridis(3, option = "viridis")

# data
data = read.csv("../data/Medstat_hospitalizations_aggregated/hosp_buffer_8000.csv") |>
  mutate(date = as.Date(date),
         stratum_dow = as.factor(stratum_dow),
         y64 = a014y + a1564y,
         o64 = a6574y + a7584y + a85plusy,
         f_id_binary = ifelse(f_id >= 72, 1, 0),
         station_year = paste0(station, year))


# index to include only stratum that have hosp counts
ind_dow = tapply(data$all, data$stratum_dow, sum)

# MMT function
source("functions/findmin.R")

# define the maximum lag distance we account for
maxlago <- 3

# subpopulations
groups_id = c("all", "mal", "fem", "y64", "o64", "cvd", "resp", "inf", "uri", "ment" )

# letters for plotting
letters = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")

# colors
colors <- viridis(3, option = "viridis")
foehn_col = "green4";temp_col = "gold2"

#----


### Crossbasis definition
#----

# temperature
cb.temp <- crossbasis(data$temp,
                      lag=21,
                      argvar=list(fun="ns", knots = quantile(data$temp, c(.5,.9), na.rm=TRUE)),
                      arglag=list(fun="ns", knots = logknots(21,3)),
                      group = data$station)

# modifier functions
modif     <- cb.temp * data$f_id_bin
modif_rev <- cb.temp * data$f_id_bin_rev

#----



### Figure 3a
#----

# define subpopulation
colvar  =  "all"

# define formula for both model
formula1  <- as.formula(paste0(colvar, "~ cb.temp + modif"))
formula2 <- as.formula(paste0(colvar, "~ cb.temp + modif_rev"))

# model with and without foehn
mod_modif     <-gnm(formula1, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
mod_modif_rev <- gnm(formula2, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

# prediction to calculate MHT
pred_modif   <- crosspred(cb.temp, mod_modif, cen = 20)

# get min value for centering
min1 <- findmin(cb.temp, pred_modif, from=quantile(data$temp, .1),to=quantile(data$temp, .9))
  print(min1)

# predict with the MHT
pred_modif_new     <- crosspred(cb.temp, mod_modif, cen = min1, by = .1)
pred_modif_rev_new <- crosspred(cb.temp, mod_modif_rev, cen = min1, by = .1)


# save the plot
png("figures/Figure3a.png", width = 1350, height =1200, res = 300)
par(mfrow=c(1,1),
    mar = c(3,3,.5,.5),
    mgp = c(1.8, .5, 0))


plot(pred_modif_new,
     "overall",
     ylab = "relative risk",
     xlab = "temperature [\u00B0C]",
     col = foehn_col,
     ci.arg = list(col = alpha(colour = foehn_col, .15)),
     lwd = 2,
     main ="",
     ylim = c(0.7,2.5),
     cex.axis = 0.7,
     cex.lab = 0.7)


lines(pred_modif_rev_new,
      "overall",
      col = temp_col,
      ci = "area",
      ci.arg = list(col = alpha(colour = temp_col, .15)),
      lwd = 2)

abline(v =quantile(data$temp, .99), col = "black", lty = 2)

legend("top", ncol = 1, legend = c("temperature on foehn days (Model 4)", "temperature on non-foehn days (Model 4)"), col = c("green4", "gold2"),
       bty = "n", lwd=c(2,2), cex = 0.7)


dev.off()

#----





























