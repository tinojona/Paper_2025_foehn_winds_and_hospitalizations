################################################################################
### RQ2: the foehn wind heat interaction
################################################################################



### Data
#----

rm(list=ls())

# libraries
library(dlnm);library(splines);library(ggplot2);library(viridis);library(gnm);
library(dplyr);library(knitr);library(lmtest); library(grid);library(gridBase)
library(ggstance)

# colors
colors <- viridis(3, option = "viridis")

# data
data = read.csv("/Volumes/FS/_ISPM/CCH/Tino/master_thesis/data/Medstat_hospitalizations_aggregated/hosp_buffer_8000.csv") |>
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

# Cumulative relative risk from temperature exposure (Model 3) for all-cause
# hospitalizations with 95% confidence interval.

png("output/figures/Figure3a.png", width = 1800, height = 1800, res = 300)
par(
  mar = c(3.5, 3.5, 1,1),
  mgp = c(2, 0.7, 0))


mod <- gnm(all ~ cb.temp, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
pred <- crosspred(cb.temp, mod, cumul=FALSE, cen = 0)
min1 <- findmin(cb.temp,pred,from=quantile(data$temp, .1),to=quantile(data$temp, .9))
prednew <- crosspred(cb.temp, mod, cumul=FALSE, cen = min1)

plot(prednew,
     "overall",
     col = colors[1],
     ci = "area",
     ci.arg = list(col = alpha(colour = colors[1], 0.25)),
     xlab = "temperature [\u00B0C]",
     ylab = "relative risk",
     lwd = 2,
     main = "",
     cex.axis = 1,
     cex.lab = 1,
     ylim = c(0.8,2.3)
)


legend("top", ncol = 1,
      legend = c("temperature \n(Model 3)"),
      col = c(colors[1]), bty = "n", lwd=c(3), cex = 1.5)

dev.off()


#----

### Figure 3b
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
png("output/figures/Figure3b.png", width = 1350, height =1200, res = 300)
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
abline(v =quantile(data$temp, .01), col = "black", lty = 2)

legend("top", ncol = 1, legend = c("temperature on foehn days (Model 4)", "temperature on non-foehn days (Model 4)"), col = c("green4", "gold2"),
       bty = "n", lwd=c(2,2), cex = 0.7)


dev.off()

#----


### Table for Figure 3c
#----

# Figure 3c as a Table

# empty data frame
table_estimates = data.frame(all = rep(NA, 2),
                             mal = rep(NA, 2),
                             fem = rep(NA, 2),
                             y64 = rep(NA, 2),
                             o64  = rep(NA, 2),
                             cvd = rep(NA, 2),
                             resp = rep(NA, 2),
                             inf = rep(NA, 2),
                             uri = rep(NA, 2),
                             ment = rep(NA, 2),
                             row.names = c("Model_4_foehn_days", "Model_4_non_foehn_days"))

# loop through all subpopulations
for (i in 1:length(groups_id)) {

  # select subpopulation
  colvar  =  groups_id[i]

  # formula
  formula1  <- as.formula(paste0(colvar, "~ cb.temp + modif"))
  formula2 <- as.formula(paste0(colvar, "~ cb.temp + modif_rev"))

  # model with and without foehn
  mod_modif     <-gnm(formula1, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
  mod_modif_rev <- gnm(formula2, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

  # prediction with and without foehn
  pred_modif   <- crosspred(cb.temp, mod_modif, cen = 20, cumul=FALSE)

  # get min value of both predictions and use for centering
  min1 <- findmin(cb.temp,pred_modif,from=quantile(data$temp, .1),to=quantile(data$temp, .9))

  # predict with new min values
  pred_modif_new     <- crosspred(cb.temp, mod_modif, cen = min1, cumul=FALSE, by = .1)
  pred_modif_rev_new <- crosspred(cb.temp, mod_modif_rev, cen = min1, cumul=FALSE, by = .1)

  # extract prediction for value 24.7 and save it
  table_estimates[1,i] = paste0( sprintf("%.3f",round(pred_modif_new$allRRfit["24.7"],digits=3)),
                                 " [",
                                 sprintf("%.3f",round(pred_modif_new$allRRlow["24.7"], digits=3)),
                                 "-",
                                 sprintf("%.3f",round(pred_modif_new$allRRhigh["24.7"], digits = 3)), "]")

  table_estimates[2,i] = paste0( sprintf("%.3f",round(pred_modif_rev_new$allRRfit["24.7"],digits=3)),
                                 " [",
                                 sprintf("%.3f",round(pred_modif_rev_new$allRRlow["24.7"], digits=3)),
                                 "-",
                                 sprintf("%.3f",round(pred_modif_rev_new$allRRhigh["24.7"], digits = 3)), "]")



}

table_estimates = t(table_estimates)

#----



### Figure 3c
#----

# Cumulative relative risk (Model 4) for subpopulations at 24.7 °C with 95% confidence
# intervals.

# empty dataframe
table_estimates = data.frame(categories = c("all","all", "mal", "mal","fem","fem","y64","y64","o64","o64",
                                            "cvd","cvd","resp","resp","inf","inf","uri","uri","ment","ment"),
                             model = rep(c("temp + foehn", "temp - foehn"), 10),
                             pred = rep(NA, 20),
                             CI_low = rep(NA, 20),
                             CI_high = rep(NA, 20))


# loop through all specifications
for (i in 1:nrow(table_estimates)) {

  # extract current subpopulation
  colvar  =  table_estimates$categories[i]

  # uneven numbers is foehn scenario
  if( i %% 2 != 0) {
    formula1 = as.formula(paste0(colvar, "~ cb.temp + modif"))

    # even number is non-foehn scenario
  } else {
    formula1 = as.formula(paste0(colvar, "~ cb.temp + modif_rev"))
  }

  # GET UNIFORM MMT
  mod_modif     <-gnm(as.formula(paste0(colvar, "~ cb.temp + modif")), data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
  pred_modif   <- crosspred(cb.temp, mod_modif, cen = 20, cumul=FALSE)
  min1 <- findmin(cb.temp,pred_modif,from=quantile(data$temp, .1),to=quantile(data$temp, .9))

  # model with and without foehn
  mod_modif     <- gnm(formula1, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

  # predict with new min values
  pred_modif_new  <- crosspred(cb.temp, mod_modif, cen = min1, cumul=FALSE, by = .1)

  # extract relative risk and confidence intervals
  table_estimates$pred[i]    = pred_modif_new$allRRfit["24.7"]
  table_estimates$CI_low[i]  = pred_modif_new$allRRlow["24.7"]
  table_estimates$CI_high[i] = pred_modif_new$allRRhigh["24.7"]
}


table_estimates <- table_estimates |>
  mutate(categories = c("all", "all", "male", "male", "female", "female","<65 years", "<65 years", ">64 years", ">64 years", "circulatory", "circulatory", "respiratory", "respiratory", "infectious", "infectious", "genitourinary", "genitourinary" ,"mental", "mental"),
         categories = factor(categories, levels = rev(unique(categories))),
         model = factor(model, levels = c("temp + foehn", "temp - foehn")))




# save figure
png("figures/Figure3c.png", width = 2000, height =2500, res = 300)
par(
  mar = c(3.5, 6, 1,1),
  mgp = c(2.1, 0.7, 0))

# create empty plot
plot(1, type = "n", xlim = c(0.5, 2.5), ylim = c(0.5, length(table_estimates$categories) + 0.5),
     xlab = "relative risk", ylab = "", yaxt = "n", xaxt="n", bty = "n",
     cex.axis = 1, cex.lab = 1)

# create dashed 1-line
abline(v = 1, lty = 2, col = "black")

# define vertical offsets for different models
offset <- - 0.25
vertical_shift <- 0.3

# add error bars
for(i in 1:length(table_estimates$pred)) {

  # assign color based on foehn presence
  col <- ifelse(table_estimates$model[i] == "temp + foehn", foehn_col, temp_col)

  # adjust the vertical position by combining the shift and offset
  vertical_position <- length(table_estimates$categories) - i + 1 - vertical_shift

  # apply model-specific offset
  if (table_estimates$model[i] == "temp + foehn") {
    vertical_position <- vertical_position + offset
  } else {
    vertical_position <- vertical_position - offset
  }

  # draw the error bars
  segments(table_estimates$CI_low[i], vertical_position, table_estimates$CI_high[i], vertical_position, col = col, lwd = 2)
  # draw the points
  points(table_estimates$pred[i], vertical_position, pch = 16, col = col, cex = 1)
}

# add continuous y-axis with ticks every second label
axis(2, at = seq(1, length(table_estimates$categories), by = 2),
     labels = rev(table_estimates$categories)[seq(1, length(table_estimates$categories), by = 2)],
     las = 1,cex.axis = 1)

# add x-axis
axis(1, at = seq(0.5,2.5,0.5), cex.axis = 1)

# add borders at the x and y axis
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[3], border = "black", lwd = 1)  # Bottom border
rect(par("usr")[1], par("usr")[3], par("usr")[1], par("usr")[4], border = "black", lwd = 2)


dev.off()

#----


### Figure 3
#----

# run the code of Figure 3c!

# (a) Cumulative relative risk from temperature exposure (Model 3) for all-cause
# hospitalizations with 95% confidence interval. (b) Cumulative relative risk from
# temperature exposure (Model 4) for all-cause hospitalizations divided into foehn and
# non-foehn days with 95% confidence intervals. The dotted line indicates the temperature
# corresponding to the 99th percentile of the temperature distribution (24.7 °C). (c)
# Cumulative relative risk (Model 4) for subpopulations at 24.7 °C with 95% confidence
# intervals.


# save the plot
# png("output/figures/Figure3.png", width = 1650, height =1900, res = 300)
# pdf("output/figures/Figure3.pdf", width = 5.7, height =7)

layout_matrix <- matrix(c(1, 1, 2, 2,
                          1, 1, 4, 4,
                          3, 3, 4, 4,
                          3, 3, 4, 4),
                        nrow = 4, byrow = TRUE)

# Set up the layout
graphics::layout(mat = layout_matrix)

par(
  mar = c(3.5, 3.5, 1,1),
  mgp = c(2, 0.7, 0))

# PLOT 1 TOP LEFT : TEMPERATURE ONLY
cb.temp <- crossbasis(data$temp,
                      lag=21,
                      argvar=list(fun="ns", knots = quantile(data$temp, c(.5,.9), na.rm=TRUE)),
                      arglag=list(fun="ns", knots = logknots(21,3)),
                      group = data$station)
mod <- gnm(all ~ cb.temp, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
pred <- crosspred(cb.temp, mod, cumul=FALSE, cen = 0)
min1 <- findmin(cb.temp,pred,from=quantile(data$temp, .1),to=quantile(data$temp, .9))
prednew <- crosspred(cb.temp, mod, cumul=FALSE, cen = min1)

plot(prednew,
     "overall",
     col = colors[1],
     ci = "area",
     ci.arg = list(col = alpha(colour = colors[1], 0.25)),
     xlab = "temperature [\u00B0C]",
     ylab = "relative risk",
     lwd = 2,
     main = "",
     cex.axis = 1,
     cex.lab = 1,
     ylim = c(0.8,2.3)
)

text(24, 2 , labels = "(a)", pos = 4, cex = 1)



# PLOT 2 : LEGEND
plot.new()

par(
  mar = c(3.5, 6, 1,1),
  mgp = c(2, 0.7, 0))

legend("bottom", ncol = 1,
       legend = c("temperature \n(Model 3)",
                  "temperature on foehn days  \n(Model 4)",
                  "temperature on non-foehn days  \n(Model 4)"),
       col = c(colors[1], "green4", "gold2"), bty = "n", lwd=c(3,3,3), cex = 1, ,y.intersp = 1.7)


# PLOT 3 BOTTOM LEFT : TEMPERATURE ON FOEHN AND NON-FOEHN DAYS
par(
  mar = c(3.5, 3.5, 1,1),
  mgp = c(2, 0.7, 0))

colvar  =  groups_id[1]

formula1  <- as.formula(paste0(colvar, "~ cb.temp + modif"))
formula2 <- as.formula(paste0(colvar, "~ cb.temp + modif_rev"))

# model with and without foehn
mod_modif     <-gnm(formula1, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

# prediction with and without foehn
pred_modif   <- crosspred(cb.temp, mod_modif, cen = 20)

min1 <- findmin(cb.temp,pred_modif,from=quantile(data$temp, .1),to=quantile(data$temp, .9))

pred_modif_new     <- crosspred(cb.temp, mod_modif, cen = min1, by = .1)
pred_modif_rev_new <- crosspred(cb.temp, mod_modif_rev, cen = min1, by = .1)

plot(pred_modif_new,              ## cumulative exposure
     "overall",
     ylab = "relative risk",
     xlab = "temperature [\u00B0C]",
     col = "white",
     ci.arg = list(col = alpha(colour = "white", .15)),
     lwd = 2,
     main ="",
     ylim = c(0.8,2.3),
     cex.axis = 1,
     cex.lab = 1)

abline(v =quantile(data$temp, .99), col = "black", lty = 2)

lines(pred_modif_new,
      "overall",
      col = foehn_col,
      ci = "area",
      ci.arg = list(col = alpha(colour = foehn_col, .15)),
      lwd = 2)


lines(pred_modif_rev_new,
      "overall",
      col = temp_col,
      ci = "area",
      ci.arg = list(col = alpha(colour = temp_col, .15)),
      lwd = 2)

text(25, 2.15 , labels = "(b)", pos = 4, cex = 1)


# PLOT 4 BOTTOM RIGHT : Lollipop FIGURE
par(
  mar = c(3.5, 4.5, 1,1),
  mgp = c(2.1, 0.7, 0))

# create empty plot
plot(1, type = "n", xlim = c(0.5, 2.5), ylim = c(0.5, length(table_estimates$categories) + 0.5),
     xlab = "relative risk", ylab = "", yaxt = "n", xaxt="n", bty = "n",
     cex.axis = 1, cex.lab = 1)

# create dashed 1-line
abline(v = 1, lty = 2, col = "black")

# define vertical offsets for different models
offset <- - 0.25
vertical_shift <- 0.3

# add error bars
for(i in 1:length(table_estimates$pred)) {

  # assign color based on foehn presence
  col <- ifelse(table_estimates$model[i] == "temp + foehn", foehn_col, temp_col)

  # adjust the vertical position by combining the shift and offset
  vertical_position <- length(table_estimates$categories) - i + 1 - vertical_shift

  # apply model-specific offset
  if (table_estimates$model[i] == "temp + foehn") {
    vertical_position <- vertical_position + offset
  } else {
    vertical_position <- vertical_position - offset
  }

  # draw the error bars
  segments(table_estimates$CI_low[i], vertical_position, table_estimates$CI_high[i], vertical_position, col = col, lwd = 2)
  # draw the points
  points(table_estimates$pred[i], vertical_position, pch = 16, col = col, cex = 1)
}

# add continuous y-axis with ticks every second label
axis(2, at = seq(1, length(table_estimates$categories), by = 2),
     labels = rev(table_estimates$categories)[seq(1, length(table_estimates$categories), by = 2)],
     las = 1,cex.axis = 1)

# add x-axis
axis(1, at = seq(0.5,2.5,0.5), cex.axis = 1)

# add borders at the x and y axis
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[3], border = "black", lwd = 1)  # Bottom border
rect(par("usr")[1], par("usr")[3], par("usr")[1], par("usr")[4], border = "black", lwd = 2)



dev.off()

#----



### sFigure 5
#----

# Cumulative relative risk of the interaction between foehn winds and
# temperature with 95% confidence interval for all-cause hospitalizations.

# all
colvar  =  "all"

# formula
formula1 <- as.formula(paste0(colvar, "~ cb.temp + modif_rev"))

# model with foehn wind
mod_modif <- gnm(formula1, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

# prediction of the interaction term
pred_modif <- crosspred(modif_rev, mod_modif, cen = 15.5, by = 0.1)



# save the figure
png("output/figures/sFigure5.png", width = 1000, height =1000, res = 300)
par(mfrow=c(1,1),
    mar = c(3,3,.5,.5),
    mgp = c(1.8, .5, 0))

plot(pred_modif,              ## cumulative exposure
     "overall",
     ylab = "relative risk",
     xlab = "temperature [\u00B0C]",
     col = foehn_col,
     ci.arg = list(col = alpha(colour = foehn_col, .15)),
     # col = "grey2",
     # ci.arg = list(col = alpha(colour = "grey2", .15)),
     lwd = 2,
     main ="",
     ylim = c(0.7,2.5),
     # ylim = c(0.9,1.8),
     cex.axis = 0.7,
     cex.lab = 0.7)

abline(v=24.7, lty = 2)
abline(v=-8.9, lty = 2)

est_CI <- paste0(
  round(pred_modif$allRRfit["24.7"], digits = 3),
  " (CI: ",
  round(pred_modif$allRRlow["24.7"], digits = 3),
  "; " ,
  round(pred_modif$allRRhigh["24.7"], digits = 3),
  ")")

text(26, 2, labels = est_CI, pos = 2, cex = 0.9)

est_CI <- paste0(
  round(pred_modif$allRRfit["-8.9"], digits = 3),
  " (CI: ",
  round(pred_modif$allRRlow["-8.9"], digits = 3),
  "; " ,
  round(pred_modif$allRRhigh["-8.9"], digits = 3),
  ")")

text(26, 1.5, labels = est_CI, pos = 2, cex = 0.9)

dev.off()


#----



### sFigure 6
#----

# Cumulative relative risk with 95% confidence intervals for (a) male, (b)
# female, (c) 64 years and younger, (d) older than 64 years, (e) circulatory, (f) respiratory,
# (g) infectious, (h) genitourinary, (i) mental hospitalization with a binary foehn wind
# intensity threshold value of 72 which corresponds to 6 hours of full foehn wind. The
# green line shows the temperature hospitalization association when foehn winds were
# present, the yellow line when foehn winds were absent.

# save figure
png("output/figures/sFigure6.png", width = 1800, height =1800, res = 300)
par(mfrow=c(3,3),
    mar = c(3,3,.5,.5),
    mgp = c(1.8, .5, 0))

# loop though all subpopulations
for (i in 1:length(groups_id[2:10])) {

  colvar  =  groups_id[i+1]

  formula1  <- as.formula(paste0(colvar, "~ cb.temp + modif"))
  formula2 <- as.formula(paste0(colvar, "~ cb.temp + modif_rev"))

  # model with and without foehn
  mod_modif     <-gnm(formula1, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
  mod_modif_rev <- gnm(formula2, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

  # prediction with and without foehn
  pred_modif   <- crosspred(cb.temp, mod_modif, cen = 20, cumul=FALSE)


  # get min value of both predictions and use for centering
  min1 <- findmin(cb.temp,pred_modif,from=quantile(data$temp, .1),to=quantile(data$temp, .9))

  # predict with new min values
  pred_modif_new     <- crosspred(cb.temp, mod_modif, cen = min1, cumul=FALSE, by = .1)
  pred_modif_rev_new <- crosspred(cb.temp, mod_modif_rev, cen = min1, cumul=FALSE, by = .1)

  plot(pred_modif_new,              ## cumulative exposure
       "overall",
       ylab = "relative risk",
       xlab = "temperature [\u00B0C]",
       col = foehn_col,
       ci.arg = list(col = alpha(colour = foehn_col, .2)),
       lwd = 2,
       main ="",
       ylim = c(0.7,3)
  )


  lines(pred_modif_rev_new,           ## cumulative exposure
        "overall",
        col = temp_col,
        ci = "area",
        ci.arg = list(col = alpha(colour = temp_col, .2)),
        lwd = 2)

  text(25, 3, labels = paste0("(", letters[i], ")"), pos = 2)
}

dev.off()

#----




### sFigure 7
#----

# Cumulative relative risk of the interaction between foehn wind and temperature
# with 95% confidence interval for (a) male, (b) female, (c) 64 years and younger,
# (d) older than 64, (e) circulatory, (f) respiratory, (g) infectious, (h) genitourinary,
# mental hospitalizations. The relative risk and its confidence interval shown inside each
# correspond to the exposure of 24.7°C.

png("output/figures/sFigure7.png", width = 1800, height =1800, res = 300)
par(mfrow=c(3,3),
    mar = c(3,3,.5,.5),
    mgp = c(1.8, .5, 0))

for (i in 1:length(groups_id[2:10])) {

  # select subpopulation
  colvar  =  groups_id[i+1]

  # GET UNIFORM MMT
  mod_modif     <-gnm(as.formula(paste0(colvar, "~ cb.temp + modif")), data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
  pred_modif   <- crosspred(cb.temp, mod_modif, cen = 20, cumul=FALSE)
  min1 <- findmin(cb.temp,pred_modif,from=quantile(data$temp, .1),to=quantile(data$temp, .9))

  # formula for foehn and non-foehn days
  formula1  <- as.formula(paste0(colvar, "~ cb.temp + modif_rev"))

  # model with foehn
  mod_modif     <- gnm(formula1, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

  # prediction with foehn
  pred_modif   <- crosspred(modif_rev, mod_modif, cen = min1, by = 0.1)

  plot(pred_modif,              ## cumulative exposure
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

  text(-15, 2.4, labels = paste0("(", letters[i], ")"), pos = 2)

  abline(v=24.7, lty = 2)
  abline(v=-8.9, lty = 2)

  est_CI <- paste0(
    round(pred_modif$allRRfit["24.7"], digits = 3),
    " \n(CI: ",
    round(pred_modif$allRRlow["24.7"], digits = 3),
    "; " ,
    round(pred_modif$allRRhigh["24.7"], digits = 3),
    ")")

  text(26, 2.3, labels = est_CI, pos = 2, cex = 0.9)

  est_CI <- paste0(
    round(pred_modif$allRRfit["-8.9"], digits = 3),
    " \n(CI: ",
    round(pred_modif$allRRlow["-8.9"], digits = 3),
    "; " ,
    round(pred_modif$allRRhigh["-8.9"], digits = 3),
    ")")

  text(26, 1.9, labels = est_CI, pos = 2, cex = 0.9)



}

dev.off()

#----



### Model 3: 3D
#----

# not in the paper

mod <- gnm(all ~ cb.temp, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
prednew <- crosspred(cb.temp, mod, cumul=FALSE, cen = 20, by = .2)

png("output/figures/X_Model3_3D.png", width = 1800, height = 1800, res = 300)
par(
  mar = c(2, 3, 1,2),
  mgp = c(2.8, 2.5, 0))

plot(prednew,
     "3d",

     main = "",
     xlab = "temperature [°C]",
     zlab = "relative risk",

     zlim = c(0.95,1.1),
     col = colors[1],
     border = "white",
     # shade = 0.1,
     lwd = .2,
     ticktype = "detailed",

     expand = .8,
     theta = 60, # rotation
     phi = 25 # elevation

)

dev.off()

#----



### sTable 7
#----

# Results (relative risks) for the sensitivity analysis on the influence of the
# observation period (1998-2019 and 2008-2019) and the foehn wind aggregation method
# (only full foehn wind) on the interaction between temperature and foehn winds
# (Model 4) with fixed minimum hospitalization temperatures for every subpopulation.
# The left most 2 columns represent the baseline risk shown in Figure 3c.

# empty data frame for sensitivity analysis
table_estimates_review = data.frame(all = rep(NA, 6),
                                    mal = rep(NA, 6),
                                    fem = rep(NA, 6),
                                    y64 = rep(NA, 6),
                                    o64  = rep(NA, 6),
                                    cvd = rep(NA, 6),
                                    resp = rep(NA, 6),
                                    inf = rep(NA, 6),
                                    uri = rep(NA, 6),
                                    ment = rep(NA, 6),
                                    row.names = c("all: temp + foehn", "all: temp - foehn", "2008-: temp + foehn", "2008-: temp - foehn", "full foehn (foehn", "full foehn (non-foehn)"))

# save MMT
MMT_subgroups = rep(NA, ncol(table_estimates_review))

# the first loop for all data
for (i in 1:length(groups_id)) {

  colvar  =  groups_id[i]

  formula1  <- as.formula(paste0(colvar, "~ cb.temp + modif"))
  formula2 <- as.formula(paste0(colvar, "~ cb.temp + modif_rev"))

  # model with and without foehn
  mod_modif     <-gnm(formula1, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
  mod_modif_rev <- gnm(formula2, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

  # prediction with and without foehn
  pred_modif   <- crosspred(cb.temp, mod_modif, cen = 20, cumul=FALSE)

  # get min value of both predictions and use for centering
  min1 <- findmin(cb.temp,pred_modif,from=quantile(data$temp, .1),to=quantile(data$temp, .9))

  # save MMT
  MMT_subgroups[i] <- min1

  print(min1)

  # predict with new min values
  pred_modif_new     <- crosspred(cb.temp, mod_modif, cen = min1, cumul=FALSE, by = .1)
  pred_modif_rev_new <- crosspred(cb.temp, mod_modif_rev, cen = min1, cumul=FALSE, by = .1)

  # extract prediction for value 24.7 and save it
  table_estimates_review[1,i] = paste0( sprintf("%.3f",round(pred_modif_new$allRRfit["24.7"],digits=3)),
                                        " [",
                                        sprintf("%.3f",round(pred_modif_new$allRRlow["24.7"], digits=3)),
                                        "-",
                                        sprintf("%.3f",round(pred_modif_new$allRRhigh["24.7"], digits = 3)), "]")

  table_estimates_review[2,i] = paste0( sprintf("%.3f",round(pred_modif_rev_new$allRRfit["24.7"],digits=3)),
                                        " [",
                                        sprintf("%.3f",round(pred_modif_rev_new$allRRlow["24.7"], digits=3)),
                                        "-",
                                        sprintf("%.3f",round(pred_modif_rev_new$allRRhigh["24.7"], digits = 3)), "]")

}



# redefine the data for 2008
data_2008 <- data |> filter(data$date >= "2008-01-01")

# index to include only stratum that have hosp counts
ind_dow_2008 = tapply(data_2008$all, data_2008$stratum_dow, sum); ind_2008 = tapply(data_2008$all, data_2008$stratum, sum)

# crossbasis temp
cb.temp_2008 <- crossbasis(data_2008$temp,
                           lag=21,
                           argvar=list(fun="ns", knots = quantile(data_2008$temp, c(.5,.9), na.rm=TRUE)),
                           arglag=list(fun="ns", knots = logknots(21,3)),
                           group = data_2008$station)



# modifier functions
modif_2008 <- cb.temp_2008 * data_2008$f_id_bin
modif_rev_2008 <- cb.temp_2008 * data_2008$f_id_bin_rev


# the second loop for post 2008 data
for (i in 1:length(groups_id)) {

  colvar  =  groups_id[i]

  formula1  <- as.formula(paste0(colvar, "~ cb.temp_2008 + modif_2008"))
  formula2 <- as.formula(paste0(colvar, "~ cb.temp_2008 + modif_rev_2008"))

  # model with and without foehn
  mod_modif     <-gnm(formula1, data = data_2008,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow_2008>0)
  mod_modif_rev <- gnm(formula2, data = data_2008,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow_2008>0)

  # predict with new min values
  pred_modif_new     <- crosspred(cb.temp, mod_modif, cen = MMT_subgroups[i], cumul=FALSE, by = .1)
  pred_modif_rev_new <- crosspred(cb.temp, mod_modif_rev, cen = MMT_subgroups[i], cumul=FALSE, by = .1)

  # extract prediction for value 24.7 and save it
  table_estimates_review[3,i] = paste0( sprintf("%.3f",round(pred_modif_new$allRRfit["24.7"],digits=3)),
                                        " [",
                                        sprintf("%.3f",round(pred_modif_new$allRRlow["24.7"], digits=3)),
                                        "-",
                                        sprintf("%.3f",round(pred_modif_new$allRRhigh["24.7"], digits = 3)), "]")

  table_estimates_review[4,i] = paste0( sprintf("%.3f",round(pred_modif_rev_new$allRRfit["24.7"],digits=3)),
                                        " [",
                                        sprintf("%.3f",round(pred_modif_rev_new$allRRlow["24.7"], digits=3)),
                                        "-",
                                        sprintf("%.3f",round(pred_modif_rev_new$allRRhigh["24.7"], digits = 3)), "]")

}





# third sensitivity for foehn wind aggregation
data_ff = read.csv("/Volumes/FS/_ISPM/CCH/Tino/master_thesis/data/Medstat_hospitalizations_aggregated/hosp_buffer_8000sensitivity_onlyfullfoehnaggregation.csv") |>
  mutate(station_date = paste0(station, date),
         date = as.Date(date),
         station = as.factor(station),
         stratum_dow = as.factor(stratum_dow),
         stratum = as.factor(stratum),
         y64 = a014y + a1564y,
         o64 = a6574y + a7584y + a85plusy,
         foehn_bin = ifelse(f_id_sens > 72, 0, 1), # binary foehn foehn_bin <- ifelse(data$f_id > 72, 0, 1);
         foehn_bin_rev = ifelse(foehn_bin == 1, 0, 1)) # foehn_bin_rev <- ifelse(foehn_bin == 1, 0, 1)

# index to include only stratum that have hosp counts
ind_dow_ff = tapply(data$all, data$stratum_dow, sum); ind_ff = tapply(data$all, data$stratum, sum)

# crossbasis temp
cb.temp_ff <- crossbasis(data_ff$temp,
                         lag=21,
                         argvar=list(fun="ns", knots = quantile(data_ff$temp, c(.5,.9), na.rm=TRUE)),
                         arglag=list(fun="ns", knots = logknots(21,3)),
                         group = data_ff$station)

# modifier functions
modif_ff <- cb.temp_ff * data_ff$foehn_bin
modif_rev_ff <- cb.temp_ff * data_ff$foehn_bin_rev

# loop through the subpops while using the MMT from the OG data
for (i in 1:length(groups_id)) {

  colvar  =  groups_id[i]

  # model the rr in the full foehn aggregated sets
  formula1  <- as.formula(paste0(colvar, "~ cb.temp_ff + modif_ff"))
  formula2  <- as.formula(paste0(colvar, "~ cb.temp_ff + modif_rev_ff"))

  # model with and without foehn
  mod_modif     <-gnm(formula1, data = data_ff,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow_ff>0)
  mod_modif_rev     <-gnm(formula2, data = data_ff,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow_ff>0)

  # prediction with and without foehn
  pred_modif_new   <- crosspred(cb.temp_ff, mod_modif, cen = MMT_subgroups[i], cumul=FALSE, by = 0.1)
  pred_modif_rev_new   <- crosspred(cb.temp_ff, mod_modif_rev, cen = MMT_subgroups[i], cumul=FALSE, by = 0.1)

  # extract prediction for value 24.7 and save it
  table_estimates_review[5,i] = paste0( sprintf("%.3f",round(pred_modif_new$allRRfit["24.7"],digits=3)),
                                         " [",
                                         sprintf("%.3f",round(pred_modif_new$allRRlow["24.7"], digits=3)),
                                         "-",
                                         sprintf("%.3f",round(pred_modif_new$allRRhigh["24.7"], digits = 3)), "]")

  table_estimates_review[6,i] = paste0( sprintf("%.3f",round(pred_modif_rev_new$allRRfit["24.7"],digits=3)),
                                         " [",
                                         sprintf("%.3f",round(pred_modif_rev_new$allRRlow["24.7"], digits=3)),
                                         "-",
                                         sprintf("%.3f",round(pred_modif_rev_new$allRRhigh["24.7"], digits = 3)), "]")

}



table_estimates_review = t(table_estimates_review)


# save table
write.csv(table_estimates_review, "output/tables/sTable7.csv")


#----


### Table for Figure 3c !cold!
#----

# Figure 3c as a Table !cold!

# empty data frame
table_estimates = data.frame(all = rep(NA, 2),
                             mal = rep(NA, 2),
                             fem = rep(NA, 2),
                             y64 = rep(NA, 2),
                             o64  = rep(NA, 2),
                             cvd = rep(NA, 2),
                             resp = rep(NA, 2),
                             inf = rep(NA, 2),
                             uri = rep(NA, 2),
                             ment = rep(NA, 2),
                             row.names = c("Model_4_foehn_days", "Model_4_non_foehn_days"))

# loop through all subpopulations
for (i in 1:length(groups_id)) {

  # select subpopulation
  colvar  =  groups_id[i]

  # formula
  formula1  <- as.formula(paste0(colvar, "~ cb.temp + modif"))
  formula2 <- as.formula(paste0(colvar, "~ cb.temp + modif_rev"))

  # model with and without foehn
  mod_modif     <-gnm(formula1, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
  mod_modif_rev <- gnm(formula2, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

  # prediction with and without foehn
  pred_modif   <- crosspred(cb.temp, mod_modif, cen = 20, cumul=FALSE)

  # get min value of both predictions and use for centering
  min1 <- findmin(cb.temp,pred_modif,from=quantile(data$temp, .1),to=quantile(data$temp, .9))

  # predict with new min values
  pred_modif_new     <- crosspred(cb.temp, mod_modif, cen = min1, cumul=FALSE, by = .1)
  pred_modif_rev_new <- crosspred(cb.temp, mod_modif_rev, cen = min1, cumul=FALSE, by = .1)

  # extract prediction for value 24.7 and save it
  table_estimates[1,i] = paste0( sprintf("%.3f",round(pred_modif_new$allRRfit["-8.9"],digits=3)),
                                 " [",
                                 sprintf("%.3f",round(pred_modif_new$allRRlow["-8.9"], digits=3)),
                                 "-",
                                 sprintf("%.3f",round(pred_modif_new$allRRhigh["-8.9"], digits = 3)), "]")

  table_estimates[2,i] = paste0( sprintf("%.3f",round(pred_modif_rev_new$allRRfit["-8.9"],digits=3)),
                                 " [",
                                 sprintf("%.3f",round(pred_modif_rev_new$allRRlow["-8.9"], digits=3)),
                                 "-",
                                 sprintf("%.3f",round(pred_modif_rev_new$allRRhigh["-8.9"], digits = 3)), "]")



}

table_estimates = t(table_estimates)

#----



### Figure 3c !cold!
#----

# Cumulative relative risk (Model 4) for subpopulations at 24.7 °C with 95% confidence
# intervals.

# empty dataframe
table_estimates = data.frame(categories = c("all","all", "mal", "mal","fem","fem","y64","y64","o64","o64",
                                            "cvd","cvd","resp","resp","inf","inf","uri","uri","ment","ment"),
                             model = rep(c("temp + foehn", "temp - foehn"), 10),
                             pred = rep(NA, 20),
                             CI_low = rep(NA, 20),
                             CI_high = rep(NA, 20))


# loop through all specifications
for (i in 1:nrow(table_estimates)) {

  # extract current subpopulation
  colvar  =  table_estimates$categories[i]

  # uneven numbers is foehn scenario
  if( i %% 2 != 0) {
    formula1 = as.formula(paste0(colvar, "~ cb.temp + modif"))

    # even number is non-foehn scenario
  } else {
    formula1 = as.formula(paste0(colvar, "~ cb.temp + modif_rev"))
  }

  # GET UNIFORM MMT
  mod_modif     <-gnm(as.formula(paste0(colvar, "~ cb.temp + modif")), data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
  pred_modif   <- crosspred(cb.temp, mod_modif, cen = 20, cumul=FALSE)
  min1 <- findmin(cb.temp,pred_modif,from=quantile(data$temp, .1),to=quantile(data$temp, .9))

  # model with and without foehn
  mod_modif     <- gnm(formula1, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

  # predict with new min values
  pred_modif_new  <- crosspred(cb.temp, mod_modif, cen = min1, cumul=FALSE, by = .1)

  # extract relative risk and confidence intervals
  table_estimates$pred[i]    = pred_modif_new$allRRfit["-8.9"]
  table_estimates$CI_low[i]  = pred_modif_new$allRRlow["-8.9"]
  table_estimates$CI_high[i] = pred_modif_new$allRRhigh["-8.9"]
}


table_estimates <- table_estimates |>
  mutate(categories = c("all", "all", "male", "male", "female", "female","<65 years", "<65 years", ">64 years", ">64 years", "circulatory", "circulatory", "respiratory", "respiratory", "infectious", "infectious", "genitourinary", "genitourinary" ,"mental", "mental"),
         categories = factor(categories, levels = rev(unique(categories))),
         model = factor(model, levels = c("temp + foehn", "temp - foehn")))




# save figure
png("output/figures/Figure3c_cold.png", width = 2000, height =2500, res = 300)
par(
  mar = c(3.5, 6, 1,1),
  mgp = c(2.1, 0.7, 0))

# create empty plot
plot(1, type = "n", xlim = c(0.5, 2.5), ylim = c(0.5, length(table_estimates$categories) + 0.5),
     xlab = "relative risk", ylab = "", yaxt = "n", xaxt="n", bty = "n",
     cex.axis = 1, cex.lab = 1)

# create dashed 1-line
abline(v = 1, lty = 2, col = "black")

# define vertical offsets for different models
offset <- - 0.25
vertical_shift <- 0.3

# add error bars
for(i in 1:length(table_estimates$pred)) {

  # assign color based on foehn presence
  col <- ifelse(table_estimates$model[i] == "temp + foehn", foehn_col, temp_col)

  # adjust the vertical position by combining the shift and offset
  vertical_position <- length(table_estimates$categories) - i + 1 - vertical_shift

  # apply model-specific offset
  if (table_estimates$model[i] == "temp + foehn") {
    vertical_position <- vertical_position + offset
  } else {
    vertical_position <- vertical_position - offset
  }

  # draw the error bars
  segments(table_estimates$CI_low[i], vertical_position, table_estimates$CI_high[i], vertical_position, col = col, lwd = 2)
  # draw the points
  points(table_estimates$pred[i], vertical_position, pch = 16, col = col, cex = 1)
}

# add continuous y-axis with ticks every second label
axis(2, at = seq(1, length(table_estimates$categories), by = 2),
     labels = rev(table_estimates$categories)[seq(1, length(table_estimates$categories), by = 2)],
     las = 1,cex.axis = 1)

# add x-axis
axis(1, at = seq(0.5,2.5,0.5), cex.axis = 1)

# add borders at the x and y axis
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[3], border = "black", lwd = 1)  # Bottom border
rect(par("usr")[1], par("usr")[3], par("usr")[1], par("usr")[4], border = "black", lwd = 2)


dev.off()

#----



### Figure 3 - only all cause temperature and modification
#----

# run the code of Figure 3c!

# (a) Cumulative relative risk from temperature exposure (Model 3) for all-cause
# hospitalizations with 95% confidence interval. (b) Cumulative relative risk from
# temperature exposure (Model 4) for all-cause hospitalizations divided into foehn and
# non-foehn days with 95% confidence intervals. The dotted line indicates the temperature
# corresponding to the 1st and 99th percentile of the temperature distribution (-8.9 and 24.7 °C).


# save the plot
png("output/figures/Figure3_only_allcause_bluered.png", width = 2400, height =1850, res = 300)
# pdf("output/figures/Figure3.pdf", width = 5.7, height =7)

# define layout
layout_matrix <- matrix(c(1, 1, 1, 1,
                          2, 2, 3, 3),
                        nrow = 2, byrow = TRUE)

graphics::layout(mat = layout_matrix, heights = c(0.2, 0.8))


# PLOT 1 TOP: LEGEND
plot.new()

par(
  mar = c(0, 0, 0,0),
  mgp = c(0, 0, 0))

legend("center", ncol = 3,
       legend = c("temperature (Model 3)",
                  "temperature on foehn days (Model 4)",
                  "temperature on non-foehn days (Model 4)"),
       # col = c(colors[1], "green4", "gold2"),
       col = c(colors[1], "brown2", "steelblue"),
       bty = "n", lwd=c(2,2,2), cex = 1.3,
       y.intersp = 0,
       x.intersp = 0.8,
       text.width = c(0.215, 0.355, 0.43)
       )







# PLOT 2 BOTTOM LEFT: TEMPERATURE ONLY
par(
  mar = c(3.5, 3.5, 1,1),
  mgp = c(2, 0.7, 0))


cb.temp <- crossbasis(data$temp,
                      lag=21,
                      argvar=list(fun="ns", knots = quantile(data$temp, c(.5,.9), na.rm=TRUE)),
                      arglag=list(fun="ns", knots = logknots(21,3)),
                      group = data$station)
mod <- gnm(all ~ cb.temp, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
pred <- crosspred(cb.temp, mod, cumul=FALSE, cen = 0)
min1 <- findmin(cb.temp,pred,from=quantile(data$temp, .1),to=quantile(data$temp, .9))
prednew <- crosspred(cb.temp, mod, cumul=FALSE, cen = min1)

plot(prednew,
     "overall",
     col = colors[1],
     ci = "area",
     ci.arg = list(col = alpha(colour = colors[1], 0.25)),
     xlab = "temperature [\u00B0C]",
     ylab = "relative risk",
     lwd = 2,
     main = "",
     cex.axis = 1.3,
     cex.lab = 1.3,
     ylim = c(0.8,2.3)
)

text(24, 2.15 , labels = "(a)", pos = 4, cex = 1.3)


# PLOT 3 BOTTOM RIGHT : TEMPERATURE ON FOEHN AND NON-FOEHN DAYS
par(
  mar = c(3.5, 3.5, 1,1),
  mgp = c(2, 0.7, 0))

colvar  =  groups_id[1]

formula1  <- as.formula(paste0(colvar, "~ cb.temp + modif"))
formula2 <- as.formula(paste0(colvar, "~ cb.temp + modif_rev"))

# model with and without foehn
mod_modif     <-gnm(formula1, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
mod_modif_rev     <-gnm(formula2, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

# prediction with and without foehn
pred_modif   <- crosspred(cb.temp, mod_modif, cen = 20)

min1 <- findmin(cb.temp,pred_modif,from=quantile(data$temp, .1),to=quantile(data$temp, .9))

pred_modif_new     <- crosspred(cb.temp, mod_modif, cen = min1, by = .1)
pred_modif_rev_new <- crosspred(cb.temp, mod_modif_rev, cen = min1, by = .1)

plot(pred_modif_new,              ## cumulative exposure
     "overall",
     ylab = "relative risk",
     xlab = "temperature [\u00B0C]",
     col = "white",
     ci.arg = list(col = alpha(colour = "white", .15)),
     lwd = 2,
     main ="",
     ylim = c(0.8,2.3),
     cex.axis = 1.3,
     cex.lab = 1.3)

abline(v =quantile(data$temp, .99), col = "black", lty = 2)
abline(v =quantile(data$temp, .01), col = "black", lty = 2)

lines(pred_modif_new,
      "overall",
      # col = foehn_col,
      col = "brown2",
      ci = "area",
      # ci.arg = list(col = alpha(colour = foehn_col, .15)),
      ci.arg = list(col = alpha(colour = "brown2", .15)),
      lwd = 2)


lines(pred_modif_rev_new,
      "overall",
      # col = temp_col,
      col = "steelblue",
      ci = "area",
      # ci.arg = list(col = alpha(colour = temp_col, .15)),
      ci.arg = list(col = alpha(colour = "steelblue", .15)),
      lwd = 2)

text(25, 2.15 , labels = "(b)", pos = 4, cex = 1.3)


dev.off()

#----



### Figure 4 - temperature modification for cold and heat (two plots!)
#----

# Cumulative relative risk (Model 4) for subpopulations at 24.7 °C with 95% confidence
# intervals.

# empty dataframe
table_estimates = data.frame(categories = c("all","all", "mal", "mal","fem","fem","y64","y64","o64","o64",
                                            "cvd","cvd","resp","resp","inf","inf","uri","uri","ment","ment"),
                             model = rep(c("temp + foehn", "temp - foehn"), 10),
                             pred = rep(NA, 20),
                             CI_low = rep(NA, 20),
                             CI_high = rep(NA, 20))

table_estimates_cold = data.frame(categories = c("all","all", "mal", "mal","fem","fem","y64","y64","o64","o64",
                                            "cvd","cvd","resp","resp","inf","inf","uri","uri","ment","ment"),
                             model = rep(c("temp + foehn", "temp - foehn"), 10),
                             pred = rep(NA, 20),
                             CI_low = rep(NA, 20),
                             CI_high = rep(NA, 20))


# loop through all specifications
for (i in 1:nrow(table_estimates)) {

  # extract current subpopulation
  colvar  =  table_estimates$categories[i]

  # uneven numbers is foehn scenario
  if( i %% 2 != 0) {
    formula1 = as.formula(paste0(colvar, "~ cb.temp + modif"))

    # even number is non-foehn scenario
  } else {
    formula1 = as.formula(paste0(colvar, "~ cb.temp + modif_rev"))
  }

  # GET UNIFORM MMT
  mod_modif     <-gnm(as.formula(paste0(colvar, "~ cb.temp + modif")), data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)
  pred_modif   <- crosspred(cb.temp, mod_modif, cen = 20, cumul=FALSE)
  min1 <- findmin(cb.temp,pred_modif,from=quantile(data$temp, .1),to=quantile(data$temp, .9))

  # model with and without foehn
  mod_modif     <- gnm(formula1, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

  # predict with new min values
  pred_modif_new  <- crosspred(cb.temp, mod_modif, cen = min1, cumul=FALSE, by = .1)

  # extract relative risk and confidence intervals
  table_estimates$pred[i]    = pred_modif_new$allRRfit["24.7"]
  table_estimates$CI_low[i]  = pred_modif_new$allRRlow["24.7"]
  table_estimates$CI_high[i] = pred_modif_new$allRRhigh["24.7"]

  # extract relative risk and confidence intervals
  table_estimates_cold$pred[i]    = pred_modif_new$allRRfit["-8.9"]
  table_estimates_cold$CI_low[i]  = pred_modif_new$allRRlow["-8.9"]
  table_estimates_cold$CI_high[i] = pred_modif_new$allRRhigh["-8.9"]
}


table_estimates <- table_estimates |>
  mutate(categories = c("all", "all", "male", "male", "female", "female","<65 years", "<65 years", ">64 years", ">64 years", "circulatory", "circulatory", "respiratory", "respiratory", "infectious", "infectious", "genitourinary", "genitourinary" ,"mental", "mental"),
         categories = factor(categories, levels = rev(unique(categories))),
         model = factor(model, levels = c("temp + foehn", "temp - foehn")))

table_estimates_cold <- table_estimates_cold |>
  mutate(categories = c("all", "all", "male", "male", "female", "female","<65 years", "<65 years", ">64 years", ">64 years", "circulatory", "circulatory", "respiratory", "respiratory", "infectious", "infectious", "genitourinary", "genitourinary" ,"mental", "mental"),
         categories = factor(categories, levels = rev(unique(categories))),
         model = factor(model, levels = c("temp + foehn", "temp - foehn")))













# save figure
png("output/figures/Figure4.png", width = 1700, height =1550, res = 300)

layout_matrix <- matrix(c(1, 1, 1, 2, 2,
                          3, 3, 3, 4, 4),
                        nrow = 2, byrow = TRUE)

graphics::layout(mat = layout_matrix, heights = c(0.15, 0.85))

# Set up the layout
# graphics::layout(mat = layout_matrix)



# LEGENDS
par(mar = c(0, 0, 0,0))
plot.new()


legend(x = 0.2, y = 1 , ncol = 1,
       legend = c("cold on foehn days  (Model 4)",
                  "cold on non-foehn days  (Model 4)"),
       col = c("steelblue", "steelblue"), bty = "n", lwd=c(2,2), lty = c(3,1), cex = c(.9,.9),
       y.intersp = 1.5,
       x.intersp = 1,
       # text.width = c(0.45, 0.55),
       seg.len = 1.5
)

par(mar = c(0, 0, 0,0))
plot.new()



legend(x = -0.04, y = 1 , ncol = 1,
       legend = c("heat on foehn days  (Model 4)",
                  "heat on non-foehn days  (Model 4)"),
       col = c("brown2", "brown2"), bty = "n", lwd=c(2,2), lty = c(3,1), cex = c(.9,.9),
       y.intersp = 1.5,
       x.intersp = 1,
       seg.len = 1.5
       # text.width = c(0.45, 0.55)
)


#  COLD PLOT
par(
  mar = c(3.5, 6.5, 0,3),
  mgp = c(2.1, 0.7, 0))


# define vertical offsets for different models
offset <- - 0.2
vertical_shift <- 0.33

spacing <- 0.6
n_cat <- length(table_estimates_cold$categories)

plot(1, type = "n", xlim = c(0.5, 3),
     # ylim = c(0.5, length(table_estimates_cold$categories) * 1.5 + 0.5),
     ylim = c(0.5, n_cat * spacing - 0.8),
     xlab = "relative risk", ylab = "", yaxt = "n", xaxt="n", bty = "n",
     cex.axis = 1, cex.lab = 1)

# create dashed 1-line
abline(v = 1, lty = 2, col = "black")

# add error bars
for(i in 1:length(table_estimates_cold$pred)) {

  # assign color based on foehn presence
  lty <- ifelse(table_estimates_cold$model[i] == "temp + foehn", 3, 1)
  col <- "steelblue"

  vertical_position <- (n_cat - i + 1) * spacing - vertical_shift

  # and apply offset
  if (table_estimates_cold$model[i] == "temp + foehn") {
    vertical_position <- vertical_position + offset
  } else {
    vertical_position <- vertical_position - offset
  }

  # draw the error bars
  segments(table_estimates_cold$CI_low[i], vertical_position, table_estimates_cold$CI_high[i], vertical_position, col = col, lwd = 2, lty = lty)
  # draw the points
  points(table_estimates_cold$pred[i], vertical_position, pch = 16, col = col, cex = 1)
}

axis(2,
     at = seq(1, n_cat, by = 2) * spacing,
     labels = rev(table_estimates_cold$categories)[seq(1, n_cat, by = 2)],
     las = 1, cex.axis = 1)

# add x-axis
axis(1, at = seq(0.5,2.5,0.5), cex.axis = 1)

# add borders at the x and y axis
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[3], border = "black", lwd = 1)  # Bottom border
rect(par("usr")[1], par("usr")[3], par("usr")[1], par("usr")[4], border = "black", lwd = 2)

text(2.5, 10, labels = paste0("(a)"), pos = 2)


# HEAT PLOT

par(mar = c(3.5, 0, 0,0.3),
    mgp = c(2.1, 0.7, 0))


plot(1, type = "n", xlim = c(0.5, 3),
     # ylim = c(0.5, length(table_estimates$categories) * 1.5 + 0.5),
     ylim = c(0.5, n_cat * spacing - 0.8),
     xlab = "relative risk", ylab = "", yaxt = "n", xaxt="n", bty = "n",
     cex.axis = 1, cex.lab = 1)

# create dashed 1-line
abline(v = 1, lty = 2, col = "black")


# add error bars
for(i in 1:length(table_estimates$pred)) {

  # assign color based on foehn presence
  lty <- ifelse(table_estimates_cold$model[i] == "temp + foehn", 3, 1)
  col <- "brown2"

  vertical_position <- (n_cat - i + 1) * spacing - vertical_shift

  # and apply offset
  if (table_estimates_cold$model[i] == "temp + foehn") {
    vertical_position <- vertical_position + offset
  } else {
    vertical_position <- vertical_position - offset
  }

  # draw the error bars
  segments(table_estimates$CI_low[i], vertical_position, table_estimates$CI_high[i], vertical_position, col = col, lwd = 2, lty = lty)
  # draw the points
  points(table_estimates$pred[i], vertical_position, pch = 16, col = col, cex = 1)
}

axis(2,
     at = seq(1, n_cat, by = 2) * spacing,
     labels = FALSE,
     las = 1, cex.axis = 1)

# add x-axis
axis(1, at = seq(0.5,2.5,0.5), cex.axis = 1)

# add borders at the x and y axis
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[3], border = "black", lwd = 1)  # Bottom border
rect(par("usr")[1], par("usr")[3], par("usr")[1], par("usr")[4], border = "black", lwd = 2)

text(2.5, 10, labels = paste0("(b)"), pos = 2)

dev.off()

#----


