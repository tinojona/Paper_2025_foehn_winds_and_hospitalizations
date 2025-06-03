################################################################################
### RQ1: the foehn wind hospitalization association
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
groups_id = colnames(data)[c(3,9,10,26,27, 13, 14, 11, 15, 12)]

#----


### Crossbasis definition
#----

# temperature
cb.temp <- crossbasis(data$temp,
                      lag=21,
                      argvar=list(fun="ns", knots = quantile(data$temp, c(.5,.9), na.rm=TRUE)),
                      arglag=list(fun="ns", knots = logknots(21,3)),
                      group = data$station)

# foehn wind
cb.foehn <- crossbasis(data$f_id,
                       lag = 3,
                       argvar = list(fun="lin"),
                       arglag = list(fun="integer"),
                       group = data$station)

#----


### Table 2
#----

# Cumulative relative risk for hospitalizations by cause for the effect of foehn
# winds (Model 1) and foehn winds adjusted for daily mean temperature (Model 2) across
# 3 days of lag for an increase of 6h full foehn wind equivalent. 95% confidence intervals
# are given in square brackets.

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
                             row.names = c("RR foehn [CI]", "RR foehn + temp [CI]"))

for (i in 1:length(groups_id)) {

  colvar = groups_id[i]

  # formula, only foehn
  formula1 <- as.formula(paste(colvar, "~ cb.foehn"))

  # model
  mod_nm1 <- gnm(formula1, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

  # prediction
  pred_nm1 <- crosspred(cb.foehn, mod_nm1, at=0:288, cumul=FALSE, cen = 0)

  # extract prediction for value 72 and save it
  table_estimates[1,i] = paste0( sprintf("%.3f",round(pred_nm1$allRRfit["72"],digits=3)),
                                 " [",
                                 sprintf("%.3f",round(pred_nm1$allRRlow["72"], digits=3)),
                                 "-",
                                 sprintf("%.3f",round(pred_nm1$allRRhigh["72"], digits = 3)), "]")



  # formula, foehn + temp
  formula2 <- as.formula(paste(colvar, "~ cb.foehn + cb.temp"))

  # model
  mod_nm2 <- gnm(formula2, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

  # prediction
  pred_nm2 <- crosspred(cb.foehn, mod_nm2, at=0:288, cumul=FALSE, cen = 0)

  # extract prediction for value 72 and save it
  table_estimates[2,i] = paste0( sprintf("%.3f",round(pred_nm2$allRRfit["72"],digits = 3)),
                                 " [",
                                 sprintf("%.3f",round(pred_nm2$allRRlow["72"],digits=3)),
                                 "-",
                                 sprintf("%.3f",round(pred_nm2$allRRhigh["72"], digits = 3)), "]")



}

# switch cols to rows
table_estimates <- t(table_estimates)

# assign proper col and rownames
colnames(table_estimates) <- c("RR Model 1 [CI]", "RR Model 2 [CI]")
rownames(table_estimates) <- c("all", "male", "female", "<65 years", ">64 years", "cardiovascular", "respiratory", "infectious", "genitourinary", "mental")

# save table
write.csv(table_estimates, file = "output/tables/Table2.csv")

#----



### sTable 4
#----

# Results for the sensitivity analysis of different buffer sizes on all cause
# hospitalizations with 95% confidence intervals.

# empty table
sensitivity_table = data.frame(
  Model1= rep(NA,2),
  Model2=rep(NA,2))

rownames(sensitivity_table) = c("5km", "10km")

sizes = c(5000,10000)

for(i in 1:length(sizes)){

  # extract buffer size
  buffer = sizes[i]

  # read data
  data_sens = read.csv(paste0("/Volumes/FS/_ISPM/CCH/Tino/master_thesis/data/Medstat_hospitalizations_aggregated/hosp_buffer_", buffer, ".csv")) |>
    mutate(date = as.Date(date),
           stratum_dow = as.factor(stratum_dow),
           y64 = a014y + a1564y,
           o64 = a6574y + a7584y + a85plusy,
           f_id_binary = ifelse(f_id >= 72, 1, 0),
           station_year = paste0(station, year))

  # index to include only stratum that have hosp counts
  ind_dow_sens = tapply(data_sens$all, data_sens$stratum_dow, sum);

  # define crossbasis
  cb.temp <- crossbasis(data_sens$temp,
                        lag=21,
                        argvar=list(fun="ns", knots = quantile(data_sens$temp, c(.5,.9), na.rm=TRUE)),
                        arglag=list(fun="ns", knots = logknots(21,3)),
                        group = data_sens$station)

  cb.foehn <- crossbasis(data_sens$f_id,
                         lag = 3,
                         argvar = list(fun="lin"),
                         arglag = list(fun="integer"),
                         group = data_sens$station)

  # formula
  formula1 <- as.formula("all ~ cb.foehn")

  # model
  mod_nm1 <- gnm(formula1, data = data_sens,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow_sens>0)

  # prediction
  pred_nm1 <- crosspred(cb.foehn, mod_nm1, at=0:288, cumul=FALSE, cen = 0)

  # extract prediction for value 72 and save it
  sensitivity_table[i,1] = paste0( sprintf("%.3f",round(pred_nm1$allRRfit["72"],digits=3)),
                                   " [",
                                   sprintf("%.3f",round(pred_nm1$allRRlow["72"], digits=3)),
                                   "-",
                                   sprintf("%.3f",round(pred_nm1$allRRhigh["72"], digits = 3)), "]")

  # formula, foehn + temp
  formula2 <- as.formula("all ~ cb.foehn + cb.temp")

  # model
  mod_nm2 <- gnm(formula2, data = data_sens,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow_sens>0)

  # prediction
  pred_nm2 <- crosspred(cb.foehn, mod_nm2, at=0:288, cumul=FALSE, cen = 0)

  # extract prediction for value 72 and save it
  sensitivity_table[i,2] = paste0( sprintf("%.3f",round(pred_nm2$allRRfit["72"],digits = 3)),
                                   " [",
                                   sprintf("%.3f",round(pred_nm2$allRRlow["72"],digits=3)),
                                   "-",
                                   sprintf("%.3f",round(pred_nm2$allRRhigh["72"], digits = 3)), "]")

}


# save table
write.csv(sensitivity_table, file = "output/tables/sTable4.csv")

#----


### sTable 5
#----

# Results for the sensitivity analysis on different lag periods of all cause
# hospitalizations with 95% confidence intervals.

sensitivity_table = data.frame(
  Model1= rep(NA,4),
  Model2=rep(NA,4))

rownames(sensitivity_table) = c(1, 5, 10, 15)

for(i in 1:nrow(sensitivity_table)){

  # define lag period
  lag_sens <- rownames(sensitivity_table)[i]
  print(lag_sens)

  cb.foehn <- crossbasis(data$f_id,
                         lag = as.integer(lag_sens),
                         argvar = list(fun="lin"),
                         arglag = list(fun="integer"),
                         group = data$station)

  # formula
  formula1 <- as.formula("all ~ cb.foehn")

  # model
  mod_nm1 <- gnm(formula1, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

  # prediction
  pred_nm1 <- crosspred(cb.foehn, mod_nm1, at=0:288, cumul=FALSE, cen = 0)

  # extract prediction for value 72 and save it
  sensitivity_table[i,1] = paste0( sprintf("%.3f",round(pred_nm1$allRRfit["72"],digits=3)),
                                   " [",
                                   sprintf("%.3f",round(pred_nm1$allRRlow["72"], digits=3)),
                                   "-",
                                   sprintf("%.3f",round(pred_nm1$allRRhigh["72"], digits = 3)), "]")

  # formula, foehn + temp
  formula2 <- as.formula("all ~ cb.foehn + cb.temp")

  # model
  mod_nm2 <- gnm(formula2, data = data,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow>0)

  # prediction
  pred_nm2 <- crosspred(cb.foehn, mod_nm2, at=0:288, cumul=FALSE, cen = 0)

  # extract prediction for value 72 and save it
  sensitivity_table[i,2] = paste0( sprintf("%.3f",round(pred_nm2$allRRfit["72"],digits = 3)),
                                   " [",
                                   sprintf("%.3f",round(pred_nm2$allRRlow["72"],digits=3)),
                                   "-",
                                   sprintf("%.3f",round(pred_nm2$allRRhigh["72"], digits = 3)), "]")

}

write.csv(sensitivity_table, file = "output/tables/sTable5.csv")

#----


### sTable 6
#----

# Sensitivity analysis on the foehn wind aggregation method: only full foehn wind considered

data_ff = read.csv("/Volumes/FS/_ISPM/CCH/Tino/master_thesis/data/Medstat_hospitalizations_aggregated/hosp_buffer_8000sensitivity_onlyfullfoehnaggregation.csv") |>
  mutate(station_date = paste0(station, date),
         date = as.Date(date),
         station = as.factor(station),
         stratum_dow = as.factor(stratum_dow),
         stratum = as.factor(stratum),
         y64 = a014y + a1564y,
         o64 = a6574y + a7584y + a85plusy)

# index to include only stratum that have hosp counts
ind_dow_ff = tapply(data_ff$all, data_ff$stratum_dow, sum)


cb.temp_ff<- crossbasis(data_ff$temp,
                        lag=21,
                        argvar=list(fun="ns", knots = quantile(data_ff$temp, c(.5,.9), na.rm=TRUE)),
                        arglag=list(fun="ns", knots = logknots(21,3)),
                        group = data_ff$station)

cb.foehn_ff <- crossbasis(data_ff$f_id_sens,
                          lag = 3,
                          argvar = list(fun="lin"),
                          arglag = list(fun="integer"),
                          group = data_ff$station)


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
                             row.names = c("RR Model 1 [CI]", "RR Model 2 [CI]"))

for (i in 1:length(groups_id)) {

  colvar = groups_id[i]

  # formula, only foehn
  formula1 <- as.formula(paste(colvar, "~ cb.foehn_ff"))

  # model
  mod_nm1 <- gnm(formula1, data = data_ff,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow_ff>0)

  # prediction
  pred_nm1 <- crosspred(cb.foehn_ff, mod_nm1, at=0:288, cumul=FALSE, cen = 0)

  # extract prediction for value 72 and save it
  table_estimates[1,i] = paste0( sprintf("%.3f",round(pred_nm1$allRRfit["72"],digits=3)),
                                 " [",
                                 sprintf("%.3f",round(pred_nm1$allRRlow["72"], digits=3)),
                                 "-",
                                 sprintf("%.3f",round(pred_nm1$allRRhigh["72"], digits = 3)), "]")



  # formula, foehn + temp
  formula2 <- as.formula(paste(colvar, "~ cb.foehn_ff + cb.temp_ff"))

  # model
  mod_nm2 <- gnm(formula2, data = data_ff,  family=quasipoisson(), eliminate=stratum_dow, subset=ind_dow_ff>0)

  # prediction
  pred_nm2 <- crosspred(cb.foehn_ff, mod_nm2, at=0:288, cumul=FALSE, cen = 0)

  # extract prediction for value 72 and save it
  table_estimates[2,i] = paste0( sprintf("%.3f",round(pred_nm2$allRRfit["72"],digits = 3)),
                                 " [",
                                 sprintf("%.3f",round(pred_nm2$allRRlow["72"],digits=3)),
                                 "-",
                                 sprintf("%.3f",round(pred_nm2$allRRhigh["72"], digits = 3)), "]")



}

# transpose the table
table_estimates_transposed = t(table_estimates)

# assign proper rownames
rownames(table_estimates_transposed) <- c("all", "male", "female", "<65 years", ">64 years", "cardiovascular", "respiratory", "infectious", "urinary", "mental")

# save table
write.csv(table_estimates_transposed, "output/tables/sTable6.csv")

#----











