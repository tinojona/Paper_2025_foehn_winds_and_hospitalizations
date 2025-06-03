################################################################################
# FOEHN WIND Crossbasis
################################################################################
# - determine the best crossbasis parameters to model foehn wind

# Conclusion
# - 8 km buffer shows the strongest response to foehn winds
# - linear exposure resonse has the best qAIC
# - integer lag response using the 3 day lag period has the best qAIC


### Data
#----

rm(list=ls())

# libraries
library(dlnm);library(splines);library(ggplot2);library(viridis);library(gnm);library(tidyverse)

# data
data = read.csv("../data/Medstat_hospitalizations_aggregated/hosp_buffer_8000.csv") |>
  mutate(date = as.Date(date),
         stratum_dow = as.factor(stratum_dow))

# index to include only stratum that have hosp counts
ind_dow = tapply(data$all, data$stratum_dow, sum)

# q-AIC computation
source("functions/qAIC_determination.R")

# define the maximum lag distance we account for
maxlago <- 3

#----


### ARGVAR ARGLAG DEFINITION
#----

# two lists of argvar and arglag arguments
v_var <- list(list(fun = "ns", knots = quantile(data$f_id, c(.8, .9), na.rm=TRUE),Boundary=range(data$f_id)),
              list(fun = "ns", knots = quantile(data$f_id, c(.8, .9, .95), na.rm=TRUE),Boundary=range(data$f_id)),
              list(fun = "strata", breaks = equalknots(data$f_id, nk = 3)),
              list(fun = "strata", breaks = equalknots(data$f_id, nk = 4)),
              list(fun = "strata", breaks = equalknots(data$f_id, nk = 5)),
              list(fun = "ns", knots = equalknots(data$f_id, nk=2) ,Boundary=range(data$f_id)),
              list(fun = "ns", knots = equalknots(data$f_id, nk=3) ,Boundary=range(data$f_id)),
              list(fun = "ns", knots = equalknots(data$f_id, nk=4) ,Boundary=range(data$f_id)),
              list(fun = "lin")
)

v_lag <- list(list(fun = "integer"),
              list(fun = "strata", breaks = 1),
              list(fun = "ns", knots = 1),
              list(fun = "ns", knots = c(1,2))
)

#----


### qAIC DETERMINATION
#----

# create an empty matrix to store the qAIC
qaic_tab <- matrix(NA,
                   nrow = length(v_var),
                   ncol=length(v_lag),
                   dimnames = list(c(v_var), c(v_lag)))


# determine the qAIC for each combination
for (i in 1:length(v_var)){

  # extract variable function
  argvar = v_var[[i]]

  for (j in 1:length(v_lag)) {

    #  extract lag function
    arglag = v_lag[[j]]

    # crossbasis
    cb.f_id <- crossbasis(data$f_id,
                          lag = maxlago,
                          argvar = argvar,
                          arglag = arglag,
                          group = data$station)

    # model
    mod <- gnm(all ~ cb.f_id,
               data = data,
               eliminate = stratum_dow,
               subset = ind_dow > 0,
               family = quasipoisson())

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

#----

