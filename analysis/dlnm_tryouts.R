################################################################################
# Package DLNM exploration
# - load the data and correct NAs (how?) ONLY USE COMPLETE.CASES()
# - perform crossbasis function for  foehn TODO and temp (done by Ana)
# - calculate model
# - predict using the model
# - plot the prediction in all possible ways

# Notes for plots that could be nice:
# - Overall effect of foehn: foehn vs RR with all lags included
# - effect of foehn on lag day 1-3
# - convert RR to percentages (easier to read)


# packages ####
library(dlnm)
library(splines)
library(ggplot2)
library(viridis)
library(gnm)
#####


# load data (and shorten it for tryouts) ####
rm(list=ls())
data = read.csv("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/data/MedStat_aggregated/hosp_buffer_5000.csv")

# HANDLE NAs HERE SET TO 0
data$f_id[is.na(data$f_id)] <- 0
data_short = data[1:5114,]
#####








# crossbasis functions for foehn and temp ####
cb.foehn <- crossbasis(data$f_id,                       # input variable
                        lag=3,                                # maximum lag to be considered
                        argvar=list(fun="lin"),  # defines how you want to model the nonlinear association between the exposure variable and the outcome.
                        arglag=list(fun="integer"))   # This list defines how the effect of the exposure changes over time (the lag structure)


### CROSSBASIS FOEHN #####
## lag = 3 is fixed

# linear
cb.foehn <- crossbasis(data$f_id,
                       lag=3,
                       argvar=list(fun="lin"),
                       arglag=list(fun="integer"))


#####

### MODEL ######

mod_nm <- gnm(all ~ cb.foehn, data = data,  family=quasipoisson(), eliminate=stratum, subset=ind>0)

# mod_lm <- glm(all ~ cb.foehn + ns(X, 7*14) + dow,  family=quasipoisson(), data_short)
#####

### PREDICTION ####
pred_nm <- crosspred(cb.foehn,
                     mod_nm,
                     at=0:288,
                     cumul=FALSE,
                     cen = 0)

#######









cb.temp <- crossbasis(data$temp,                       # input variable
                        lag=21,                               # maximum lag to be considered
                        argvar=list(fun="ns", knots = quantile(data$temp, c(.1,.75,.9), na.rm=TRUE)), # defines how you want to model the nonlinear association between the exposure variable and the outcome.
                        arglag=list(fun="ns", knots = logknots(21,3)))

# for multiple stations the group function cuts so that there is no continuous lag from
# one station to the other; the same for temp
# cb.foehn_all <- crossbasis(data$f_id,
#                            lag=5,
#                            argvar=list(fun="poly", degree = 2),
#                            arglag=list(fun="poly",degree = 2),
#                            group = data$station)
# cb.temp_all <- crossbasis(data$temp,
#                            lag=5,
#                            argvar=list(fun="poly", degree = 2),
#                            arglag=list(fun="poly",degree = 2),
#                            group = data$station)

#####
data$date = as.Date(data$date)
data$year = as.factor(format(data$date, "%Y"))
data$month = as.factor(format(data$date, "%m"))
data$dow = as.factor(data$dow)

data$stratum = with(data, factor(paste(station, year, month, dow, sep="-")))

ind = tapply(data$all, data$stratum, sum)

# model for impact on all hospitalisations   ####
model <- gnm(all ~ cb.foehn + cb.temp, data = data,  family=quasipoisson(), eliminate=stratum, subset=ind>0)

model2 <- glm(all ~ cb.foehn + cb.temp + ns(X, 7*14) + dow,  family=quasipoisson(), data_short)
# with all data
# model_all <- glm(all ~ cb.foehn_all + cb.temp_all + ns(X, 7*14) + dow,
#              family=quasipoisson(), data)

#####


# predict the impact of foehn on hosp ####
pred <- crosspred(cb.foehn,                    # basis
                  model,                       # model based on crossbasis functions
                  at=0:288,                   # This specifies whether to make predictions by each individual lag or across cumulative lags. can also be T or F
                  cumul=FALSE,                 # cumulative effect of the exposure across all lags, instead of separating the effects by each lag.
                  cen = 0)                     # the value that the effect will be compared to
# all data
# pred <- crosspred(cb.foehn_all,
#                   model_all,
#                   at=0:288,
#                   bylag=0.2,
#                   cumul=FALSE,
#                   cen = 0)
#####



# Different plots for visualization ####

plot(pred,                 ## Overall cumulative exposure-response
     "overall",
     col = 2,
     xlab = "Exposure",
     ylab = "Cumulative Response",
     lwd = 2,
     main = "Overall cumulative exposure-response")
# grid()



plot(pred,              ## exposure at specific exposure increase
     "slices",                       # Show exposure-response curves for specific lags
     var = 50,                       # For a specific exposure level
     #ci = "lines",                   # Use lines for confidence intervals
     col = 2,                        # Curve color
     ci.arg = list(col = rgb(0.5,0.5,0.5,0.3)),    # Uncomment if you want to customize CI color
     xlab = "Exposure (Foehn)",      # Custom x-axis label
     ylab = "Relative Risk (RR)",    # Custom y-axis label
     main = "Exposure-Response at Exposure Increase of 50",  # Custom title
     lwd = 2,                        # Line width for the curve
     ylim = c(0.9, 1.1)              # Control y-axis limits to focus on RR range
)


plot(pred,              ## exposure at specific lag
     "slices",                       # Show exposure-response curves for specific lags
     lag  = 0,                       # For a specific lg level
     ci = "area",                   # Use lines for confidence intervals
     col = 2,                        # Curve color
     ci.arg = list(col = rgb(0.5,0.5,0.5,0.3)),    # Uncomment if you want to customize CI color
     xlab = "Exposure (Foehn)",      # Custom x-axis label
     ylab = "Relative Risk (RR)",    # Custom y-axis label
     main = "Exposure-Response at Lag of 0",  # Custom title
     lwd = 3,                        # Line width for the curve
     ylim = c(0.85, 1.15)            # Control y-axis limits to focus on RR range
)



plot(pred,            ## 3D Plot
     "3d",
     main ="3-D Plot",
     xlab = "Exposure (Foehn)",
     ylab = "Lag (Days)",
     zlab = "Relative Risk",
     theta = 25,                    # The azimuthal angle (horizontal rotation)
     phi = 30,                      # The colatitude (vertical rotation).
     col = "skyblue2",              # color of the 3d surface
     border = rgb(0.5,0.5,0.5,0.5), # color of the borders
     #ticktype = "simple",           # type of grid on surface, alt.: "detailed"
     #nticks = 10                    # number of ticks
     r = 1.5
     )


# tryouts for hashed background
plot(pred,              ## exposure at specific exposure increase
     "slices",                       # Show exposure-response curves for specific lags
     var = 50,                       # For a specific exposure level
     #ci = "lines",                   # Use lines for confidence intervals
     col = 2,                        # Curve color
     ci.arg = list(density = 20, col = 4 ,angle = -45),    # Uncomment if you want to customize CI color
     xlab = "Exposure (Foehn)",      # Custom x-axis label
     ylab = "Relative Risk (RR)",    # Custom y-axis label
     main = "Exposure-Response at Exposure Increase of 50",  # Custom title
     lwd = 2,                        # Line width for the curve
     ylim = c(0.9, 1.1)              # Control y-axis limits to focus on RR range
)
lines(pred, ci = "area", var = 20, ci.arg = list(density = 20, col = 1), lwd = 2)

## I AM NOT HAPPY WITH MY CONTOUR PLOTS
plot(pred, "contour",               ## contour plot
     main = "Exposure-Lag-Response Contour Plot",  # Title
     xlab = "Exposure (Foehn)",                   # Custom x-axis label
     ylab = "Lag (Days)",                         # Custom y-axis label
     label = TRUE
     )



# Extract the matrix of predicted values for the contour plot
response_surface <- pred$matfit # Predicted matrix (response)

# Define the axes for exposure (x) and lag (y)
exposure <- pred$predvar
lag <- 0:25

# Generate the contour plot with labels
contour(exposure, lag, response_surface,
        main = "Exposure-Lag-Response Contour Plot",  # Title
        xlab = "Exposure (Foehn)",                   # Custom x-axis label
        ylab = "Lag (Days)",                         # Custom y-axis label
        nlevels = 10,                                # Number of contour levels
        col = heat.colors(10),                       # Color of the contour lines
        lwd = 2,                                     # Line width
        labcex = 0.8)




     #labels = TRUE)

     #col = heat.colors(10))

# dev.off()
#####




