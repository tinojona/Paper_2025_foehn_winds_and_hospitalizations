################################################################################
# Package DLNM exploration
# - load the data and correct NAs (how?)
# - perform crossbasis function for  foehn and temp
# - calculate model
# - predict using the model
# - plot the prediction in all possible ways

# Notes for plots that could be nice:
# - Overall effect of foehn: foehn vs RR with all lags included
# - effect of foehn on lag day 1-3
# - convert RR to percentages (easier to read)
#


# packages ####
library(dlnm)
library(splines)
library(ggplot2)
library(viridis)
#####


# load data (and shorten it for tryouts) ####
rm(list=ls())
data = read.csv("C:/Users/tinos/Documents/Master - Climate Science/3 - Master Thesis/master_thesis/data/MedStat_aggregated/hosp_buffer_5000.csv")

# HANDLE NAs HERE SET TO 0
data$f_id[is.na(data$f_id)] <- 0
data_short = data[1:5114,]
#####



# crossbasis functions for foehn and temp ####
cb.foehn <- crossbasis(data_short$f_id,                       # input variable
                        lag=5,                                # maximum lag to be considered
                        argvar=list(fun="poly", degree = 2),  # defines how you want to model the nonlinear association between the exposure variable and the outcome.
                        arglag=list(fun="poly",degree = 2))   # This list defines how the effect of the exposure changes over time (the lag structure)

cb.temp <- crossbasis(data_short$temp,                       # input variable
                        lag=5,                               # maximum lag to be considered
                        argvar=list(fun="poly", degree = 2), # defines how you want to model the nonlinear association between the exposure variable and the outcome.
                        arglag=list(fun="poly",degree = 2))

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


# model for impact on all hospitalisations   ####
model <- glm(all ~ cb.foehn + cb.temp + ns(X, 7*14) + dow,
              family=quasipoisson(), data_short)

# with all data
# model_all <- glm(all ~ cb.foehn_all + cb.temp_all + ns(X, 7*14) + dow,
#              family=quasipoisson(), data)

#####


# predict the impact of foehn on hosp ####
pred <- crosspred(cb.foehn,                    # basis
                  model,                       # model based on crossbasis functions
                  at=0:288,                    # values of the exposure at which to make predictions
                  bylag=0.2,                   # This specifies whether to make predictions by each individual lag or across cumulative lags. can also be T or F
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
     ylim = c(0.75,2),
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




