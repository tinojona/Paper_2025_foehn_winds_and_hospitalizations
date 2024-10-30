################################################################################
# Here I run my app to compute different plots and data for my interaction analysis



# libraries
library(shiny)
library(bslib)


# source ui and server
source("app/ui.R")
source("app/server.R")


# allocate the functions
shinyApp(ui = ui, server = server)


# run the app
runApp("app")
