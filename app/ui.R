# in ui.R

library(bslib)

ui <- bslib::page_sidebar(

  # App title
  title = "Foehn as a modifier of the temperature hospitalization association",

  # Buffer radius selection (in meters)
  sidebar = sidebar(
    helpText("The buffer radius determines the MedStat regions that will be included in the further analysis"),

    selectInput("buffer_radius",
                label = "Choose a buffer radius [m]",
                choices = list(2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000),
                selected = 8000),

    helpText("From experience: 8000m buffer radius shows the largest effect of foehn")
  ,

  # Checkbox for generating descriptive statistics
  checkboxInput("generate_stats", label = "Should I generate descriptive statistics?", value = FALSE)),

  # Main content area for plot output
  card(
    card_header("Figure of the selected MedStat regions"),
    plotOutput("selected_radius_plot")  # Buffer radius plot
  ),

  # Card for descriptive statistics (histograms)
  card(
    card_header("Descriptive Statistics - Histograms"),
    conditionalPanel(
      condition = "input.generate_stats == true",  # Display only if checkbox is checked
      plotOutput("histograms")  # Output placeholder for histograms
    )
  )
)
