




server <- function(input, output) {



  output$selected_radius <- renderText({
    paste0("You have selected this", input$buffer_radius)
  })









}
