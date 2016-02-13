
# This is the user-interface definition of a Shiny web application.


library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("ECG Simulation"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    wellPanel(
            tags$style(type="text/css", '#leftPanel { width:200px; float:left;}'),
            id = "leftPanel",
            sliderInput("rate",
                        "Heart Rate",
                        min = 60,
                        max = 120,
                        value = 72),
            img(src="ecgaxis.png", width=160),
            sliderInput("qaxis",
                        "QRS axis (degrees)",
                        min=-90,
                        max=270,
                        value=60),
            sliderInput("paxis",
                        "P axis (degrees)",
                        min=-90,
                        max=270,
                        value=60),
            sliderInput("taxis",
                        "T axis (degrees)",
                        min=-90,
                        max=270,
                        value=60),
            
            checkboxInput("show_p", "Show P Wave", value=TRUE),
            checkboxInput("show_q", "Show QRS Wave", value=TRUE),
            checkboxInput("show_t", "Show T Wave", value=TRUE),
            checkboxInput("show_u", "Show U Wave", value=TRUE)
    ),

    # Show a plot of the generated distribution
    mainPanel(
            tabsetPanel(
            tabPanel("12-Leads", plotOutput("display12leads")),
            tabPanel("I", plotOutput("display_lead_I")),
            tabPanel("II", plotOutput("display_lead_II")),
            tabPanel("III", plotOutput("display_lead_III")),
            tabPanel("aVR", plotOutput("display_lead_aVR")),
            tabPanel("aVL", plotOutput("display_lead_aVL")),
            tabPanel("aVF", plotOutput("display_lead_aVF")),
            tabPanel("V1", plotOutput("display_lead_V1")),
            tabPanel("V2", plotOutput("display_lead_V2")),
            tabPanel("V3", plotOutput("display_lead_V3")),
            tabPanel("V4", plotOutput("display_lead_V4")),
            tabPanel("V5", plotOutput("display_lead_V5")),
            tabPanel("V6", plotOutput("display_lead_V6"))
##            ,tabPanel("DEBUG", textOutput("display_debug"))
            )
    )
  )
))
