
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

###
### Load libraries
################################################################################
library(shiny)

###
### Actual UI code
################################################################################
shinyUI(fluidPage(

    ## Application title (More space if not used).
    titlePanel(title = NULL, windowTitle = "PS Trimming in Three Groups"),

    ## Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel = sidebarPanel(
            width = 2,

            ## Scaling shiny plots to window height
            ## https://stackoverflow.com/questions/26782041/scaling-shiny-plots-to-window-height
            tags$head(tags$style("#three_trim_plot{height:90vh !important;}")),

            ## Sample size
            tags$div(class = "header", checked = NA,
                     tags$h5("Sample size")),
            sliderInput("n",
                        NULL,
                        min = 1,
                        max = 50000,
                        value = 10000),

            ## Dirichlet(alpha)
            ## adding the new div tag to the sidebar
            ## https://shiny.rstudio.com/articles/html-tags.html
            tags$div(class = "header", checked = NA,
                     tags$h5("Concentration")),
            sliderInput("alpha_overall",
                        NULL,
                        min = 0.01,
                        max = 20,
                        value = 2),
            tags$div(class = "header", checked = NA,
                     tags$h5("Relative groups sizes")),
            ## https://stackoverflow.com/questions/20637248/shiny-4-small-textinput-boxes-side-by-side
            splitLayout(numericInput("alpha0",
                                     NULL,
                                     min = 0.01,
                                     max = 100,
                                     value = 1),
                        numericInput("alpha1",
                                     NULL,
                                     min = 0.01,
                                     max = 100,
                                     value = 1),
                        numericInput("alpha2",
                                     NULL,
                                     min = 0.01,
                                     max = 100,
                                     value = 1)),

            ## Opacity
            tags$div(class = "header", checked = NA,
                     tags$h5("Opacity for scatter plot")),
            splitLayout(sliderInput("overall_opacity",
                                    "Overall",
                                    min = 0,
                                    max = 1,
                                    value = 1),
                        sliderInput("trimmed_opacity",
                                    "Trimmed",
                                    min = 0,
                                    max = 1,
                                    value = 1)),

            ## Which score to plot
            tags$div(class = "header", checked = NA,
                     tags$h5("Which score?")),
            selectInput("plot_pi",
                        NULL,
                        choices = c("Propensity" = FALSE,
                                    "Preference" = TRUE)),

            ## What type of plot
            tags$div(class = "header", checked = NA,
                     tags$h5("Plot density?")),
            selectInput("plot_density",
                        NULL,
                        choices = c("No" = FALSE,
                                    "Yes" = TRUE)),

            ## Facetting
            tags$div(class = "header", checked = NA,
                     tags$h5("By group?")),
            selectInput("facet",
                        NULL,
                        choices = c("No" = FALSE,
                                    "Yes" = TRUE)),

            ## Trimming control
            tags$div(class = "header", checked = NA,
                     tags$h5("Multinomial trimming thresholds")),
            sliderInput("thres_crump_multi",
                        "Crump (PS scale)",
                        min = 0,
                        max = 0.3,
                        value = 0.07),
            sliderInput("thres_sturmer_multi",
                        "Sturmer (Quantile scale)",
                        min = 0,
                        max = 0.3,
                        value = 0.03),
            sliderInput("thres_walker_multi",
                        "Walker (Preference scale)",
                        min = 0,
                        max = 0.3,
                        value = 0.2)

        ),

        ## Show a plot of the generated distribution
        mainPanel = mainPanel(width = 10,
                              plotOutput(outputId = "three_trim_plot",
                                         width = "auto",
                                         height = "700px"),
                              textOutput(outputId = "prevalence")
        )
    )
))
