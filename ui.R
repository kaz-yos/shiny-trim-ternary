
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

            ## What to plot
            tags$div(class = "header", checked = NA,
                     tags$h5("Which score?")),
            selectInput("plot_pi",
                        NULL,
                        choices = c("Propensity" = FALSE,
                                    "Preference" = TRUE)),

            ## Sample size
            tags$div(class = "header", checked = NA,
                     tags$h5("Sample size")),
            sliderInput("n",
                        NULL,
                        min = 1,
                        max = 10000,
                        value = 100),

            ## Dirichlet(alpha)
            ## adding the new div tag to the sidebar
            ## https://shiny.rstudio.com/articles/html-tags.html
            tags$div(class = "header", checked = NA,
                     tags$h5("Dirichlet parameters")),
            sliderInput("alpha0",
                        NULL,
                        min = 0.01,
                        max = 20,
                        value = 2),
            sliderInput("alpha1",
                        NULL,
                        min = 0.01,
                        max = 20,
                        value = 2),
            sliderInput("alpha2",
                        NULL,
                        min = 0.01,
                        max = 20,
                        value = 2),

            ## Opacity
            tags$div(class = "header", checked = NA,
                     tags$h5("Trim opacity")),
            sliderInput("show_trimmed",
                        NULL,
                        min = 0,
                        max = 1,
                        value = 0.1),

            ## Trimming control
            tags$div(class = "header", checked = NA,
                     tags$h5("Trimming thresholds")),
            sliderInput("thres_crump_multi",
                        "Crump, multi",
                        min = 0,
                        max = 1,
                        value = 0.07),
            sliderInput("thres_sturmer_multi",
                        "Sturmer, multi",
                        min = 0,
                        max = 1,
                        value = 0.03),
            sliderInput("thres_walker_multi",
                        "Walker, multi",
                        min = 0,
                        max = 1,
                        value = 0.2),
            sliderInput("thres_crump_pair",
                        "Crump, pair",
                        min = 0,
                        max = 1,
                        value = 0.1),
            sliderInput("thres_sturmer_pair",
                        "Sturmer, pair",
                        min = 0,
                        max = 1,
                        value = 0.05),
            sliderInput("thres_walker_pair",
                        "Walker, pair",
                        min = 0,
                        max = 1,
                        value = 0.3),

            ## Pairwise trimming rule control
            tags$div(class = "header", checked = NA,
                     tags$h5("Pairwise rule")),
            selectInput("all_three",
                        NULL,
                        choices = c("All three" = TRUE,
                                    "Two at a time" = FALSE))
        ),

        ## Show a plot of the generated distribution
        mainPanel = mainPanel(width = 10,
                              plotOutput(outputId = "three_trim_plot",
                                         width = "auto",
                                         height = "700px")
        )
    )
))
