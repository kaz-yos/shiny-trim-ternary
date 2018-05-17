
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

###
### Load libraries
################################################################################
library(shiny)

source("../shiny_trim_empirical_ternary/server_functions.R")


###
### Actual server code
################################################################################
## input comes from the UI.
shinyServer(function(input, output) {

    ## Create a reactive dataset.
    ## https://shiny.rstudio.com/tutorial/lesson6/
    reactive_data_fun <- reactive({
        generate_dirichlet_plot_data(n = input$n,
                                     alpha = c(input$alpha0, input$alpha1, input$alpha2),
                                     add_margin = FALSE)
    })


    ## Add plot object named "three_trim_plot" to the output.
    output$three_trim_plot <- renderPlot({

        ## Plot using the reactive dataset.
        plot_3d_together(data                = reactive_data_fun(),
                         overall_opacity     = input$overall_opacity,
                         trimmed_opacity     = input$trimmed_opacity,
                         thres_crump_multi   = input$thres_crump_multi,
                         thres_sturmer_multi = input$thres_sturmer_multi,
                         thres_walker_multi  = input$thres_walker_multi,
                         thres_crump_pair    = input$thres_crump_pair,
                         thres_sturmer_pair  = input$thres_sturmer_pair,
                         thres_walker_pair   = input$thres_walker_pair,
                         all_three           = input$all_three,
                         plot_pi             = input$plot_pi,
                         plot_density        = input$plot_density,
                         facet               = input$facet)

        ## Add title including prevalence.
        if (input$facet == FALSE) {
            total_alpha <- input$alpha0 + input$alpha1 + input$alpha2
            overall_title <- sprintf("Prevalence: %.2f:%.2f:%.2f",
                                     input$alpha0 / total_alpha,
                                     input$alpha1 / total_alpha,
                                     input$alpha2 / total_alpha)
            title(main = overall_title, outer = TRUE, line = -1)
        }

    })

})
