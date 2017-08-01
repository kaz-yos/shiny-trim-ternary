
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

###
### Load libraries
################################################################################
library(shiny)
## Added
library(tidyverse)
library(assertthat)
library(MCMCpack)
library(plot3D)
library(ggtern)


###
### Define functions
################################################################################
###  Data generator
## Multinomial treatment generator.
rmultinom1 <- function(n, prob) {
    assert_that(is.vector(prob))
    mat_group_x_n <- rmultinom(n = n, size = 1, prob = prob)
    group_names <- rownames(mat_group_x_n)
    if (is.null(group_names)) {
        group_names <- as.character(seq_len(nrow(mat_group_x_n)) - 1)
    }
    ## Return group names chosen
    sapply(seq_len(ncol(mat_group_x_n)), function(i) {
        group_names[mat_group_x_n[,i] == 1]
    })
}
## Version that work on a matrix of probability vectors.
rmultinom1_prob_mat <- function(prob_mat) {
    ## Assume n x group matrix of probability.
    n <- nrow(prob_mat)

    ## Loop over rows of prob_mat matrix.
    ## Each row is a vector of probabilities.
    sapply(seq_len(nrow(prob_mat)), function(i) {
        prob_vec <- prob_mat[i,]
        rmultinom1(n = 1, prob = prob_vec)
    })
}

## Generate n Dirichlet(alpha) data points.
generate_dirichlet_data <- function(n, alpha) {
    ## Propensity score vector from Dirichlet(alpha)
    ps_mat <- rdirichlet(n = n, alpha = alpha)
    colnames(ps_mat) <- as.character(seq_len(ncol(ps_mat)) - 1)
    ## Sample treatment assignment based on PS1
    A <- rmultinom1_prob_mat(ps_mat)
    ## Common part
    df <- data_frame(id  = seq_len(n),
                     type = "raw",
                     A   = A)
    ## Name PS appropriately
    df_ps <- as.data.frame(ps_mat)
    colnames(df_ps) <- paste0("ps", colnames(ps_mat))
    ## Combine
    bind_cols(df, df_ps)
}

## Generate n Dirichlet(alpha) data points with or without margins.
generate_dirichlet_plot_data <- function(n, alpha, add_margin = TRUE) {
    ## Get data frame
    df <- generate_dirichlet_data(n = n, alpha = alpha)

    if (add_margin) {
        ## Groups
        group_names <- as.character(seq_len(length(alpha)) - 1)
        n_groups <- length(alpha)

        margin_df <- lapply(group_names, function(group) {
            ## Subset to a group
            df_group <- df[df$A == group,]
            ## Extract group-specific PS
            ps_group <- df_group[,paste0("ps",group)]
            ## Zeroall PS
            df_group[,Filter(f = function(elt) {grepl("^ps", elt)}, x = names(df_group))] <- 0
            ## Recover
            df_group[,paste0("ps",group)] <- ps_group
            ## Change type
            df_group$type <- group
            ## Return df
            df_group
        }) %>% bind_rows

        return(bind_rows(df, margin_df))

    } else {

        return(df)
    }
}


###  Trimming functions
## Function to check being in a range
in_range <- function(x, l, u) {
    l <= x & x <= u
}

###   Crump
## Receive three PS, return "keep" indicator
trim_crump_multi <- function(A, ps0, ps1, ps2, thres = 0.1) {
    ## All of them has to be in [thres, 1.0].
    ## Treatment vector is not used.
    as.numeric(ps0 >= thres & ps1 >= thres & ps2 >= thres)
}

trim_crump_pair <- function(A, ps0, ps1, ps2, thres = 0.1, all_three = TRUE) {

    l <- thres
    u <- 1 - thres

    ## Pairwise PS
    ps1v0 <- ps1 / (ps1 + ps0)
    ps2v0 <- ps2 / (ps2 + ps0)
    ps2v1 <- ps2 / (ps2 + ps1)

    ## Pick if we use all three or group specific two.
    if (all_three) {

        return(as.numeric(in_range(ps1v0, l, u) &
                          in_range(ps2v0, l, u) &
                          in_range(ps2v1, l, u)))

    } else {

        return(as.numeric((A == 0 & in_range(ps1v0, l, u) & in_range(ps2v0, l, u)) |
                          (A == 1 & in_range(ps1v0, l, u) & in_range(ps2v1, l, u)) |
                          (A == 2 & in_range(ps2v0, l, u) & in_range(ps2v1, l, u))))
    }
}

###   Sturmer
trim_sturmer_multi <- function(A, ps0, ps1, ps2, thres = 0.05) {

    ## Calculate the lower bounds as thres-quantile in respective group.
    l0 <- quantile(ps0[A == 0], prob = thres)
    l1 <- quantile(ps1[A == 1], prob = thres)
    l2 <- quantile(ps2[A == 2], prob = thres)

    ## All three must be satisfied
    as.numeric(ps0 >= l0 & ps1 >= l1 & ps2 >= l2)
}

trim_sturmer_pair <- function(A, ps0, ps1, ps2, thres = 0.05, all_three = TRUE) {

    ## Pairwise PS
    ps1v0 <- ps1 / (ps1 + ps0)
    ps2v0 <- ps2 / (ps2 + ps0)
    ps2v1 <- ps2 / (ps2 + ps1)

    ## Calculate the lower bounds as thres-quantile in respective group.
    l1v0 <- quantile(ps1v0[A == 1], prob = thres)
    l2v0 <- quantile(ps2v0[A == 2], prob = thres)
    l2v1 <- quantile(ps2v1[A == 2], prob = thres)

    ## Calculate the upper bounds as (1 - thres)-quantile in respective group.
    u1v0 <- quantile(ps1v0[A == 0], prob = 1 - thres)
    u2v0 <- quantile(ps2v0[A == 0], prob = 1 - thres)
    u2v1 <- quantile(ps2v1[A == 1], prob = 1 - thres)

    ## Pick if we use all three or group specific two.
    if (all_three) {

        return(as.numeric(in_range(ps1v0, l1v0, u1v0) &
                          in_range(ps2v0, l2v0, u2v0) &
                          in_range(ps2v1, l2v1, u2v1)))

    } else {

        return(as.numeric((A == 0 & in_range(ps1v0, l1v0, u1v0) & in_range(ps2v0, l2v0, u2v0)) |
                          (A == 1 & in_range(ps1v0, l1v0, u1v0) & in_range(ps2v1, l2v1, u2v1)) |
                          (A == 2 & in_range(ps2v0, l2v0, u2v0) & in_range(ps2v1, l2v1, u2v1))))
    }
}

###   Walker
trim_walker_multi <- function(A, ps0, ps1, ps2, thres = 0.3) {

    ## Prevalence of treatment
    p0 <- mean(A == 0)
    p1 <- mean(A == 1)
    p2 <- mean(A == 2)

    ## Preference scores
    pi0 <- (ps0 / p0) / ((ps0 / p0) + (ps1 / p1) + (ps2 / p2))
    pi1 <- (ps1 / p1) / ((ps0 / p0) + (ps1 / p1) + (ps2 / p2))
    pi2 <- (ps2 / p2) / ((ps0 / p0) + (ps1 / p1) + (ps2 / p2))

    ## All three must be satisfied
    as.numeric(pi0 >= thres & pi1 >= thres & pi2 >= thres)
}

trim_walker_pair <- function(A, ps0, ps1, ps2, thres = 0.3, all_three = TRUE) {

    ## Prevalence of treatment
    p0 <- mean(A == 0)
    p1 <- mean(A == 1)
    p2 <- mean(A == 2)
    ## Pairwise prevalence
    p1v0 <- p1 / (p1 + p0)
    p2v0 <- p2 / (p2 + p0)
    p2v1 <- p2 / (p2 + p1)

    ## Pairwise PS
    ps1v0 <- ps1 / (ps1 + ps0)
    ps2v0 <- ps2 / (ps2 + ps0)
    ps2v1 <- ps2 / (ps2 + ps1)

    ## Pairwise preference scores
    pi1v0 <- (ps1v0 / p1v0) / ((ps1v0 / p1v0) + ((1 - ps1v0) / (1 - p1v0)))
    pi2v0 <- (ps2v0 / p2v0) / ((ps2v0 / p2v0) + ((1 - ps2v0) / (1 - p2v0)))
    pi2v1 <- (ps2v1 / p2v1) / ((ps2v1 / p2v1) + ((1 - ps2v1) / (1 - p2v1)))

    ## Pick if we use all three or group specific two.
    if (all_three) {

        return(as.numeric(in_range(pi1v0, thres, 1 - thres) &
                          in_range(pi2v0, thres, 1 - thres) &
                          in_range(pi2v1, thres, 1 - thres)))

    } else {

        return(as.numeric((A == 0 & in_range(pi1v0, thres, 1 - thres) & in_range(pi2v0, thres, 1 - thres)) |
                          (A == 1 & in_range(pi1v0, thres, 1 - thres) & in_range(pi2v1, thres, 1 - thres)) |
                          (A == 2 & in_range(pi2v0, thres, 1 - thres) & in_range(pi2v1, thres, 1 - thres))))

    }
}


###  Plot functions
plot_3d_trimmng <- function(data, trim_fun,
                            show_trimmed = 0.1, all_three = TRUE, plot_pi = FALSE, facet = FALSE) {

    if (plot_pi) {
        ## Add preference scores if plotting them.
        ## Prevalence of treatment
        p0 <- mean(data$A == 0)
        p1 <- mean(data$A == 1)
        p2 <- mean(data$A == 2)
        ## Preference scores
        common_denom <- ((data$ps0 / p0) + (data$ps1 / p1) + (data$ps2 / p2))
        data$pi0 <- (data$ps0 / p0) / common_denom
        data$pi1 <- (data$ps1 / p1) / common_denom
        data$pi2 <- (data$ps2 / p2) / common_denom
    }

    ## Add keep vector
    data$keep <- trim_fun(data$A, data$ps0, data$ps1, data$ps2)

    ## Proportion remaining
    p_kept <- mean(data$keep)

    ## Plot either preference score or PS
    if (plot_pi) {

        ggtern_plot <- ggtern(filter(data, keep == 1),
                              mapping = aes(x = pi1, y = pi0, z = pi2,
                                            color = factor(A))) +
            geom_point()

    } else {

        ggtern_plot <- ggtern(filter(data, keep == 1),
                              mapping = aes(x = ps1, y = ps0, z = ps2,
                                            color = factor(A))) +
            geom_point()

    }

    ## Add those who are trimmed if any exists and asked for.
    if (nrow(filter(data, keep == 0)) > 0 & show_trimmed > 0) {

        ggtern_plot <- ggtern_plot +
            geom_point(data = filter(data, keep == 0),
                       alpha = show_trimmed)

    }

    ## Clean up
    ggtern_plot <- ggtern_plot +
        scale_L_continuous(limits = c(0,1),
                           breaks = seq(from = 0.1, to = 0.9, by = 0.1),
                           labels = rep("", 9)) +
        scale_T_continuous(limits = c(0,1),
                           breaks = seq(from = 0.1, to = 0.9, by = 0.1),
                           labels = rep("", 9)) +
        scale_R_continuous(limits = c(0,1),
                           breaks = seq(from = 0.1, to = 0.9, by = 0.1),
                           labels = rep("", 9)) +
        scale_color_discrete(guide = FALSE) +
        labs(x = "1", y = "0", z = "2")

    ## Facet if asked
    if (facet) {
        ggtern_plot <- ggtern_plot +
            facet_grid(A ~ .)
    }

    ## Return the plot object
    return(ggtern_plot)
}
## Plot 6-panels
plot_3d_together <- function(data, show_trimmed = 0.1,
                             thres_crump_multi = 0.07,
                             thres_sturmer_multi = 0.03,
                             thres_walker_multi = 0.2,
                             thres_crump_pair = 0.1,
                             thres_sturmer_pair = 0.05,
                             thres_walker_pair = 0.3,
                             all_three = TRUE,
                             plot_pi = FALSE,
                             facet = FALSE) {

    ## Multinomial
    p_crump_multi <-
        plot_3d_trimmng(data = data,
                        trim_fun = partial(trim_crump_multi, thres = thres_crump_multi),
                        show_trimmed = show_trimmed, plot_pi = plot_pi, facet = facet) +
        labs(title = "Multinomial Crump")
    p_sturmer_multi <-
        plot_3d_trimmng(data = data,
                        trim_fun = partial(trim_sturmer_multi, thres = thres_sturmer_multi),
                        show_trimmed = show_trimmed, plot_pi = plot_pi, facet = facet) +
        labs(title = "Multinomial Sturmer")
    p_walker_multi <-
        plot_3d_trimmng(data = data,
                        trim_fun = partial(trim_walker_multi, thres = thres_walker_multi),
                        show_trimmed = show_trimmed, plot_pi = plot_pi, facet = facet) +
        labs(title = "Multinomial Walker")

    ## ## Pairwise version
    ## p_crump_pair <-
    ##     plot_3d_trimmng(data = data,
    ##                     trim_fun = partial(trim_crump_pair, thres = thres_crump_pair, all_three = all_three),
    ##                     show_trimmed = show_trimmed, plot_pi = plot_pi, facet = facet) +
    ##     labs(title = "Pairwise Crump")
    ## p_sturmer_pair <-
    ##     plot_3d_trimmng(data = data,
    ##                     trim_fun = partial(trim_sturmer_pair, thres = thres_sturmer_pair, all_three = all_three),
    ##                     show_trimmed = show_trimmed, plot_pi = plot_pi, facet = facet) +
    ##     labs(title = "Pairwise Sturmer")
    ## p_walker_pair <-
    ##     plot_3d_trimmng(data = data,
    ##                     trim_fun = partial(trim_walker_pair, thres = thres_walker_pair, all_three = all_three),
    ##                     show_trimmed = show_trimmed, plot_pi = plot_pi, facet = facet) +
    ##     labs(title = "Pairwise Walker")


    ## Graph together
    ## https://cran.r-project.org/web/packages/gridExtra/vignettes/arrangeGrob.html
    ggtern::grid.arrange(p_crump_multi,
                         p_sturmer_multi,
                         p_walker_multi,
                         ## p_crump_pair,
                         ## p_sturmer_pair,
                         ## p_walker_pair,
                         ncol = 3)
}


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
                         show_trimmed        = input$show_trimmed,
                         thres_crump_multi   = input$thres_crump_multi,
                         thres_sturmer_multi = input$thres_sturmer_multi,
                         thres_walker_multi  = input$thres_walker_multi,
                         thres_crump_pair    = input$thres_crump_pair,
                         thres_sturmer_pair  = input$thres_sturmer_pair,
                         thres_walker_pair   = input$thres_walker_pair,
                         all_three           = input$all_three,
                         plot_pi             = input$plot_pi,
                         facet               = input$facet)

        ## Add title including prevalence.
        total_alpha <- input$alpha0 + input$alpha1 + input$alpha2
        overall_title <- sprintf("Prevalence: %.2f:%.2f:%.2f",
                                 input$alpha0 / total_alpha,
                                 input$alpha1 / total_alpha,
                                 input$alpha2 / total_alpha)
        title(main = overall_title, outer = TRUE, line = -1)

    })

})
