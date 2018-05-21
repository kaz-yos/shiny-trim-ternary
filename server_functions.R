################################################################################
### Function Definitions for Shiny App
##
## Created on: 2017-11-14
## Author: Kazuki Yoshida
################################################################################

## Added
library(tidyverse)
library(assertthat)
library(MCMCpack)
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

###
###  Trimming functions
################################################################################
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

        return(as.numeric(between(ps1v0, l, u) &
                          between(ps2v0, l, u) &
                          between(ps2v1, l, u)))

    } else {

        return(as.numeric((A == 0 & between(ps1v0, l, u) & between(ps2v0, l, u)) |
                          (A == 1 & between(ps1v0, l, u) & between(ps2v1, l, u)) |
                          (A == 2 & between(ps2v0, l, u) & between(ps2v1, l, u))))
    }
}

###   Sturmer
trim_sturmer_multi <- function(A, ps0, ps1, ps2, thres = 0.05) {

    ## Calculate the lower bounds as thres-quantile in respective group.
    l0 <- quantile(ps0[A == 0], prob = thres)
    l1 <- quantile(ps1[A == 1], prob = thres)
    l2 <- quantile(ps2[A == 2], prob = thres)

    ## All three must be satisfied
    result <- as.numeric(ps0 >= l0 & ps1 >= l1 & ps2 >= l2)

    attributes(result) <- list(bounds = c(l0, l1, l2))

    result
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

        return(as.numeric(between(ps1v0, l1v0, u1v0) &
                          between(ps2v0, l2v0, u2v0) &
                          between(ps2v1, l2v1, u2v1)))

    } else {

        return(as.numeric((A == 0 & between(ps1v0, l1v0, u1v0) & between(ps2v0, l2v0, u2v0)) |
                          (A == 1 & between(ps1v0, l1v0, u1v0) & between(ps2v1, l2v1, u2v1)) |
                          (A == 2 & between(ps2v0, l2v0, u2v0) & between(ps2v1, l2v1, u2v1))))
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

        return(as.numeric(between(pi1v0, thres, 1 - thres) &
                          between(pi2v0, thres, 1 - thres) &
                          between(pi2v1, thres, 1 - thres)))

    } else {

        return(as.numeric((A == 0 & between(pi1v0, thres, 1 - thres) & between(pi2v0, thres, 1 - thres)) |
                          (A == 1 & between(pi1v0, thres, 1 - thres) & between(pi2v1, thres, 1 - thres)) |
                          (A == 2 & between(pi2v0, thres, 1 - thres) & between(pi2v1, thres, 1 - thres))))

    }
}


###
###  Plot functions
################################################################################
plot_3d_trimmng <- function(data, trim_fun, thres,
                            overall_opacity = 1.0, trimmed_opacity = 1.0,
                            all_three = TRUE, plot_pi = FALSE,
                            plot_density = TRUE,
                            facet = FALSE, title = NULL, show_title = TRUE) {

    crump   <- (title == "Multinomial Crump")
    sturmer <- (title == "Multinomial Sturmer")
    walker  <- (title == "Multinomial Walker")

    ## Prevalence of treatment
    p0 <- mean(data$A == 0)
    p1 <- mean(data$A == 1)
    p2 <- mean(data$A == 2)
    ## Preference scores
    common_denom <- ((data$ps0 / p0) + (data$ps1 / p1) + (data$ps2 / p2))
    data$pi0 <- (data$ps0 / p0) / common_denom
    data$pi1 <- (data$ps1 / p1) / common_denom
    data$pi2 <- (data$ps2 / p2) / common_denom

    ## Add keep vector
    ## This has bounds attribute in Sturmer
    data$keep <- trim_fun(data$A, data$ps0, data$ps1, data$ps2, thres)

    ## Proportion remaining
    p_kept <- mean(data$keep)
    p_kept_by_group <- tapply(data$keep, data$A, mean)

    ## Create a title with proportions kept
    if (show_title) {

        title <- paste0(title, sprintf(" %.2f (%.2f/%.2f/%.2f)",
                                       p_kept,
                                       p_kept_by_group[1],
                                       p_kept_by_group[2],
                                       p_kept_by_group[3]))
    } else {

        title <- NULL

    }

    ## Plot either preference score or PS
    if (plot_pi) {

        ggtern_plot <- ggtern(filter(data, keep == 1),
                              mapping = aes(x = pi1, y = pi0, z = pi2,
                                            color = factor(A)))


    } else {

        ggtern_plot <- ggtern(filter(data, keep == 1),
                              mapping = aes(x = ps1, y = ps0, z = ps2,
                                            color = factor(A)))

    }

    ## Add density
    if (plot_density) {

        ## Otherwise plot density
        ggtern_plot <- ggtern_plot +
            stat_density_tern(data = data,
                              geom = "density_tern",
                              bins = 10)

    } else {

        ## Plot points if not plotting density
        ggtern_plot <- ggtern_plot +
            geom_point(size = 0.1, alpha = overall_opacity)

        ## Add those who are trimmed if any exists and asked for.
        if (nrow(filter(data, keep == 0)) > 0 & trimmed_opacity > 0) {

            ggtern_plot <- ggtern_plot +
                geom_point(data = filter(data, keep == 0),
                           alpha = overall_opacity * trimmed_opacity,
                           size = 0.1)

        }

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

    ## Add boundaries
    if (plot_pi) {

        ## If plotting preference score

        if (crump) {

            ## Crump needs transformation
            ggtern_plot <- ggtern_plot +
                geom_polygon(data = data_frame(ps0 = c(thres, thres, (1 - 2 * thres)),
                                               ps1 = c(thres, (1 - 2 * thres), thres),
                                               ps2 = c((1 - 2 * thres), thres, thres),
                                               pi0 = (ps0/p0) / ((ps0/p0) + (ps1/p1) + (ps2/p2)),
                                               pi1 = (ps1/p1) / ((ps0/p0) + (ps1/p1) + (ps2/p2)),
                                               pi2 = (ps2/p2) / ((ps0/p0) + (ps1/p1) + (ps2/p2))),
                             mapping = aes(color = NULL),
                             fill = NA, color = "black")
        }

        if (sturmer) {

            bounds <- attributes(data$keep)$bounds
            l0 <- bounds[1]
            l1 <- bounds[2]
            l2 <- bounds[3]

            ## Sturmer needs transformation
            ggtern_plot <- ggtern_plot +
                geom_polygon(data = data_frame(ps0 = c((1 - (l1 + l2)), l0,              l0),
                                               ps1 = c(l1,              l1,              (1 - (l0 + l2))),
                                               ps2 = c(l2,              (1 - (l0 + l1)), l2),
                                               pi0 = (ps0/p0) / ((ps0/p0) + (ps1/p1) + (ps2/p2)),
                                               pi1 = (ps1/p1) / ((ps0/p0) + (ps1/p1) + (ps2/p2)),
                                               pi2 = (ps2/p2) / ((ps0/p0) + (ps1/p1) + (ps2/p2))),
                             mapping = aes(color = NULL),
                             fill = NA, color = "black")
        }

        if (walker) {

            ## Walker is straightforward
            ggtern_plot <- ggtern_plot +
                geom_polygon(data = data_frame(pi0 = c(thres, thres, (1 - 2 * thres)),
                                               pi1 = c(thres, (1 - 2 * thres), thres),
                                               pi2 = c((1 - 2 * thres), thres, thres)),
                             mapping = aes(color = NULL),
                             fill = NA, color = "black")
        }

    } else {

        ## If plotting propensity score

        if (crump) {

            ## Crump needs is straight forward
            ggtern_plot <- ggtern_plot +
                geom_polygon(data = data_frame(ps0 = c(thres, thres, (1 - 2 * thres)),
                                               ps1 = c(thres, (1 - 2 * thres), thres),
                                               ps2 = c((1 - 2 * thres), thres, thres)),
                             mapping = aes(color = NULL),
                             fill = NA, color = "black")
        }

        if (sturmer) {

            bounds <- attributes(data$keep)$bounds
            l0 <- bounds[1]
            l1 <- bounds[2]
            l2 <- bounds[3]

            ## Sturmer needs calculation
            ggtern_plot <- ggtern_plot +
                geom_polygon(data = data_frame(ps0 = c((1 - (l1 + l2)), l0,              l0),
                                               ps1 = c(l1,              l1,              (1 - (l0 + l2))),
                                               ps2 = c(l2,              (1 - (l0 + l1)), l2)),
                             mapping = aes(color = NULL),
                             fill = NA, color = "black")
        }

        if (walker) {

            ## Walker is needs transformation
            ggtern_plot <- ggtern_plot +
                geom_polygon(data = data_frame(pi0 = c(thres, thres, (1 - 2 * thres)),
                                               pi1 = c(thres, (1 - 2 * thres), thres),
                                               pi2 = c((1 - 2 * thres), thres, thres),
                                               ps0 = (pi0 * p0) / ((pi0 * p0) + (pi1 * p1) + (pi2 * p2)),
                                               ps1 = (pi1 * p1) / ((pi0 * p0) + (pi1 * p1) + (pi2 * p2)),
                                               ps2 = (pi2 * p2) / ((pi0 * p0) + (pi1 * p1) + (pi2 * p2))),
                             mapping = aes(color = NULL),
                             fill = NA, color = "black")
        }
    }


    ## Add title
    ggtern_plot <- ggtern_plot + labs(title = title)

    ## Return the plot object
    return(ggtern_plot)
}


## Plot 6-panels
plot_3d_together <- function(data,
                             overall_opacity = 1.0,
                             trimmed_opacity = 1.0,
                             thres_crump_multi = 0.07,
                             thres_sturmer_multi = 0.03,
                             thres_walker_multi = 0.2,
                             thres_crump_pair = 0.1,
                             thres_sturmer_pair = 0.05,
                             thres_walker_pair = 0.3,
                             all_three = TRUE,
                             plot_pi = FALSE,
                             plot_density = FALSE,
                             facet = FALSE) {

    ## Multinomial
    p_crump_multi <-
        plot_3d_trimmng(data = data,
                        trim_fun = trim_crump_multi, thres = thres_crump_multi,
                        overall_opacity = overall_opacity, trimmed_opacity = trimmed_opacity,
                        plot_pi = plot_pi, plot_density = plot_density,
                        facet = facet, title = "Multinomial Crump")
    p_sturmer_multi <-
        plot_3d_trimmng(data = data,
                        trim_fun = trim_sturmer_multi, thres = thres_sturmer_multi,
                        overall_opacity = overall_opacity, trimmed_opacity = trimmed_opacity,
                        plot_pi = plot_pi, plot_density = plot_density,
                        facet = facet, title = "Multinomial Sturmer")
    p_walker_multi <-
        plot_3d_trimmng(data = data,
                        trim_fun = trim_walker_multi, thres = thres_walker_multi,
                        overall_opacity = overall_opacity, trimmed_opacity = trimmed_opacity,
                        plot_pi = plot_pi, plot_density = plot_density,
                        facet = facet, title = "Multinomial Walker")

    ## Graph together
    ## https://cran.r-project.org/web/packages/gridExtra/vignettes/arrangeGrob.html
    ggtern::grid.arrange(p_crump_multi,
                         p_sturmer_multi,
                         p_walker_multi,
                         ncol = 3)
}
