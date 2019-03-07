partitions <-
    base::c("30L+10", "30L+20", "30L+30", "30L+40", "30L+50", "30L+60",
            "30L+70")

roc_results <-
    base::data.frame()

for (i in 1:7) {
    specificity <-
        magrittr::extract2(results, i) %>%
        magrittr::use_series(training_specificities) %>%
        base::sort(decreasing = TRUE)

    sensitivity <-
        magrittr::extract2(results, i) %>%
        magrittr::use_series(training_sensitivities) %>%
        base::sort(decreasing = FALSE)

    partition <-
        magrittr::extract(partitions, i)

    set <-
        as.character("Training")

    training_results <-
        base::data.frame(specificity, sensitivity, partition, set) %>%
        dplyr::add_row(specificity = 1, sensitivity = 0, partition = partition,
                       set = set) %>%
        dplyr::arrange(sensitivity)

    specificity <-
        magrittr::extract2(results, i) %>%
        magrittr::use_series(testing_specificities) %>%
        base::sort(decreasing = TRUE)

    sensitivity <-
        magrittr::extract2(results, i) %>%
        magrittr::use_series(testing_sensitivities) %>%
        base::sort(decreasing = FALSE)

    partition <-
        magrittr::extract(partitions, i)

    set <-
        as.character("Testing")

    testing_results <-
        base::data.frame(specificity, sensitivity, partition, set) %>%
        dplyr::add_row(specificity = 1, sensitivity = 0, partition = partition,
                       set = set) %>%
        dplyr::arrange(sensitivity)

    roc_results <-
        base::rbind(roc_results, training_results, testing_results)

}

dplyr::filter(roc_results, set == "Training") %>%
    ggplot2::ggplot(ggplot2::aes(specificity, sensitivity, color = partition)) +
    ggplot2::geom_step() +
    ggplot2::scale_x_reverse() +
    ggplot2::geom_abline(intercept = 1) +
    ggplot2::theme_linedraw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
        title = "Label propagation using Gaussian Random Fields",
        subtitle = "Training, Top 5K Features by MAD, 100 Bootstrap Samples",
        x = "Specificity", y = "Sensitivity", color = NULL
    ) +
    ggplot2::scale_color_discrete(guide = ggplot2::guide_legend(nrow = 1))

dplyr::filter(roc_results, set == "Testing") %>%
    ggplot2::ggplot(ggplot2::aes(specificity, sensitivity, color = partition)) +
    ggplot2::geom_step() +
    ggplot2::scale_x_reverse() +
    ggplot2::geom_abline(intercept = 1) +
    ggplot2::theme_linedraw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
        title = "Label propagation using Gaussian Random Fields",
        subtitle = "Testing, Top 5K Features by MAD, 100 Bootstrap Samples",
        x = "Specificity", y = "Sensitivity", color = NULL
    ) +
    ggplot2::scale_color_discrete(guide = ggplot2::guide_legend(nrow = 1))
