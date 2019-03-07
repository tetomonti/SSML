partitions <-
    base::c("30L+10", "30L+20", "30L+30", "30L+40", "30L+50", "30L+60",
            "30L+70")

auc_results <-
    base::data.frame()

for (i in 1:7) {
    auc <-
        magrittr::extract2(results, i) %>%
        magrittr::use_series(training_auc)

    partition <-
        magrittr::extract(partitions, i)

    set <-
        as.character("Training")

    training_results <-
        base::data.frame(auc, partition, set)

    auc <-
        magrittr::extract2(results, i) %>%
        magrittr::use_series(testing_auc)

    partition <-
        magrittr::extract(partitions, i)

    set <-
        as.character("Testing")

    testing_results <-
        base::data.frame(auc, partition, set)

    auc_results <-
        base::rbind(auc_results, training_results, testing_results)

}

dplyr::filter(auc_results, set == "Training") %>%
    dplyr::group_by(partition) %>%
    dplyr::summarize(
        `Min.` = min(auc),
        `1st Qu.` = quantile(auc, 0.25),
        Median = median(auc),
        Mean = mean(auc),
        `3rd Qu.` = quantile(auc, 0.75),
        Max = max(auc)
    ) %>%
    dplyr::rename(Partition = partition) %>%
    knitr::kable(digits = 3)

dplyr::filter(auc_results, set == "Testing") %>%
    dplyr::group_by(partition) %>%
    dplyr::summarize(
        `Min.` = min(auc),
        `1st Qu.` = quantile(auc, 0.25),
        Median = median(auc),
        Mean = mean(auc),
        `3rd Qu.` = quantile(auc, 0.75),
        Max = max(auc)
    ) %>%
    dplyr::rename(Partition = partition) %>%
    knitr::kable(digits = 3)

dplyr::filter(auc_results, set == "Training") %>%
    ggplot2::ggplot(ggplot2::aes(x = partition, y = auc)) +
    ggplot2::stat_boxplot(geom = "errorbar") +
    ggplot2::geom_boxplot(notch = TRUE) +
    ggplot2::theme_linedraw() +
    ggplot2::labs(
        title = "Label Propagation using Gaussian Random Fields",
        subtitle = "Training, Top 5K Features by MAD, 100 Bootstrap Samples",
        x = NULL, y = "AUC"
    )

dplyr::filter(auc_results, set == "Testing") %>%
    ggplot2::ggplot(ggplot2::aes(x = partition, y = auc)) +
    ggplot2::stat_boxplot(geom = "errorbar") +
    ggplot2::geom_boxplot(notch = TRUE) +
    ggplot2::theme_linedraw() +
    ggplot2::labs(
        title = "Label Propagation using Gaussian Random Fields",
        subtitle = "Testing, Top 5K Features by MAD, 100 Bootstrap Samples",
        x = NULL, y = "AUC"
    )
