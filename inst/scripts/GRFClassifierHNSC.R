library(magrittr)

# base::load("data/miniHNSC.rda")
#
# n_unlabeled <-
#     base::c(32, 63, 95, 126, 158, 189, 221)
#
# n_bootstrap <-
#     base::seq(100)
#
# results <-
#     base::list()
#
# for (i in n_unlabeled) {
#     j <-
#         base::match(i, n_unlabeled)
#
#     training_auc <-
#         base::numeric()
#
#     training_specificities <-
#         base::list()
#
#     training_sensitivities <-
#         base::list()
#
#     testing_auc <-
#         base::numeric()
#
#     testing_specificities <-
#         base::list()
#
#     testing_sensitivities <-
#         base::list()
#
#     for (k in n_bootstrap) {
#         n_train <-
#             base::nrow(miniHNSC) %>%
#             magrittr::multiply_by(0.70) %>%
#             base::round()
#
#         i_train <-
#             base::nrow(miniHNSC) %>%
#             base::sample(n_train)
#
#         training_data <-
#             dplyr::slice(miniHNSC, i_train)
#
#         testing_data <-
#             dplyr::slice(miniHNSC, -i_train)
#
#         n_labeled <-
#             base::nrow(training_data) %>%
#             magrittr::multiply_by(0.30) %>%
#             base::round()
#
#         i_labeled <-
#             base::nrow(training_data) %>%
#             base::sample(n_labeled)
#
#         labeled_data <-
#             dplyr::slice(training_data, i_labeled) %>%
#             dplyr::select(-Grade)
#
#         labeled_labels <-
#             dplyr::slice(training_data, i_labeled) %>%
#             magrittr::use_series(Grade)
#
#         i_unlabeled <-
#             dplyr::slice(training_data, -i_labeled) %>%
#             base::nrow() %>%
#             base::seq() %>%
#             base::sample(i)
#
#         unlabeled_data <-
#             dplyr::slice(training_data, -i_labeled) %>%
#             dplyr::slice(i_unlabeled) %>%
#             dplyr::select(-Grade)
#
#         unlabeled_labels <-
#             dplyr::slice(training_data, -i_labeled) %>%
#             dplyr::slice(i_unlabeled) %>%
#             magrittr::use_series(Grade)
#
#         unlabeled_classes <-
#             base::as.integer(unlabeled_labels)
#
#         training_classifier <-
#             RSSL::GRFClassifier(X = labeled_data, y = labeled_labels,
#                                 X_u = unlabeled_data)
#
#         training_probabilities <-
#             methods::slot(training_classifier, "responsibilities")
#
#         training_seq <-
#             base::nrow(training_probabilities) %>%
#             base::seq()
#
#         unlabeled_probabilities <-
#             magrittr::extract(training_probabilities, training_seq, 1)
#
#         training_roc <-
#             pROC::roc(unlabeled_classes, unlabeled_probabilities)
#
#         training_auc[k] <-
#             pROC::auc(unlabeled_classes, unlabeled_probabilities) %>%
#             base::as.numeric()
#
#         training_specificities[[k]] <-
#             magrittr::use_series(training_roc, specificities)
#
#         training_sensitivities[[k]] <-
#             magrittr::use_series(training_roc, sensitivities)
#
#         predicted_labels <-
#             RSSL::predict(training_classifier)
#
#         predicted_data <-
#             dplyr::slice(training_data, -i_labeled) %>%
#             dplyr::slice(i_unlabeled) %>%
#             dplyr::mutate(Grade = predicted_labels)
#
#         joined_data <-
#             dplyr::slice(training_data, i_labeled) %>%
#             dplyr::bind_rows(predicted_data) %>%
#             dplyr::select(-Grade)
#
#         joined_labels <-
#             dplyr::slice(training_data, i_labeled) %>%
#             dplyr::bind_rows(predicted_data) %>%
#             magrittr::use_series(Grade)
#
#         testing_data <-
#             dplyr::slice(miniHNSC, -i_train) %>%
#             dplyr::select(-Grade)
#
#         testing_labels <-
#             dplyr::slice(miniHNSC, -i_train) %>%
#             magrittr::use_series(Grade)
#
#         testing_classes <-
#             base::as.integer(testing_labels)
#
#         testing_classifier <-
#             RSSL::GRFClassifier(X = joined_data, y = joined_labels,
#                                 X_u = testing_data)
#
#         testing_probabilities <-
#             methods::slot(testing_classifier, "responsibilities")
#
#         testing_seq <-
#             base::nrow(testing_probabilities) %>%
#             base::seq()
#
#         testing_probabilities <-
#             magrittr::extract(testing_probabilities, testing_seq, 1)
#
#         testing_roc <-
#             pROC::roc(testing_classes, testing_probabilities)
#
#         testing_auc[k] <-
#             pROC::auc(testing_classes, testing_probabilities) %>%
#             base::as.numeric()
#
#         testing_specificities[[k]] <-
#             magrittr::use_series(testing_roc, specificities)
#
#         testing_sensitivities[[k]] <-
#             magrittr::use_series(testing_roc, sensitivities)
#     }
#
#     training_specificities <-
#         purrr::reduce(training_specificities, base::cbind) %>%
#         base::apply(1, base::mean)
#
#     training_sensitivities <-
#         purrr::reduce(training_sensitivities, base::cbind) %>%
#         base::apply(1, base::mean)
#
#     testing_specificities <-
#         purrr::reduce(testing_specificities, base::cbind) %>%
#         base::apply(1, base::mean)
#
#     testing_sensitivities <-
#         purrr::reduce(testing_sensitivities, base::cbind) %>%
#         base::apply(1, base::mean)
#
#     results[[j]] <-
#         base::list(training_auc = training_auc,
#                    training_specificities = training_specificities,
#                    training_sensitivities = training_sensitivities,
#                    testing_auc = testing_auc,
#                    testing_specificities = testing_specificities,
#                    testing_sensitivities = testing_sensitivities)
#
# }
#
# base::saveRDS(results, file = "inst/extdata/GRFClassifierHNSC.rds")
#
results <-
    base::readRDS(file = "inst/extdata/GRFClassifierHNSC.rds")

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
        base::as.character("Training")

    training_results <-
        base::data.frame(auc, partition, set)

    auc <-
        magrittr::extract2(results, i) %>%
        magrittr::use_series(testing_auc)

    partition <-
        magrittr::extract(partitions, i)

    set <-
        base::as.character("Testing")

    testing_results <-
        base::data.frame(auc, partition, set)

    auc_results <-
        base::rbind(auc_results, training_results, testing_results)

}

# dplyr::filter(auc_results, set == "Training") %>%
#     dplyr::group_by(partition) %>%
#     dplyr::summarize(
#         `Min.` = base::min(auc),
#         `1st Qu.` = stats::quantile(auc, 0.25),
#         Median = stats::median(auc),
#         Mean = base::mean(auc),
#         `3rd Qu.` = stats::quantile(auc, 0.75),
#         Max = base::max(auc)
#     ) %>%
#     dplyr::rename(Partition = partition) %>%
#     knitr::kable(digits = 3)
#
dplyr::filter(auc_results, set == "Testing") %>%
    dplyr::group_by(partition) %>%
    dplyr::summarize(
        `Min.` = base::min(auc),
        `1st Qu.` = stats::quantile(auc, 0.25),
        Median = stats::median(auc),
        Mean = base::mean(auc),
        `3rd Qu.` = stats::quantile(auc, 0.75),
        Max = base::max(auc)
    ) %>%
    dplyr::rename(Partition = partition) %>%
    knitr::kable(digits = 3)

# dplyr::filter(auc_results, set == "Training") %>%
#     ggplot2::ggplot(ggplot2::aes(x = partition, y = auc)) +
#     ggplot2::stat_boxplot(geom = "errorbar") +
#     ggplot2::geom_boxplot(notch = TRUE) +
#     ggplot2::theme_linedraw() +
#     ggplot2::labs(
#         title = "Label Propagation using Gaussian Random Fields",
#         subtitle = "Training, Top 5K Features by MAD, 100 Bootstrap Samples",
#         x = NULL, y = "AUC"
#     )
#
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

ggplot2::ggsave("inst/extdata/HNSCGRF.png", plot = ggplot2::last_plot(),
                device = "png", width = 7, height = 7)

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
        base::as.character("Training")

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
        base::as.character("Testing")

    testing_results <-
        base::data.frame(specificity, sensitivity, partition, set) %>%
        dplyr::add_row(specificity = 1, sensitivity = 0, partition = partition,
                       set = set) %>%
        dplyr::arrange(sensitivity)

    roc_results <-
        base::rbind(roc_results, training_results, testing_results)

}

# dplyr::filter(roc_results, set == "Training") %>%
#     ggplot2::ggplot(ggplot2::aes(specificity, sensitivity, color = partition)) +
#     ggplot2::geom_step() +
#     ggplot2::scale_x_reverse() +
#     ggplot2::geom_abline(intercept = 1) +
#     ggplot2::theme_linedraw() +
#     ggplot2::theme(legend.position = "bottom") +
#     ggplot2::labs(
#         title = "Label propagation using Gaussian Random Fields",
#         subtitle = "Training, Top 5K Features by MAD, 100 Bootstrap Samples",
#         x = "Specificity", y = "Sensitivity", color = NULL
#     ) +
#     ggplot2::scale_color_discrete(guide = ggplot2::guide_legend(nrow = 1))
#
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

ggplot2::ggsave("inst/extdata/HNSCROC.png", plot = ggplot2::last_plot(),
                device = "png", width = 7, height = 7)
