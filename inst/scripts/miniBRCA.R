library(magrittr)

base::load("data/miniBRCA.rda")

n_unlabeled <-
    c(29, 58, 87, 116, 145, 174, 202)

n_bootstrap <-
    seq(100)

results <-
    list()

for (i in n_unlabeled) {
    j <-
        match(i, n_unlabeled)

    training_auc <-
        numeric()

    training_specificities <-
        list()

    training_sensitivities <-
        list()

    testing_auc <-
        numeric()

    testing_specificities <-
        list()

    testing_sensitivities <-
        list()

    for (k in n_bootstrap) {
        n_train <-
            base::nrow(miniBRCA) %>%
            magrittr::multiply_by(0.70) %>%
            base::round()

        i_train <-
            base::nrow(miniBRCA) %>%
            base::sample(n_train)

        training_data <-
            dplyr::slice(miniBRCA, i_train)

        testing_data <-
            dplyr::slice(miniBRCA, -i_train)

        n_labeled <-
            base::nrow(training_data) %>%
            magrittr::multiply_by(0.30) %>%
            base::round()

        i_labeled <-
            base::nrow(training_data) %>%
            base::sample(n_labeled)

        labeled_data <-
            dplyr::slice(training_data, i_labeled) %>%
            dplyr::select(-Subtype)

        labeled_labels <-
            dplyr::slice(training_data, i_labeled) %>%
            magrittr::use_series(Subtype)

        i_unlabeled <-
            dplyr::slice(training_data, -i_labeled) %>%
            base::nrow() %>%
            base::seq() %>%
            base::sample(i)

        unlabeled_data <-
            dplyr::slice(training_data, -i_labeled) %>%
            dplyr::slice(i_unlabeled) %>%
            dplyr::select(-Subtype)

        unlabeled_labels <-
            dplyr::slice(training_data, -i_labeled) %>%
            magrittr::use_series(Subtype) %>%
            magrittr::extract(i_unlabeled)

        unlabeled_classes <-
            base::as.integer(unlabeled_labels)

        training_classifier <-
            RSSL::GRFClassifier(X = labeled_data, y = labeled_labels,
                                X_u = unlabeled_data)

        training_probabilities <-
            methods::slot(training_classifier, "responsibilities")

        training_seq <-
            base::nrow(training_probabilities) %>%
            base::seq()

        unlabeled_probabilities <-
            magrittr::extract(training_probabilities, training_seq, 1)

        training_roc <-
            pROC::roc(unlabeled_classes, unlabeled_probabilities)

        training_auc[k] <-
            pROC::auc(unlabeled_classes, unlabeled_probabilities) %>%
            base::as.numeric()

        training_specificities[[k]] <-
            magrittr::use_series(training_roc, specificities)

        training_sensitivities[[k]] <-
            magrittr::use_series(training_roc, sensitivities)

        predicted_labels <-
            RSSL::predict(training_classifier)

        predicted_data <-
            dplyr::slice(training_data, -i_labeled) %>%
            dplyr::slice(i_unlabeled) %>%
            dplyr::mutate(Subtype = predicted_labels)

        joined_data <-
            dplyr::slice(training_data, i_labeled) %>%
            dplyr::bind_rows(predicted_data) %>%
            dplyr::select(-Subtype)

        joined_labels <-
            dplyr::slice(training_data, i_labeled) %>%
            dplyr::bind_rows(predicted_data) %>%
            magrittr::use_series(Subtype)

        testing_data <-
            dplyr::slice(miniBRCA, -i_train) %>%
            dplyr::select(-Subtype)

        testing_labels <-
            dplyr::slice(miniBRCA, -i_train) %>%
            magrittr::use_series(Subtype)

        testing_classes <-
            base::as.integer(testing_labels)

        testing_classifier <-
            RSSL::GRFClassifier(X = joined_data, y = joined_labels,
                                X_u = testing_data)

        testing_probabilities <-
            methods::slot(testing_classifier, "responsibilities")

        testing_seq <-
            base::nrow(testing_probabilities) %>%
            base::seq()

        testing_probabilities <-
            magrittr::extract(testing_probabilities, testing_seq, 1)

        testing_roc <-
            pROC::roc(testing_classes, testing_probabilities)

        testing_auc[k] <-
            pROC::auc(testing_classes, testing_probabilities) %>%
            base::as.numeric()

        testing_specificities[[k]] <-
            magrittr::use_series(testing_roc, specificities)

        testing_sensitivities[[k]] <-
            magrittr::use_series(testing_roc, sensitivities)
    }

    training_specificities <-
        purrr::reduce(training_specificities, base::cbind) %>%
        base::apply(1, base::mean)

    training_sensitivities <-
        purrr::reduce(training_sensitivities, base::cbind) %>%
        base::apply(1, base::mean)

    testing_specificities <-
        purrr::reduce(testing_specificities, base::cbind) %>%
        base::apply(1, base::mean)

    testing_sensitivities <-
        purrr::reduce(testing_sensitivities, base::cbind) %>%
        base::apply(1, base::mean)

    results[[j]] <-
        base::list(training_auc = training_auc,
                   training_specificities = training_specificities,
                   training_sensitivities = training_sensitivities,
                   testing_auc = testing_auc,
                   testing_specificities = testing_specificities,
                   testing_sensitivities = testing_sensitivities)

}
