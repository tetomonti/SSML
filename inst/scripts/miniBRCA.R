library(ggplot2)
library(magrittr)

base::load("data/miniBRCA.rda")

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

unlabeled_data <-
  dplyr::slice(training_data, -i_labeled) %>%
  dplyr::select(-Subtype)

unlabeled_labels <-
  dplyr::slice(training_data, -i_labeled) %>%
  magrittr::use_series(Subtype)

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

training_auc <-
  pROC::auc(unlabeled_classes, unlabeled_probabilities) %>%
  base::round(digits = 4) %>%
  base::formatC(digits = 4, format = "f")

training_text <-
  base::paste("AUC =", training_auc)

graphics::plot(training_roc, main = "Luminal A/B using Gaussian Random Fields")
graphics::text(x = 0.2, y = 0.2, training_text)

predicted_labels <-
  RSSL::predict(training_classifier)

predicted_data <-
  dplyr::slice(training_data, -i_labeled) %>%
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
  RSSL::GRFClassifier(X = joined_data, y = joined_labels, X_u = testing_data)

testing_probabilities <-
  methods::slot(testing_classifier, "responsibilities")

testing_seq <-
  base::nrow(testing_probabilities) %>%
  base::seq()

testing_probabilities <-
  magrittr::extract(testing_probabilities, testing_seq, 1)

testing_roc <-
  pROC::roc(testing_classes, testing_probabilities)

testing_auc <-
  pROC::auc(testing_classes, testing_probabilities) %>%
  base::round(digits = 4) %>%
  base::formatC(digits = 4, format = "f")

testing_text <-
  base::paste("AUC =", testing_auc)

graphics::plot(testing_roc, main = "Luminal A/B using Gaussian Random Fields")
graphics::text(x = 0.2, y = 0.2, testing_text)
