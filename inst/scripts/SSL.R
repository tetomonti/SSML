library(SSML)

data("training_1D")
data("testing_1D")

not_all_na <- function(x) {
    base::is.na(x) %>%
        base::all() %>%
        base::isFALSE()
}

not_zero_variance <- function(x) {
    base::unique(x) %>%
        base::length() %>%
        magrittr::is_greater_than(1)
}

is_not_character <- function(x) {
    base::is.character(x) %>%
        base::isFALSE()
}

is_not_na <- function(x) {
    base::is.na(x) %>%
        base::isFALSE()
}

training_colnames <-
    dplyr::select_if(training_1D, not_all_na) %>%
    dplyr::select_if(not_zero_variance) %>%
    base::colnames()

testing_colnames <-
    dplyr::select_if(testing_1D, not_all_na) %>%
    dplyr::select_if(not_zero_variance) %>%
    base::colnames()

intersecting_colnames <-
    intersect(training_colnames, testing_colnames)

training_1D <-
    training_1D[, intersecting_colnames]

testing_1D <-
    testing_1D[, intersecting_colnames]

training_1D <-
    dplyr::mutate(training_1D, gender = forcats::as_factor(gender)) %>%
    dplyr::mutate(PAM50.mRNA = forcats::as_factor(PAM50.mRNA)) %>%
    dplyr::select_if(is_not_character)

testing_1D <-
    dplyr::mutate(testing_1D, gender = forcats::as_factor(gender)) %>%
    dplyr::mutate(PAM50.mRNA = forcats::as_factor(PAM50.mRNA)) %>%
    dplyr::select_if(is_not_character)

xl <-
    dplyr::filter(training_1D, !base::is.na(PAM50.mRNA)) %>%
    dplyr::select(-PAM50.mRNA) %>%
    dplyr::mutate(gender = as.integer(gender))

xu <-
    dplyr::filter(training_1D, base::is.na(PAM50.mRNA)) %>%
    dplyr::select(-PAM50.mRNA) %>%
    dplyr::mutate(gender = base::as.integer(gender))

yl <-
    dplyr::filter(training_1D, !base::is.na(PAM50.mRNA)) %>%
    magrittr::use_series(PAM50.mRNA) %>%
    base::as.integer()

ll <-
    SSL::sslGmmEM(xl[1:72, 735:738], yl[1:72], xu[1:216, 735:738])

labeled_data <- xl[1:72, 735:744]
unlabeled_data <- xu[1:216, 735:744]
labels_vector <- yl[1:72]

save(labeled_data, file = "~/labeled_data.rda", compress = TRUE)
save(unlabeled_data, file = "~/unlabeled_data.rda", compress = TRUE)
save(labels_vector, file = "~/labels_vector.rda", compress = TRUE)
