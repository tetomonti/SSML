library(SSML)
library(SSL)

data("training_1D")
data("testing_1D")

not_all_na <- function(x) {
    is.na(x) %>%
        all() %>%
        isFALSE()
}

not_zero_variance <- function(x) {
    unique(x) %>%
        length() %>%
        magrittr::is_greater_than(1)
}

is_not_character <- function(x) {
    is.character(x) %>%
        isFALSE()
}

is_not_na <- function(x) {

    is.na(x) %>%
        isFALSE()
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

training_1D <- training_1D[, intersecting_colnames]
testing_1D <- testing_1D[, intersecting_colnames]

training_1D <-
    dplyr::mutate(training_1D, gender = forcats::as_factor(gender)) %>%
    dplyr::mutate(PAM50.mRNA = forcats::as_factor(PAM50.mRNA)) %>%
    dplyr::select_if(is_not_character)

testing_1D <-
    dplyr::mutate(testing_1D, gender = forcats::as_factor(gender)) %>%
    dplyr::mutate(PAM50.mRNA = forcats::as_factor(PAM50.mRNA)) %>%
    dplyr::select_if(is_not_character)

xl <-
    dplyr::filter(training_1D, !is.na(PAM50.mRNA)) %>%
    dplyr::select(-PAM50.mRNA) %>%
    dplyr::mutate(gender = as.integer(gender)) %>%
    dplyr::select(735:739)

xu <-
    dplyr::filter(training_1D, is.na(PAM50.mRNA)) %>%
    dplyr::select(-PAM50.mRNA) %>%
    dplyr::mutate(gender = as.integer(gender)) %>%
    dplyr::select(735:739)

yl <-
    dplyr::filter(training_1D, !is.na(PAM50.mRNA)) %>%
    magrittr::use_series(PAM50.mRNA) %>%
    base::as.integer()

xl <- xl[1:70,1:4]
xu <- xu[1:91,1:4]
yl <- yl[1:70]
l<-sslGmmEM(xl,yl,xu)
