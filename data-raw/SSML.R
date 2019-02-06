library(SSML)

extract_colData <- function(x) {
    SummarizedExperiment::colData(x) %>%
        S4Vectors::as.data.frame() %>%
        tibble::rownames_to_column()
}

extract_assay <- function(x) {
    SummarizedExperiment::assay(x) %>%
        base::t() %>%
        S4Vectors::as.data.frame() %>%
        tibble::rownames_to_column()
}

prepare_data <- function(x) {
    base::list(extract_colData, extract_assay) %>%
        purrr::invoke_map(x = x) %>%
        purrr::reduce(base::merge.data.frame)
}

TCGA <-
    prepare_data(TCGA)

LUM_A <-
    dplyr::filter(TCGA, PAM50.mRNA == "Luminal A")

LUM_B <-
    dplyr::filter(TCGA, PAM50.mRNA == "Luminal B")

nrow(LUM_A)

percent_training = 0.75
percent_labeled = 100
seed_number = 1001

base::set.seed(seed_number)

size_training <-
    base::nrow(LUM_A) %>%
    magrittr::multiply_by(percent_training) %>%
    base::round()

rows_training <-
    base::nrow(LUM_A) %>%
    base::sample(size_training)

data_training <-
    dplyr::slice(rows_training)

data_testing <-
    dplyr::slice(-rows_training)

