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

provision_data <- function(x, column_name, percent_training, percent_labeled,
                           seed_number) {
    base::set.seed(seed_number)

    size_training <-
        base::nrow(x) %>%
        magrittr::multiply_by(percent_training) %>%
        base::round()

    rows_training <-
        base::nrow(x) %>%
        base::sample(size_training)

    data_training <-
        dplyr::slice(x, rows_training)

    data_testing <-
        dplyr::slice(x, -rows_training)

    size_labeled_training <-
        base::nrow(data_training) %>%
        magrittr::multiply_by(percent_labeled) %>%
        base::round()

    size_labeled_testing <-
        base::nrow(data_testing) %>%
        magrittr::multiply_by(percent_labeled) %>%
        base::round()

    rows_labels_training <-
        base::nrow(data_training) %>%
        base::sample(size_labeled_training)

    rows_labels_testing <-
        base::nrow(data_testing) %>%
        base::sample(size_labeled_testing)

    data_training[-rows_labels_training, column_name] <- NA
    # data_testing[-rows_labels_testing, column_name] <- NA

    list(data_training, data_testing)
}

join_data <- function(x, y) {
    x_one <-
        magrittr::extract2(x, 1)

    x_two <-
        magrittr::extract2(x, 2)

    y_one <-
        magrittr::extract2(y, 1)

    y_two <-
        magrittr::extract2(y, 2)

    data_training <-
        dplyr::bind_rows(x_one, y_one)

    data_testing <-
        dplyr::bind_rows(x_two, y_two)

    list(data_training, data_testing)
}

save_data <- function(x, set_id) {
    data_one <-
        magrittr::extract2(x, 1)

    data_two <-
        magrittr::extract2(x, 2)

    name_one <-
        base::paste0("training_", set_id)

    name_two <-
        base::paste0("testing_", set_id)

    file_one <-
        base::paste0("./data/", name_one, ".rda")

    file_two <-
        base::paste0("./data/", name_two, ".rda")

    base::assign(name_one, data_one)
    base::assign(name_two, data_two)

    base::save(list = name_one, file = file_one, compress = TRUE)
    base::save(list = name_two, file = file_two, compress = TRUE)
}

BRCA <-
    prepare_data(BRCA)

LUM_A <-
    dplyr::filter(BRCA, PAM50.mRNA == "Luminal A") %>%
    provision_data("PAM50.mRNA", 0.70, 1.00, 1001)

LUM_B <-
    dplyr::filter(BRCA, PAM50.mRNA == "Luminal B") %>%
    provision_data("PAM50.mRNA", 0.70, 1.00, 1001)

join_data(LUM_A, LUM_B) %>%
    save_data("1A")

LUM_A <-
    dplyr::filter(BRCA, PAM50.mRNA == "Luminal A") %>%
    provision_data("PAM50.mRNA", 0.70, 0.25, 1001)

LUM_B <-
    dplyr::filter(BRCA, PAM50.mRNA == "Luminal B") %>%
    provision_data("PAM50.mRNA", 0.70, 0.25, 1001)

join_data(LUM_A, LUM_B) %>%
    save_data("1D")
