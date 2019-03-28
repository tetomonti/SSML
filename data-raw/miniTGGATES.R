library(magrittr)

load("~/SSML/inst/extdata/TGGATES.rda")

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

TGGATES <-
    prepare_data(TGGATES) %>%
    tibble::column_to_rownames() %>%
    dplyr::select(23, 29:11272) %>%
    dplyr::select_if(not_all_na) %>%
    dplyr::select_if(not_zero_variance)

top_features <-
    dplyr::select(TGGATES, -1) %>%
    base::apply(2, stats::mad) %>%
    base::sort(decreasing = TRUE, index.return = TRUE) %>%
    magrittr::extract2(2) %>%
    magrittr::add(1) %>%
    magrittr::extract(1:5000)

miniTGGATES <-
    dplyr::select(TGGATES, Liver.Cancer, top_features) %>%
    dplyr::rename(Carcinogenic = Liver.Cancer) %>%
    dplyr::mutate(Carcinogenic = forcats::as_factor(Carcinogenic))

usethis::use_data(miniTGGATES)
