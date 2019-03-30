library(magrittr)

HNSC <-
    curatedTCGAData::curatedTCGAData(diseaseCode = "HNSC",
                                     assays = "RNASeq2GeneNorm",
                                     dry.run = FALSE)

colData <-
    MultiAssayExperiment::colData(HNSC)

sampleMap <-
    MultiAssayExperiment::sampleMap(HNSC)

colData <-
    colData[sampleMap$primary, ]

S4Vectors::rownames(colData) <-
    base::as.character(sampleMap$colname)

HNSC <-
    HNSC[[1]]

SummarizedExperiment::assayNames(HNSC) <-
    base::as.character("RNASeq2GeneNorm")

SummarizedExperiment::colData(HNSC) <-
    colData

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

HNSC <-
    prepare_data(HNSC) %>%
    dplyr::filter(!base::is.na(pathologic_stage)) %>%
    dplyr::mutate(pathologic_stage = dplyr::case_when(
        pathologic_stage == "stage i" ~ "Low",
        pathologic_stage == "stage ii" ~ "Low",
        pathologic_stage == "stage iii" ~ "High",
        pathologic_stage == "stage iva" ~ "High",
        pathologic_stage == "stage ivb" ~ "High",
        pathologic_stage == "stage ivc" ~ "High"
    )) %>%
    dplyr::filter(stringr::str_detect(rowname, "11A", negate = TRUE)) %>%
    tibble::column_to_rownames() %>%
    dplyr::select(7, 1445:21946) %>%
    dplyr::select_if(not_all_na) %>%
    dplyr::select_if(not_zero_variance)

top_features <-
    dplyr::select(HNSC, -1) %>%
    base::apply(2, stats::mad) %>%
    base::sort(decreasing = TRUE, index.return = TRUE) %>%
    magrittr::extract2(2) %>%
    magrittr::add(1) %>%
    magrittr::extract(1:5000)

miniHNSC <-
    dplyr::select(HNSC, pathologic_stage, top_features) %>%
    dplyr::rename(Grade = pathologic_stage) %>%
    dplyr::mutate(Grade = forcats::as_factor(Grade))

usethis::use_data(miniHNSC)
