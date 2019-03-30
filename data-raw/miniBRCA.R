library(magrittr)

BRCA <-
  curatedTCGAData::curatedTCGAData(diseaseCode = "BRCA",
                                   assays = "RNASeq2GeneNorm",
                                   dry.run = FALSE)

colData <-
  MultiAssayExperiment::colData(BRCA)

sampleMap <-
  MultiAssayExperiment::sampleMap(BRCA)

colData <-
  colData[sampleMap$primary, ]

S4Vectors::rownames(colData) <-
  base::as.character(sampleMap$colname)

BRCA <-
  BRCA[[1]]

SummarizedExperiment::assayNames(BRCA) <-
  base::as.character("RNASeq2GeneNorm")

SummarizedExperiment::colData(BRCA) <-
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

BRCA <-
  prepare_data(BRCA) %>%
  dplyr::filter(!base::is.na(PAM50.mRNA)) %>%
  dplyr::filter(PAM50.mRNA != "Basal-like") %>%
  dplyr::filter(PAM50.mRNA != "HER2-enriched") %>%
  dplyr::filter(PAM50.mRNA != "Normal-like") %>%
  dplyr::filter(stringr::str_detect(rowname, "11A", negate = TRUE)) %>%
  dplyr::filter(stringr::str_detect(rowname, "11B", negate = TRUE)) %>%
  tibble::column_to_rownames() %>%
  dplyr::select(2674, 2685:23185) %>%
  dplyr::select_if(not_all_na) %>%
  dplyr::select_if(not_zero_variance)

top_features <-
    dplyr::select(BRCA, -1) %>%
    base::apply(2, stats::mad) %>%
    base::sort(decreasing = TRUE, index.return = TRUE) %>%
    magrittr::extract2(2) %>%
    magrittr::add(1) %>%
    magrittr::extract(1:5000)

miniBRCA <-
    dplyr::select(BRCA, PAM50.mRNA, top_features) %>%
    dplyr::rename(Subtype = PAM50.mRNA) %>%
    dplyr::mutate(Subtype = forcats::as_factor(Subtype))

usethis::use_data(miniBRCA)
