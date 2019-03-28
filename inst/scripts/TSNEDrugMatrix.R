library(magrittr)

Carcinogenicity <-
    dplyr::select(SSML::miniDrugMatrix, 1)

DrugMatrix <-
    dplyr::select(SSML::miniDrugMatrix, -1) %>%
    Rtsne::Rtsne() %>%
    magrittr::use_series(Y) %>%
    base::as.data.frame() %>%
    dplyr::rename(X = V1) %>%
    dplyr::rename(Y = V2) %>%
    dplyr::bind_cols(Carcinogenicity) %>%
    dplyr::rename(Carcinogenicity = CARCINOGENICITY) %>%
    dplyr::filter(!base::is.na(Carcinogenicity))

ggplot2::ggplot(DrugMatrix, ggplot2::aes(X, Y, shape = Carcinogenicity,
                                         color = Carcinogenicity)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(title = "DrugMatrix TSNE",
                  subtitle = "Top 5K Features by MAD", x = NULL, y = NULL)
