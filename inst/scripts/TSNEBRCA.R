library(magrittr)

Subtype <-
    dplyr::select(SSML::miniBRCA, 1)

BRCA <-
    dplyr::select(SSML::miniBRCA, -1) %>%
    Rtsne::Rtsne() %>%
    magrittr::use_series(Y) %>%
    base::as.data.frame() %>%
    dplyr::rename(X = V1) %>%
    dplyr::rename(Y = V2) %>%
    dplyr::bind_cols(Subtype)

ggplot2::ggplot(BRCA, ggplot2::aes(X, Y, shape = Subtype, color = Subtype)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(title = "BRCA TSNE", subtitle = "Top 5K Features by MAD",
                  x = NULL, y = NULL)
