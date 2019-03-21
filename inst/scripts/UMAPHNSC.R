library(magrittr)

Grade <-
    dplyr::select(SSML::miniHNSC, 1)

HNSC <-
    dplyr::select(SSML::miniHNSC, -1) %>%
    umap::umap() %>%
    magrittr::use_series(layout) %>%
    base::as.data.frame() %>%
    dplyr::rename(X = V1) %>%
    dplyr::rename(Y = V2) %>%
    dplyr::bind_cols(Grade)

ggplot2::ggplot(HNSC, ggplot2::aes(X, Y, shape = Grade, color = Grade)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(title = "HNSC UMAP", subtitle = "Top 5K Features by MAD",
                  x = NULL, y = NULL)
