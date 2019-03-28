library(magrittr)

Grade <-
    dplyr::select(SSML::miniHNSC, 1)

HNSC <-
    dplyr::select(SSML::miniHNSC, -1) %>%
    Rtsne::Rtsne() %>%
    magrittr::use_series(Y) %>%
    base::as.data.frame() %>%
    dplyr::rename(X = V1) %>%
    dplyr::rename(Y = V2) %>%
    dplyr::bind_cols(Grade)

ggplot2::ggplot(HNSC, ggplot2::aes(X, Y, shape = Grade, color = Grade)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(title = "HNSC tSNE", subtitle = "Top 5K Features by MAD",
                  x = NULL, y = NULL)

ggplot2::ggsave("inst/extdata/HNSCtSNE.png", plot = ggplot2::last_plot(),
                device = "png", width = 5, height = 5)
