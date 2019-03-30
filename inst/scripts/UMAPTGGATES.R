library(magrittr)

Carcinogenic <-
    dplyr::select(SSML::miniTGGATES, 1)

TGGATES <-
    dplyr::select(SSML::miniTGGATES, -1) %>%
    umap::umap() %>%
    magrittr::use_series(layout) %>%
    base::as.data.frame() %>%
    dplyr::rename(X = V1) %>%
    dplyr::rename(Y = V2) %>%
    dplyr::bind_cols(Carcinogenic) %>%
    dplyr::filter(!base::is.na(Carcinogenic))

ggplot2::ggplot(TGGATES, ggplot2::aes(X, Y, shape = Carcinogenic,
                                         color = Carcinogenic)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(title = "TG-GATES UMAP",
                  subtitle = "Top 5K Features by MAD", x = NULL, y = NULL)

ggplot2::ggsave("inst/extdata/TGGATESUMAP.png", plot = ggplot2::last_plot(),
                device = "png", width = 5, height = 5)
