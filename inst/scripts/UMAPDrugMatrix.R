library(magrittr)

Carcinogenic <-
    dplyr::select(SSML::miniDrugMatrix, 1)

DrugMatrix <-
    dplyr::select(SSML::miniDrugMatrix, -1) %>%
    umap::umap() %>%
    magrittr::use_series(layout) %>%
    base::as.data.frame() %>%
    dplyr::rename(X = V1) %>%
    dplyr::rename(Y = V2) %>%
    dplyr::bind_cols(Carcinogenic) %>%
    dplyr::filter(!base::is.na(Carcinogenic))

ggplot2::ggplot(DrugMatrix, ggplot2::aes(X, Y, shape = Carcinogenic,
                                         color = Carcinogenic)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(title = "DrugMatrix UMAP",
                  subtitle = "Top 5K Features by MAD", x = NULL, y = NULL)

ggplot2::ggsave("inst/extdata/DrugMatrixUMAP.png", plot = ggplot2::last_plot(),
                device = "png", width = 7, height = 7)
