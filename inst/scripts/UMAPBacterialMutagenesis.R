library(magrittr)

Genotoxic <-
    dplyr::select(SSML::miniBacterialMutagenesis, 1)

BacterialMutagenesis <-
    dplyr::select(SSML::miniBacterialMutagenesis, -1) %>%
    umap::umap() %>%
    magrittr::use_series(layout) %>%
    base::as.data.frame() %>%
    dplyr::rename(X = V1) %>%
    dplyr::rename(Y = V2) %>%
    dplyr::bind_cols(Genotoxic) %>%
    dplyr::filter(!base::is.na(Genotoxic))

ggplot2::ggplot(BacterialMutagenesis, ggplot2::aes(X, Y, shape = Genotoxic, color = Genotoxic)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(title = "BacterialMutagenesis UMAP", subtitle = "Top 5K Features by MAD",
                  x = NULL, y = NULL)

ggplot2::ggsave("inst/extdata/BacterialMutagenesisUMAP.png", plot = ggplot2::last_plot(),
                device = "png", width = 7, height = 7)
