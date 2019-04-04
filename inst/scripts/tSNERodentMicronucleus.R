library(magrittr)

Genotoxic <-
    dplyr::select(SSML::miniRodentMicronucleus, 1)

RodentMicronucleus <-
    dplyr::select(SSML::miniRodentMicronucleus, -1) %>%
    Rtsne::Rtsne() %>%
    magrittr::use_series(Y) %>%
    base::as.data.frame() %>%
    dplyr::rename(X = V1) %>%
    dplyr::rename(Y = V2) %>%
    dplyr::bind_cols(Genotoxic) %>%
    dplyr::filter(!base::is.na(Genotoxic))

ggplot2::ggplot(RodentMicronucleus, ggplot2::aes(X, Y, shape = Genotoxic,
                                     color = Genotoxic)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(title = "RodentMicronucleus tSNE",
                  subtitle = "Top 5K Features by MAD", x = NULL, y = NULL)

ggplot2::ggsave("inst/extdata/RodentMicronucleustSNE.png", plot = ggplot2::last_plot(),
                device = "png", width = 7, height = 7)
