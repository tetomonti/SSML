library(Biobase)
library(magrittr)

load("~/monti/inst/extdata/CarcinogenomeProject.rda")

phenoData(CarcinogenomeProject) %>%
    slot("data") %>%
    as.data.frame()

featureData(CarcinogenomeProject) %>%
    slot("data") %>%
    as.data.frame()
