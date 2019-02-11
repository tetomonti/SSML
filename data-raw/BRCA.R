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

usethis::use_data(BRCA)
