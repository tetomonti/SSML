if (!require("magrittr", character.only = TRUE)) {
    BiocManager::install("magrittr")
    require("magrittr", character.only = TRUE)
}

BRCA <-
    curatedTCGAData::curatedTCGAData(diseaseCode = "BRCA",
                                     assays = "RNASeq2GeneNorm",
                                     dry.run = FALSE)

colData <-
    MultiAssayExperiment::colData(BRCA)

sampleMap <-
    MultiAssayExperiment::sampleMap(BRCA)

primary <-
    use_series(sampleMap, primary)

colname <-
    use_series(sampleMap, colname)

colData <-
    extract(colData, primary, NULL)

S4Vectors::rownames(colData) <-
    as.character(sampleMap$colname)

BRCA <-
    BRCA[[1]]

SummarizedExperiment::assayNames(BRCA) <-
    as.character("RNASeq2GeneNorm")

SummarizedExperiment::colData(BRCA) <-
    colData


