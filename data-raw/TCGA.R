TCGA <- curatedTCGAData::curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2GeneNorm", dry.run = FALSE)
colData <- MultiAssayExperiment::colData(TCGA)
sampleMap <- MultiAssayExperiment::sampleMap(TCGA)
colData <- colData[sampleMap$primary, ]
S4Vectors::rownames(colData) <- as.character(sampleMap$colname)
TCGA <- TCGA[[1]]
SummarizedExperiment::assayNames(TCGA) <- as.character("RNASeq2GeneNorm")
SummarizedExperiment::colData(TCGA) <- colData
usethis::use_data(TCGA)
