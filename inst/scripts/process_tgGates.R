##Load in data and clean phenotype information
##Keep only highest dose and time of exposure

library(Biobase)
library(data.table)
tggates <- readRDS("/projectnb/cp2018ssml/datasets/tggates/Rat.Liver.in_vivo.repeat.annotated.RDS")
pheno_tggates <- pData(tggates)
#keep only highest dose batches
pheno_tg_topdose <-pheno_tggates[pheno_tggates$DOSE_LEVEL=='High',]
#only keep numeric information
pheno_tg_topdose$SACRIFICE_PERIOD <- as.numeric(gsub("[^\\d]+", "", pheno_tg_topdose$SACRIFICE_PERIOD, perl=TRUE))
#make datatable
dt_pheno <- as.data.table(pheno_tg_topdose)
#get top sacrifice period
pheno_toptime <- dt_pheno[dt_pheno[, .I[SACRIFICE_PERIOD == max(SACRIFICE_PERIOD)], by=Chemical]$V1]
#match the expression data to the phenotype data
expression <- as.data.frame(t(exprs(tggates)[,as.character(pheno_toptime$BARCODE)]))
#take average of replicates
expression$chem <- pheno_toptime$Chemical
expression <- aggregate(. ~ chem, expression, mean)

#match pheno file to expression file
pheno <- as.data.frame(pheno_toptime)
pheno <- subset(pheno,!duplicated(pheno$Chemical))
pheno <- pheno[order(pheno$Chemical),]
rownames(pheno) <- pheno$Chemical


#fix expression file and add rownames as chemicals
expression <- as.matrix(t(expression[,2:ncol(expression)]))
colnames(expression) <- pheno$Chemical

#make new ExpressionSet
gates <- ExpressionSet(expression, phenoData = AnnotatedDataFrame(pheno))
save(gates, file='/projectnb/cp2018ssml/workdir/tggates_condensed.rda')







