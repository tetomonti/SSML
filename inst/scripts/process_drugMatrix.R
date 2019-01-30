

##Load in data and clean phenotype information
##Keep only highest dose and time of exposure


library(Biobase)
library(data.table)

load("/projectnb/cp2018ssml/datasets/drugmatrix/data_curated_mich_ensg.RData")
ph <- pData(data_curated)
##change format of dose and duration to numeric only
ph$DOSE <- as.numeric(gsub("[^\\d]+", "", ph$DOSE, perl=TRUE))
ph$DURATION <- as.numeric(gsub("[^\\d]+", "", ph$DURATION, perl=TRUE))
#turn into datatable 
group <- as.data.table(ph)
group$rownames <- rownames(ph)
#take max of duration
new_ph <- group[group[, .I[DURATION == max(DURATION)], by=chem]$V1]
group2 <- as.data.table(new_ph)
#take max of dose
pheno <- group2[group2[, .I[DOSE == max(DOSE)], by=chem]$V1]
pheno <- pheno[!is.na(pheno$chem),]
#match expression to pheno data
exp <- as.data.frame(t(exprs(data_curated)[,pheno$rownames]))
#take average across chemicals for expression
exp$chem <- pheno$chem
expression <- aggregate(. ~ chem, exp, mean)
#match phenodata to expression data
pheno <- as.data.frame(pheno)
pheno <- pheno[,1:ncol(pheno)-1]
pheno <- subset(pheno,!duplicated(pheno$chem))
pheno <- pheno[order(pheno$chem),]
#fix expression and pheno to add chemical names as identifiers 
expression <- as.matrix(t(expression[,2:ncol(expression)]))
colnames(expression) <- pheno$chem
rownames(pheno) <- pheno$chem
#create expressionset and save
drugmat <- ExpressionSet(expression, phenoData = AnnotatedDataFrame(pheno))
save(drugmat, file='/projectnb/cp2018ssml/workdir/drugmatrix_condensed.rda')
