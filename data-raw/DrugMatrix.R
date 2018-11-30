load(x)

ph <- pData(data_curated)
ph$DOSE <- as.numeric(gsub("[^\\d]+", "", ph$DOSE, perl=TRUE))
ph$DURATION <- as.numeric(gsub("[^\\d]+", "", ph$DURATION, perl=TRUE))
group <- as.data.table(ph)
group$rownames <- rownames(ph)
new_ph <- group[group[, .I[DURATION == max(DURATION)], by=chem]$V1]
group2 <- as.data.table(new_ph)
pheno <- group2[group2[, .I[DOSE == max(DOSE)], by=chem]$V1]
pheno <- pheno[!is.na(pheno$chem),]
exp <- as.data.frame(t(exprs(data_curated)[,pheno$rownames]))
exp$chem <- pheno$chem
expression <- aggregate(. ~ chem, exp, mean)
pheno <- as.data.frame(pheno)
pheno <- pheno[,1:ncol(pheno)-1]
pheno <- subset(pheno,!duplicated(pheno$chem))
expression <- as.matrix(t(expression[,2:ncol(expression)]))
colnames(expression) <- pheno$chem
rownames(pheno) <- pheno$chem
drugmat <- ExpressionSet(expression, phenoData = AnnotatedDataFrame(pheno))
save(drugmat, file='/projectnb/cp2018ssml/workdir/drugmatrix_condensed.rda')
