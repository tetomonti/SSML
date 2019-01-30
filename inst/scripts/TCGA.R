library(Biobase)
library(e1071)
library(RSSL)
library(pROC)

# Load data
tcga <- readRDS(paste("/Users/JimmyKiely/Documents/Monti Lab/",
                      "BRCA_2018_11_01_subtyped_DESeq2_log_eSet.rds", sep=""))

# Take only pheno data
ph <- phenoData(tcga)
pheno <- ph[ph$pam50=="LumB" | ph$pam50=="LumA",]

# Create a binary variable for pam50 (0 = LumA, 1 = LumB)
pheno$pam50_binary <- as.factor(as.numeric(pheno$pam50=="LumB"))

# Split data into 70% training, 30% test
smp_size <- floor(0.70 * nrow(pheno))
set.seed(123)
train_tcga <- sample(seq_len(nrow(pheno)), size = smp_size)

# Split the phenotype data by test and training
train_pheno <- pheno[train_tcga, ]
test_pheno <- pheno[-train_tcga, ]

# Split expression data by phenotype split
train_exp <- t(exprs(tcga)[,rownames(train_pheno)])
test_exp <- t(exprs(tcga)[,rownames(test_pheno)])

# Supervised SVM with linear kernel ('svm')
model <- svm(x=train_exp, y=train_pheno$pam50_binary, scale=F, kernel="linear")
pred <- predict(model, newdata = test_exp)

# ROC curve for supervised 'svm'
roc <- roc(as.numeric(as.character(test_pheno$pam50_binary)), 
           as.numeric(as.character(pred)))
plot(smooth(roc, method='fitdistr'), main="ROC Plot for LumA vs LumB Supervised 'svm' TCGA")
text(x=0.2, y=0.2, labels=paste("AUC:",auc(roc)))

# Indices for randomly unlabeled data, setting 25% unlabeled
set.seed(123)
max_ind <- length(train_exp[,1]) 
percent_unlabeled <- 0.25
random_ind <- sample(1:max_ind, floor(percent_unlabeled * max_ind))

# Splitting data into labeled and unlabeled
labeled <- train_exp[-random_ind,]
labels <- train_pheno$pam50_binary[-random_ind]
unlabeled <- train_exp[random_ind,]

# Semi-supervised SVM with linear kernel
model_2 <- WellSVM(X=labeled, y=labels, X_u=unlabeled, x_center=F, gamma=0)
pred_2 <- predict(model_2, newdata = test_exp)

# ROC curve for WellSVM
roc_2 <- roc(as.numeric(as.character(test_pheno$pam50_binary)), as.numeric(as.character(pred_2)))
plot(smooth(roc_2, method='fitdistr'), main="ROC Plot for LumA vs LumB Semi-Supervised 'WellSVM' TCGA")
text(x=0.2, y=0.2, labels=paste("AUC:",auc(roc_2)))
