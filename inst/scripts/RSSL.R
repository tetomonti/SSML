### Loading libraries
library(Biobase)
library(e1071)
library(RSSL)
library(pROC)
library(dplyr)

# Loading in training and testing 1D data
training <- readRDS("/Users/JimmyKiely/Documents/Monti Lab/Datasets/training_1D.rds")
testing <- readRDS("/Users/JimmyKiely/Documents/Monti Lab/Datasets/testing_1D.rds")

# Get labels and convert to binary -- Luminal A= 0 and Luminal B=1 
train_ss_labels <- as.factor(as.numeric(training$PAM50.mRNA=="Luminal B"))

# Get only expression data
train_ss_exp <- training[,2686:ncol(training)]

# Get top 5000 genes by median absolute deviation
colmad <- apply(train_ss_exp, 2, mad)
top5000 <- tail(sort(colmad),5000)

# Subset top 5000 genes in training
train_ss_exp_5000 <- train_ss_exp[,colnames(train_ss_exp) %in% names(top5000) ]

# Subset top 5000 genes in testing
test_exp <- testing[,2686:ncol(testing)]
test_exp_5000 <- test_exp[,colnames(test_exp) %in% names(top5000) ]

# Get labels and convert to binary -- Luminal A= 0 and Luminal B=1 
test_labels <- as.factor(as.numeric(testing$PAM50.mRNA=="Luminal B"))

# Remove all NA samples (fully labeled dataset)
training_fl <- training[!is.na(training$PAM50.mRNA),]

# Get only expression
train_fl_exp <-  training_fl[,2686:ncol(training_fl)]

# Get top 5000 genes
train_fl_exp_5000 <- train_fl_exp[,colnames(train_fl_exp) %in% names(top5000) ]

# Get lables in binary
train_fl_labels <- as.factor(as.numeric(training_fl$PAM50.mRNA=="Luminal B"))

# Supervised SVM with linear kernel ('svm')
model_fl <- svm(x=train_fl_exp_5000, y=train_fl_labels, scale=T)
pred_fl <- predict(model_fl, newdata = test_exp_5000)

# ROC curve
roc_fl <- roc(as.numeric(as.character(test_labels)), 
              as.numeric(as.character(pred_fl)))
plot(smooth(roc_fl, method='fitdistr'), main="ROC Plot for LumA vs LumB Supervised 'svm' TCGA")
text(x=0.2, y=0.2, labels=paste("AUC:",auc(roc_fl)))

# Grabbing unlabeled expression 5000 data
unlabeled_exp_5000 <- train_ss_exp_5000[which(is.na(train_ss_labels)),]

# Semi-supervised SVM with linear kernel
model_ss <- WellSVM(X=train_fl_exp_5000, y=train_fl_labels, X_u=unlabeled_exp_5000, scale=T, gamma=0)
pred_ss <- predict(model_ss, newdata = test_exp_5000)

# ROC curve for WellSVM
roc_ss <- roc(as.numeric(as.character(test_labels)), as.numeric(as.character(pred_ss)))
plot(smooth(roc_ss, method='fitdistr'), main="ROC Plot for LumA vs LumB Semi-Supervised 'WellSVM' TCGA")
text(x=0.2, y=0.2, labels=paste("AUC:",auc(roc_ss)))
