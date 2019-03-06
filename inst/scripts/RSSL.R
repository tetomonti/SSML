library(tidyverse)
library(SummarizedExperiment)
library(RSSL)
library(pROC)
args = commandArgs(trailingOnly=TRUE)

# Loading in data
load("/projectnb/cp2018ssml/datasets/tcga/BRCA.rda")

# Functions for splitting training/testing, and labeled/unlabeled
extract_colData <- function(x) {
    SummarizedExperiment::colData(x) %>%
        S4Vectors::as.data.frame() %>%
        tibble::rownames_to_column()
}

extract_assay <- function(x) {
    SummarizedExperiment::assay(x) %>%
        base::t() %>%
        S4Vectors::as.data.frame() %>%
        tibble::rownames_to_column()
}

prepare_data <- function(x) {
    base::list(extract_colData, extract_assay) %>%
        purrr::invoke_map(x = x) %>%
        purrr::reduce(base::merge.data.frame)
}

provision_data <- function(x, column_name, percent_training, percent_labeled,
                           seed_number) {
    base::set.seed(seed_number)
    
    size_training <-
        base::nrow(x) %>%
        magrittr::multiply_by(percent_training) %>%
        base::round()
    
    rows_training <-
        base::nrow(x) %>%
        base::sample(size_training)
    
    data_training <-
        dplyr::slice(x, rows_training)
    
    data_testing <-
        dplyr::slice(x, -rows_training)
    
    size_labeled_training <-
        base::nrow(data_training) %>%
        magrittr::multiply_by(percent_labeled) %>%
        base::round()
    
    size_labeled_testing <-
        base::nrow(data_testing) %>%
        magrittr::multiply_by(percent_labeled) %>%
        base::round()
    
    rows_labels_training <-
        base::nrow(data_training) %>%
        base::sample(size_labeled_training)
    
    rows_labels_testing <-
        base::nrow(data_testing) %>%
        base::sample(size_labeled_testing)
    
    data_training[-rows_labels_training, column_name] <- NA
    #data_testing[-rows_labels_testing, column_name] <- NA
    
    base::list(data_training, data_testing)
}

join_data <- function(x, y) {
    x_one <-
        magrittr::extract2(x, 1)
    
    x_two <-
        magrittr::extract2(x, 2)
    
    y_one <-
        magrittr::extract2(y, 1)
    
    y_two <-
        magrittr::extract2(y, 2)
    
    data_training <-
        dplyr::bind_rows(x_one, y_one)
    
    data_testing <-
        dplyr::bind_rows(x_two, y_two)
    
    base::list(data_training, data_testing)
}

# Format data
BRCA <- prepare_data(BRCA)
random_seed <- as.integer(args[1])

LUM_A <-
    dplyr::filter(BRCA, PAM50.mRNA == "Luminal A") %>%
    provision_data("PAM50.mRNA", 0.70, 0.30, random_seed)

LUM_B <-
    dplyr::filter(BRCA, PAM50.mRNA == "Luminal B") %>%
    provision_data("PAM50.mRNA", 0.70, 0.30, random_seed)

x <- join_data(LUM_A, LUM_B)

training <- magrittr::extract2(x, 1)

testing <- magrittr::extract2(x, 2)

#get top 5000 genes by MAD in training
training_exp <- training[,2686:ncol(training)]
colmad_training <- apply(training_exp, 2, mad)
top5000_training <- tail(sort(colmad_training),5000)

##get test data
#only expression
test_exp <- testing[,2686:ncol(testing)]
## get labels and convert to binary -- Luminal A= 0 and Luminal B=1 
test_labels <- as.factor(as.numeric(testing$PAM50.mRNA=="Luminal B"))
#subset top 5000 genes in testing
test_exp_5000 <- test_exp[,colnames(test_exp) %in% names(top5000_training)]

#get fully data with labels and unlabeled data
labeled_training <- training[!is.na(training$PAM50.mRNA),]
unlabeled_training <- training[is.na(training$PAM50.mRNA),]

#split unlabeled data into 7 groups
ul_training_splits <- split(unlabeled_training, cut(1:nrow(unlabeled_training), 7))

#get data for the fully labeled set 
train_fl_exp <-  labeled_training[,2686:ncol(labeled_training)]
#get top 5000 genes
train_fl_exp_5000 <- train_fl_exp[,colnames(train_fl_exp) %in% names(top5000_training) ]
#get lables in binary
train_fl_labels <- as.factor(as.numeric(labeled_training$PAM50.mRNA=="Luminal B"))

# top 5000 edited training, labeled and unlabeled
training_final <- labeled_training[,2686:ncol(labeled_training)]
training_final <- training_final[,colnames(training_final) %in% names(top5000_training)]
training_labels <- as.factor(as.numeric(labeled_training$PAM50.mRNA=="Luminal B"))


# svm for fl + no unlabeled
model_fl <- SVM(X=training_final, y=training_labels, scale=T)
pred_fl <- predict(model_fl, newdata = test_exp_5000)
d_values <- decisionvalues(model_fl, test_exp_5000)
roc_fl <- roc(test_labels, d_values)
saveRDS(model_fl, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/", args[1],"_ss0_training_model.rds", sep=""))
saveRDS(pred_fl, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/",args[1],"_ss0_prediction.rds", sep=""))

pdf(paste("/projectnb/cp2018ssml/workdir/jimmy/roc/",args[1],"_fl_ROC.pdf", sep=""))
plot(roc_fl, main=paste("ROC Plot for LumA vs LumB in semi-supervised FL", sep=""))
text(x=0.2, y=0.2, labels=paste("AUC:", auc(roc_fl)))
dev.off()
auc_vec <- auc(roc_fl)

# ss svm for the fl + 1/7 of the data 
unlabeled_ss1 <- ul_training_splits[[1]][,2686:ncol(ul_training_splits[[1]])]
ss1 <- unlabeled_ss1[,colnames(unlabeled_ss1) %in% names(top5000_training)]
model_ss1 <- WellSVM(X=training_final, y=training_labels, X_u=ss1, scale=T, gamma=0)
pred_ss1 <- predict(model_ss1, newdata = test_exp_5000, probs=TRUE)
d_values_1 <- decisionvalues(model_ss1, test_exp_5000)
roc_1 <- roc(test_labels, d_values_1)
saveRDS(model_ss1, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/", args[1],"_ss1_training_model.rds", sep=""))
saveRDS(pred_ss1, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/",args[1],"_ss1_prediction.rds", sep=""))

pdf(paste("/projectnb/cp2018ssml/workdir/jimmy/roc/",args[1],"_ss1_ROC.pdf", sep=""))
plot(roc_1, main=paste("ROC Plot for LumA vs LumB in semi-supervised ss1", sep=""))
text(x=0.2, y=0.2, labels=paste("AUC:", auc(roc_1)))
dev.off()
auc_vec <- c(auc_vec, auc(roc_1))

# ss svm for the fl + 2/7 of the data 
unlabeled_ss2 <- rbind(unlabeled_ss1, ul_training_splits[[2]][,2686:ncol(ul_training_splits[[2]])])
ss2 <- unlabeled_ss2[,colnames(unlabeled_ss2) %in% names(top5000_training)]
model_ss2 <- WellSVM(X=training_final, y=training_labels, X_u=ss2, scale=T, gamma=0)
pred_ss2 <- predict(model_ss2, newdata = test_exp_5000)
d_values_2 <- decisionvalues(model_ss2, test_exp_5000)
roc_2 <- roc(test_labels, d_values_2)
saveRDS(model_ss2, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/", args[1],"_ss2_training_model.rds", sep=""))
saveRDS(pred_ss2, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/",args[1],"_ss2_prediction.rds", sep=""))

pdf(paste("/projectnb/cp2018ssml/workdir/jimmy/roc/",args[1],"_ss2_ROC.pdf", sep=""))
plot(roc_2, main=paste("ROC Plot for LumA vs LumB in semi-supervised ss2", sep=""))
text(x=0.2, y=0.2, labels=paste("AUC:", auc(roc_2)))
dev.off()
auc_vec <- c(auc_vec, auc(roc_2))

# ss svm for the fl + 3/7 of the data 
unlabeled_ss3 <- rbind(unlabeled_ss2, ul_training_splits[[3]][,2686:ncol(ul_training_splits[[3]])])
ss3 <- unlabeled_ss3[,colnames(unlabeled_ss3) %in% names(top5000_training)]
model_ss3 <- WellSVM(X=training_final, y=training_labels, X_u=ss3, scale=T, gamma=0)
pred_ss3 <- predict(model_ss3, newdata = test_exp_5000)
d_values_3 <- decisionvalues(model_ss3, test_exp_5000)
roc_3 <- roc(test_labels, d_values_3)
saveRDS(model_ss3, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/", args[1],"_ss3_training_model.rds", sep=""))
saveRDS(pred_ss3, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/",args[1],"_ss3_prediction.rds", sep=""))

pdf(paste("/projectnb/cp2018ssml/workdir/jimmy/roc/",args[1],"_ss3_ROC.pdf", sep=""))
plot(roc_3, main=paste("ROC Plot for LumA vs LumB in semi-supervised ss3", sep=""))
text(x=0.2, y=0.2, labels=paste("AUC:", auc(roc_3)))
dev.off()
auc_vec <- c(auc_vec, auc(roc_3))

# ss svm for the fl + 4/7 of the data 
unlabeled_ss4 <- rbind(unlabeled_ss3, ul_training_splits[[4]][,2686:ncol(ul_training_splits[[4]])])
ss4 <- unlabeled_ss4[,colnames(unlabeled_ss4) %in% names(top5000_training)]
model_ss4 <- WellSVM(X=training_final, y=training_labels, X_u=ss4, scale=T, gamma=0)
pred_ss4 <- predict(model_ss4, newdata = test_exp_5000)
d_values_4<- decisionvalues(model_ss4, test_exp_5000)
roc_4 <- roc(test_labels, d_values_4)
saveRDS(model_ss4, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/", args[1],"_ss4_training_model.rds", sep=""))
saveRDS(pred_ss4, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/",args[1],"_ss4_prediction.rds", sep=""))

pdf(paste("/projectnb/cp2018ssml/workdir/jimmy/roc/",args[1],"_ss4_ROC.pdf", sep=""))
plot(roc_4, main=paste("ROC Plot for LumA vs LumB in semi-supervised ss4", sep=""))
text(x=0.2, y=0.2, labels=paste("AUC:", auc(roc_4)))
dev.off()
auc_vec <- c(auc_vec, auc(roc_4))

# ss svm for the fl + 5/7 of the data 
unlabeled_ss5 <- rbind(unlabeled_ss4, ul_training_splits[[5]][,2686:ncol(ul_training_splits[[5]])])
ss5 <- unlabeled_ss5[,colnames(unlabeled_ss5) %in% names(top5000_training)]
model_ss5 <- WellSVM(X=training_final, y=training_labels, X_u=ss5, scale=T, gamma=0)
pred_ss5 <- predict(model_ss5, newdata = test_exp_5000)
d_values_5 <- decisionvalues(model_ss5, test_exp_5000)
roc_5 <- roc(test_labels, d_values_5)
saveRDS(model_ss5, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/", args[1],"_ss5_training_model.rds", sep=""))
saveRDS(pred_ss5, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/",args[1],"_ss5_prediction.rds", sep=""))

pdf(paste("/projectnb/cp2018ssml/workdir/jimmy/roc/",args[1],"_ss5_ROC.pdf", sep=""))
plot(roc_5, main=paste("ROC Plot for LumA vs LumB in semi-supervised ss5", sep=""))
text(x=0.2, y=0.2, labels=paste("AUC:", auc(roc_5)))
dev.off()
auc_vec <- c(auc_vec, auc(roc_5))

# ss svm for the fl + 6/7 of the data 
unlabeled_ss6 <- rbind(unlabeled_ss5, ul_training_splits[[6]][,2686:ncol(ul_training_splits[[6]])])
ss6 <- unlabeled_ss6[,colnames(unlabeled_ss6) %in% names(top5000_training)]
model_ss6 <- WellSVM(X=training_final, y=training_labels, X_u=ss6, scale=T, gamma=0)
pred_ss6 <- predict(model_ss6, newdata = test_exp_5000)
d_values_6 <- decisionvalues(model_ss6, test_exp_5000)
roc_6 <- roc(test_labels, d_values_6)
saveRDS(model_ss6, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/", args[1],"_ss6_training_model.rds", sep=""))
saveRDS(pred_ss6, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/",args[1],"_ss6_prediction.rds", sep=""))

pdf(paste("/projectnb/cp2018ssml/workdir/jimmy/roc/",args[1],"_ss6_ROC.pdf", sep=""))
plot(roc_6, main=paste("ROC Plot for LumA vs LumB in semi-supervised ss6", sep=""))
text(x=0.2, y=0.2, labels=paste("AUC:", auc(roc_6)))
dev.off()
auc_vec <- c(auc_vec, auc(roc_6))

# ss svm for the fl + 7/7 of the data 
unlabeled_ss7 <- rbind(unlabeled_ss6, ul_training_splits[[7]][,2686:ncol(ul_training_splits[[7]])])
ss7 <- unlabeled_ss7[,colnames(unlabeled_ss7) %in% names(top5000_training)]
model_ss7 <- WellSVM(X=training_final, y=training_labels, X_u=ss7, scale=T, gamma=0)
pred_ss7 <- predict(model_ss7, newdata = test_exp_5000)
d_values_7 <- decisionvalues(model_ss7, test_exp_5000)
roc_7 <- roc(test_labels, d_values_7)
saveRDS(model_ss7, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/", args[1],"_ss7_training_model.rds", sep=""))
saveRDS(pred_ss7, file = paste("/projectnb/cp2018ssml/workdir/jimmy/models/",args[1],"_ss7_prediction.rds", sep=""))

pdf(paste("/projectnb/cp2018ssml/workdir/jimmy/roc/",args[1],"_ss7_ROC.pdf", sep=""))
plot(roc_7, main=paste("ROC Plot for LumA vs LumB in semi-supervised ss7", sep=""))
text(x=0.2, y=0.2, labels=paste("AUC:", auc(roc_7)))
dev.off()
auc_vec <- c(auc_vec, auc(roc_7))

saveRDS(auc_vec, file = paste("/projectnb/cp2018ssml/workdir/jimmy/","auc_vector_", args[1], ".rds", sep=""))
