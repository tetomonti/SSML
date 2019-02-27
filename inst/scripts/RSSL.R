library(tidyverse)
library(SummarizedExperiment)
library(ssc)
library(RSSL)

load("/Users/JimmyKiely/Documents/Monti Lab/BRCA.rda")

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


BRCA <- prepare_data(BRCA)

LUM_A <-
    dplyr::filter(BRCA, PAM50.mRNA == "Luminal A") %>%
    provision_data("PAM50.mRNA", 0.70, 0.30, 1001)

LUM_B <-
    dplyr::filter(BRCA, PAM50.mRNA == "Luminal B") %>%
    provision_data("PAM50.mRNA", 0.70, 0.30, 1001)

x <- join_data(LUM_A, LUM_B)

training <- magrittr::extract2(x, 1)

testing <- magrittr::extract2(x, 2)

#Get top 5000 genes by MAD in training
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
#______________________________________________________________________________________________
# matrix for storing final results
final_auc_results <- matrix(NA, nrow=8, ncol=10)
rownames(final_auc_results) <- c("FL", "1/7 UL", "2/7 UL", "3/7 UL", "4/7 UL", "5/7 UL",
                                 "6/7 UL", "7/7 UL")

# svm for fl + no unlabeled
model_fl <- SVM(X=training_final, y=training_labels, scale=T)
pred_fl <- predict(model_fl, newdata = test_exp_5000)
d_values <- decisionvalues(model_fl, test_exp_5000)
roc <- roc(test_labels, d_values)
plot(smooth(roc), main='ROC SVM Fully Labeled')
text(x=0.2, y=0.2, labels=paste("AUC:",auc(roc)))

# ss svm for the fl + 1/7 of the data 
unlabeled_ss1 <- ul_training_splits[[1]][,2686:ncol(ul_training_splits[[1]])]
unlabeled_ss1 <- unlabeled_ss1[,colnames(unlabeled_ss1) %in% names(top5000_training)]
model_ss1 <- WellSVM(X=training_final, y=training_labels, X_u=unlabeled_ss1, scale=T, gamma=0)
pred_ss1 <- predict(model_ss1, newdata = test_exp_5000, probs=TRUE)
d_values_1 <- decisionvalues(model_ss1, test_exp_5000)
roc_1 <- roc(test_labels, d_values_1)
plot(smooth(roc_1), main='ROC SVM Labeled + 1/7 Unlabeled')
text(x=0.2, y=0.2, labels=paste("AUC:",auc(roc_1)))

# ss svm for the fl + 2/7 of the data 
unlabeled_ss2 <- ul_training_splits[[2]][,2686:ncol(ul_training_splits[[2]])]
unlabeled_ss2 <- unlabeled_ss2[,colnames(unlabeled_ss2) %in% names(top5000_training)]
model_ss2 <- WellSVM(X=training_final, y=training_labels, X_u=unlabeled_ss2, scale=T, gamma=0)
pred_ss2 <- predict(model_ss2, newdata = test_exp_5000)
d_values_2 <- decisionvalues(model_ss2, test_exp_5000)
roc_2 <- roc(test_labels, d_values_2)
plot(smooth(roc_2), main='ROC SVM Labeled + 2/7 Unlabeled')
text(x=0.2, y=0.2, labels=paste("AUC:",auc(roc_2)))

# ss svm for the fl + 3/7 of the data 
unlabeled_ss3 <- ul_training_splits[[3]][,2686:ncol(ul_training_splits[[3]])]
unlabeled_ss3 <- unlabeled_ss3[,colnames(unlabeled_ss3) %in% names(top5000_training)]
model_ss3 <- WellSVM(X=training_final, y=training_labels, X_u=unlabeled_ss3, scale=T, gamma=0)
pred_ss3 <- predict(model_ss3, newdata = test_exp_5000)
d_values_3 <- decisionvalues(model_ss3, test_exp_5000)
roc_3 <- roc(test_labels, d_values_3)
plot(smooth(roc_3), main='ROC SVM Labeled + 3/7 Unlabeled')
text(x=0.2, y=0.2, labels=paste("AUC:",auc(roc_3)))

# ss svm for the fl + 4/7 of the data 
unlabeled_ss4 <- ul_training_splits[[4]][,2686:ncol(ul_training_splits[[4]])]
unlabeled_ss4 <- unlabeled_ss4[,colnames(unlabeled_ss4) %in% names(top5000_training)]
model_ss4 <- WellSVM(X=training_final, y=training_labels, X_u=unlabeled_ss4, scale=T, gamma=0)
pred_ss4 <- predict(model_ss4, newdata = test_exp_5000)
d_values_4<- decisionvalues(model_ss4, test_exp_5000)
roc_4 <- roc(test_labels, d_values_4)
plot(smooth(roc_4), main='ROC SVM Labeled + 4/7 Unlabeled')
text(x=0.2, y=0.2, labels=paste("AUC:",auc(roc_4)))

# ss svm for the fl + 5/7 of the data 
unlabeled_ss5 <- ul_training_splits[[5]][,2686:ncol(ul_training_splits[[5]])]
unlabeled_ss5 <- unlabeled_ss5[,colnames(unlabeled_ss5) %in% names(top5000_training)]
model_ss5 <- WellSVM(X=training_final, y=training_labels, X_u=unlabeled_ss5, scale=T, gamma=0)
pred_ss5 <- predict(model_ss5, newdata = test_exp_5000)
d_values_5 <- decisionvalues(model_ss5, test_exp_5000)
roc_5 <- roc(test_labels, d_values_5)
plot(smooth(roc_5), main='ROC SVM Labeled + 5/7 Unlabeled')
text(x=0.2, y=0.2, labels=paste("AUC:",auc(roc_5)))

# ss svm for the fl + 6/7 of the data 
unlabeled_ss6 <- ul_training_splits[[6]][,2686:ncol(ul_training_splits[[6]])]
unlabeled_ss6 <- unlabeled_ss6[,colnames(unlabeled_ss6) %in% names(top5000_training)]
model_ss6 <- WellSVM(X=training_final, y=training_labels, X_u=unlabeled_ss6, scale=T, gamma=0)
pred_ss6 <- predict(model_ss6, newdata = test_exp_5000)
d_values_6 <- decisionvalues(model_ss6, test_exp_5000)
roc_6 <- roc(test_labels, d_values_6)
plot(smooth(roc_6), main='ROC SVM Labeled + 6/7 Unlabeled')
text(x=0.2, y=0.2, labels=paste("AUC:",auc(roc_6)))

# ss svm for the fl + 7/7 of the data 
unlabeled_ss7 <- ul_training_splits[[7]][,2686:ncol(ul_training_splits[[7]])]
unlabeled_ss7 <- unlabeled_ss7[,colnames(unlabeled_ss7) %in% names(top5000_training)]
model_ss7 <- WellSVM(X=training_final, y=training_labels, X_u=unlabeled_ss7, scale=T, gamma=0)
pred_ss7 <- predict(model_ss7, newdata = test_exp_5000)
d_values_7 <- decisionvalues(model_ss7, test_exp_5000)
roc_7 <- roc(test_labels, d_values_7)
plot(smooth(roc_7), main='ROC SVM Labeled + All Unlabeled')
text(x=0.2, y=0.2, labels=paste("AUC:",auc(roc_7)))

