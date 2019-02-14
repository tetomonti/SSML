##load data
library(Biobase)
library(ssc)
library(randomForest)

##note arg1 = input file name for training, arg2 = input file name for testing, arg[3] = name of split (i.e. 1D) for output files
args = commandArgs(trailingOnly=TRUE)
training <- load(paste("/projectnb/cp2018ssml/datasets/tcga/", args[1],".rda", sep=""))
testing <- load(paste("/projectnb/cp2018ssml/datasets/tcga/", args[2],".rda", sep=""))

#get only expression data
train_ss_exp <- training[,2686:ncol(training)]
## get labels and convert to binary -- Luminal A= 0 and Luminal B=1 
train_ss_labels <- as.factor(as.numeric(training$PAM50.mRNA=="Luminal B"))

#get top 5000 genes by median absolute deviation
colmad <- apply(train_ss_exp, 2, mad)
top5000 <- tail(sort(colmad),5000)

#subset top 5000 genes in training
train_ss_exp_5000 <- train_ss_exp[,colnames(train_ss_exp) %in% names(top5000) ]

##get test data
#only expression
test_exp <- testing[,2686:ncol(testing)]
#subset top 5000 genes in testing
test_exp_5000 <- test_exp[,colnames(test_exp) %in% names(top5000) ]
## get labels and convert to binary -- Luminal A= 0 and Luminal B=1 
test_labels <- as.factor(as.numeric(testing$PAM50.mRNA=="Luminal B"))

#remove all NA samples (fully labeled dataset)
training_fl <- training[!is.na(training$PAM50.mRNA),]
#get only expression
train_fl_exp <-  training_fl[,2686:ncol(training_fl)]
#get top 5000 genes
train_fl_exp_5000 <- train_fl_exp[,colnames(train_fl_exp) %in% names(top5000) ]
#get lables in binary
train_fl_labels <- as.factor(as.numeric(training_fl$PAM50.mRNA=="Luminal B"))

#set up training for ssc self training (with unlabeled data)
model1 <- selfTraining(x = train_ss_exp_5000, y = train_ss_labels,
                    learner = randomForest,
                    pred = "predict", pred.pars = list(type="prob"))
#get predictions
pred <- predict(model1, as.matrix(test_exp_5000))

#save semi-supervised results
saveRDS(model1, file = paste("/projectnb/cp2018ssml/workdir/rebecca/",arg[3],"_ss_training_model.rds", sep=""))
saveRDS(pred, file = paste("/projectnb/cp2018ssml/workdir/rebecca/",arg[3],"_ss_prediction_model.rds", sep=""))

#run fully supervised rf
model_fl <- randomForest(x=train_fl_exp, y=train_fl_labels)
pred_fl <- predict(model, newdata = test_exp)

#save predictions and model 
saveRDS(model_fl, file = paste("/projectnb/cp2018ssml/workdir/rebecca/",arg[3],"fl_training_model.rds", sep=""))
saveRDS(pred_fl, file = paste("/projectnb/cp2018ssml/workdir/rebecca/1D_fl_prediction_model.rds", sep=""))
