# Print filename
print(c("File used: ", commandArgs(trailingOnly=FALSE)[4]));
# Pull the species_name from the command line argument

cmd_args = commandArgs(trailingOnly=TRUE);
species_name = paste(cmd_args[1], sep = "");

print(c("SPECIES NAME: ", species_name)); 


# Random forest package.
library(randomForest)
# verification package is needed for the roc.area function used in class.sum.
library(verification)


#################
# This is the class.sum code that R. Cutler wrote to calculate some useful performance metrics.
#################
class.sum=function(x){
  pcc=100*sum(diag(x))/sum(x)
  spec=100*x[1,1]/sum(x[1,])
  sens=100*x[2,2]/sum(x[2,])
  youden=(spec + sens - 100)
  return(cbind(c("Percent Correctly Classified = ","Specificity = ","Sensitivity = ","Youden Index ="),c(pcc,spec,sens,youden)))
}

kappa=function(x){
  n=sum(x)
  pobs=(x[1,1]+x[2,2])/n
  pexp=(sum(x[1,])*sum(x[,1])+sum(x[2,])*sum(x[,2]))/n^2
  kappa=(pobs-pexp)/(1-pexp)
  return(kappa)
}

#################
# Hyperparameters to grid search through during model optimization.
#################
mtry = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
ntree = c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800)

# initialize mean_Youden_list.
k = 1
# mean_Youden_list returns the mean (of 5 crossvalidation folds) for each set of hyperparameters.
# Each model is optimized based on max Youden Index = Sensitivity + Specificity - 1.
mean_Youden_list=list()

#################
# The 5-fold cross validation.
# Runs the model on the 5 separate train/test datasets and output max mean Youden Index.
# Each dataset is on a fixed 70/30 train/test split with prevalence preserved between train and test datasets.
#################

for (i in mtry){
  for (j in ntree){
  try({
    
    #################
    # train 1
    #################
    
    train1 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "1_TRAIN.csv", sep=""));
    val1 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "1_TEST.csv", sep=""));
    
    predictors1 = train1[1:11]
    response1 = as.factor(train1[[species_name]])
    upsampled_data1 = upSample(predictors1, response1, list=FALSE, yname="response1")
    table(upsampled_data1$response1) 
    response1<-as.factor(upsampled_data1$response1) 
    
    taxa_val1 <- as.factor(val1[[species_name]])
    
    set.seed(422)
    RF_model1 = randomForest(response1 ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM,
                            data=upsampled_data1,
                            mtry = i,
                            ntree = j)
    
    # This code is used to predict on the validation dataset and determine final model performance.
    RF_model.prediction1 = predict(RF_model1, val1, type = "response")
    RF_model.prob1 = predict(RF_model1, val1, type = "prob")
    # Extract the probabilities of presence (i.e. column of probabilities with header 1)
    RF_model.prob1_presence <- RF_model.prob1[,'1']
    
    # Validation confusion matrix.
    confusion_matrix1=table(taxa_val1, RF_model.prediction1)
    
    numeric_taxa_val1 <- as.numeric(as.character(taxa_val1))
    
    # Performance metrics.
    metric1 <- class.sum(confusion_matrix1)
    Youden1 <- metric1[4,2]
    Youden1 <- as.numeric(Youden1)
    AUC1 <- roc.area(numeric_taxa_val1, RF_model.prob1_presence)$A
    AUC1 <-as.numeric(AUC1)
    Sens1 <-metric1[3,2]
    Sens1 <- as.numeric(Sens1)
    PCC1 <- metric1[1,2]
    PCC1 <- as.numeric(PCC1)
    Kappa1<-kappa(confusion_matrix1)
    
    #################
    # train 2
    #################
  
    train2 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "2_TRAIN.csv", sep=""));
    val2 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "2_TEST.csv", sep=""));
    
    predictors2 = train2[1:11]
    response2 = as.factor(train2[[species_name]])
    upsampled_data2 = upSample(predictors2, response2, list=FALSE, yname="response2")
    table(upsampled_data2$response2) 
    response2<-as.factor(upsampled_data2$response2) 
    
    taxa_val2 <- as.factor(val2[[species_name]])

    set.seed(422)
    RF_model2 = randomForest(response2 ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM,
                            data=upsampled_data2, 
                            mtry = i,
                            ntree = j)
    
    # This code is used to predict on the validation dataset and determine final model performance.
    RF_model.prediction2 = predict(RF_model2, val2, type = "response")
    RF_model.prob2 = predict(RF_model2, val2, type = "prob")
    # Extract the probabilities of presence (i.e. column of probabilities with header 1)
    RF_model.prob2_presence <- RF_model.prob2[,'1']
    
    # Validation confusion matrix.
    confusion_matrix2=table(taxa_val2, RF_model.prediction2)
    
    numeric_taxa_val2 <- as.numeric(as.character(taxa_val2))
    
    # Performance metrics.
    metric2 <- class.sum(confusion_matrix2)
    Youden2 <- metric2[4,2]
    Youden2 <- as.numeric(Youden2)
    AUC2 <- roc.area(numeric_taxa_val2, RF_model.prob2_presence)$A
    AUC2 <-as.numeric(AUC2)
    Sens2 <-metric2[3,2]
    Sens2 <- as.numeric(Sens2)
    PCC2 <- metric2[1,2]
    PCC2 <- as.numeric(PCC2)
    Kappa2<-kappa(confusion_matrix2)
    
    #################
    # train 3
    #################
  
    train3 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "3_TRAIN.csv", sep=""));
    val3 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "3_TEST.csv", sep=""));
    
    predictors3 = train3[1:11]
    response3 = as.factor(train3[[species_name]])
    upsampled_data3 = upSample(predictors3, response3, list=FALSE, yname="response3")
    table(upsampled_data3$response3) 
    response3<-as.factor(upsampled_data3$response3) 

    taxa_val3 <- as.factor(val3[[species_name]])
    
    set.seed(422)
    RF_model3 = randomForest(response3 ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM,
                            data=upsampled_data3, 
                            mtry = i,
                            ntree = j)
    
    # This code is used to predict on the validation dataset and determine final model performance.
    RF_model.prediction3 = predict(RF_model3, val3, type = "response")
    RF_model.prob3 = predict(RF_model3, val3, type = "prob")
    # Extract the probabilities of presence (i.e. column of probabilities with header 1)
    RF_model.prob3_presence <- RF_model.prob3[,'1']
    
    # Validation confusion matrix.
    confusion_matrix3=table(taxa_val3, RF_model.prediction3)
    
    numeric_taxa_val3 <- as.numeric(as.character(taxa_val3))
    
    # Performance metrics.
    metric3 <- class.sum(confusion_matrix3)
    Youden3 <- metric3[4,2]
    Youden3 <- as.numeric(Youden3)
    AUC3 <- roc.area(numeric_taxa_val3, RF_model.prob3_presence)$A
    AUC3 <-as.numeric(AUC3)
    Sens3 <-metric3[3,2]
    Sens3 <- as.numeric(Sens3)
    PCC3 <- metric3[1,2]
    PCC3 <- as.numeric(PCC3)
    Kappa3<-kappa(confusion_matrix3)

    #################
    # train 4
    #################
  
    train4 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "4_TRAIN.csv", sep=""));
    val4 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "4_TEST.csv", sep=""));
    
    predictors4 = train4[1:11]
    response4 = as.factor(train4[[species_name]])
    upsampled_data4 = upSample(predictors4, response4, list=FALSE, yname="response4")
    table(upsampled_data4$response4) 
    response4<-as.factor(upsampled_data4$response4) 

    taxa_val4 <- as.factor(val4[[species_name]])
    
    set.seed(422)
    RF_model4 = randomForest(response4 ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM,
                            data=upsampled_data4, 
                            mtry = i,
                            ntree = j)
    
    # This code is used to predict on the validation dataset and determine final model performance.
    RF_model.prediction4 = predict(RF_model4, val4, type = "response")
    RF_model.prob4 = predict(RF_model4, val4, type = "prob")
    # Extract the probabilities of presence (i.e. column of probabilities with header 1)
    RF_model.prob4_presence <- RF_model.prob4[,'1']
    
    # Validation confusion matrix.
    confusion_matrix4=table(taxa_val4, RF_model.prediction4)
    
    numeric_taxa_val4 <- as.numeric(as.character(taxa_val4))
    
    # Performance metrics.
    metric4 <- class.sum(confusion_matrix4)
    Youden4 <- metric4[4,2]
    Youden4 <- as.numeric(Youden4)
    AUC4 <- roc.area(numeric_taxa_val4, RF_model.prob4_presence)$A
    AUC4 <-as.numeric(AUC4)
    Sens4 <-metric4[3,2]
    Sens4 <- as.numeric(Sens4)
    PCC4 <- metric4[1,2]
    PCC4 <- as.numeric(PCC4)
    Kappa4<-kappa(confusion_matrix4)
    
    #################
    # train 5
    #################
  
    train5 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "5_TRAIN.csv", sep=""));
    val5 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "5_TEST.csv", sep=""));
    
    predictors5 = train5[1:11]
    response5 = as.factor(train5[[species_name]])
    upsampled_data5 = upSample(predictors5, response5, list=FALSE, yname="response5")
    table(upsampled_data5$response5) 
    response5<-as.factor(upsampled_data5$response5) 

    taxa_val5 <- as.factor(val5[[species_name]])
    
    set.seed(422)
    RF_model5 = randomForest(response5 ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM,
                            data=upsampled_data5, 
                            mtry = i,
                            ntree = j)
    
    # This code is used to predict on the validation dataset and determine final model performance.
    RF_model.prediction5 = predict(RF_model5, val5, type = "response")
    RF_model.prob5 = predict(RF_model5, val5, type = "prob")
    # Extract the probabilities of presence (i.e. column of probabilities with header 1)
    RF_model.prob5_presence <- RF_model.prob5[,'1']
    
    # Validation confusion matrix.
    confusion_matrix5=table(taxa_val5,RF_model.prediction5)

    numeric_taxa_val5 <- as.numeric(as.character(taxa_val5))
    
    # Performance metrics.
    metric5 <- class.sum(confusion_matrix5)
    Youden5 <- metric5[4,2]
    Youden5 <- as.numeric(Youden5)
    AUC5 <- roc.area(numeric_taxa_val5, RF_model.prob5_presence)$A
    AUC5 <-as.numeric(AUC5)
    Sens5 <-metric5[3,2]
    Sens5 <- as.numeric(Sens5)
    PCC5 <- metric5[1,2]
    PCC5 <- as.numeric(PCC5)
    Kappa5<-kappa(confusion_matrix5)
    
    #################
    # Mean Youden Index per loop.
    #################
    
    Mean_Youden = (Youden1+Youden2+Youden3+Youden4+Youden5)/5
    mean_Youden_list[k]=Mean_Youden
    k = k + 1
    
    Mean_Kappa = (Kappa1+Kappa2+Kappa3+Kappa4+Kappa5)/5
    Mean_AUC = (AUC1+AUC2+AUC3+AUC4+AUC5)/5
    Mean_Sens = (Sens1+Sens2+Sens3+Sens4+Sens5)/5
    Mean_PCC = (PCC1+PCC2+PCC3+PCC4+PCC5)/5
    
    print(cbind(c("mtry","ntrees"),c(i,j)))
    print("Mean Youden")
    print(Mean_Youden)
    print("Mean Kappa")
    print(Mean_Kappa)
    print("Mean AUC")
    print(Mean_AUC)
    print("Mean Sens")
    print(Mean_Sens)
    print("Mean PCC")
    print(Mean_PCC)
    print(".")
    print(".")
    }, silent = FALSE)
  }  
}  

# Prints the mean Youden Index value for each set of hyperparameters (i.e. models).
# Max Youden Index value in this list corresponds to optimal model and set of hyperparameters. 
print("mean_Youden_list");
print(mean_Youden_list)

# Returns the highest Youden Index value from the 5-fold cross validation.
# Then query the raw output text file to find corresponding hyperparameters if desired (or scroll through R-console output).
print("max(unlist(mean_Youden_list))");
max(unlist(mean_Youden_list))
