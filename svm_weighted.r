# Print filename
print(c("File used: ", commandArgs(trailingOnly=FALSE)[4]));
# Pull the species_name from the command line argument

cmd_args = commandArgs(trailingOnly=TRUE);
species_name = paste(cmd_args[1], sep = "");

print(c("SPECIES NAME: ", species_name)); 


# The core SVM package.
library(e1071)
# verification package is needed for the roc.area function used in class.sum.
library(verification)

#################
# This is the class.sum and kappa code that R. Cutler wrote to calculate some useful performance metrics.
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
gamma = c(0.001, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3 , 0.33, 0.36, 0.39, 0.43, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67, 0.7 , 0.73, 0.76, 0.79, 0.82, 0.85, 0.88, 0.91, 0.94, 0.97, 1.0 , 1.03, 1.06, 1.09, 1.12, 1.15, 1.18, 1.21, 1.24, 1.27, 1.3 , 1.33, 1.36, 1.39, 1.42, 1.46, 1.49, 1.52, 1.55, 1.58, 1.61, 1.64, 1.67, 1.7 , 1.73, 1.76, 1.79, 1.82, 1.85, 1.88, 1.91, 1.94, 1.97, 2.0 , 2.03, 2.06, 2.09, 2.12, 2.15, 2.18, 2.21, 2.24, 2.27, 2.3 , 2.33, 2.36, 2.39, 2.42, 2.45, 2.49, 2.52, 2.55, 2.58, 2.61, 2.64, 2.67, 2.7 , 2.73, 2.76, 2.79, 2.82, 2.85, 2.88, 2.91, 2.94, 2.97, 3.0)
cost = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0, 15.1, 15.2, 15.3, 15.4, 15.5, 15.6, 15.7, 15.8, 15.9, 16.0, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 17.0, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.8, 17.9, 18.0, 18.1, 18.2, 18.3, 18.4, 18.5, 18.6, 18.7, 18.8, 18.9, 19.0, 19.1, 19.2, 19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9, 20.0)

# Initialize mean_Youden_list.
k = 1
# mean_Youden_list returns the mean (of 5 crossvalidation folds) for each set of hyperparameters.
# Each model is optimized based on max Youden Index = Sensitivity + Specificity - 1.
mean_Youden_list=list()

#################
# The 5-fold cross validation.
# Runs the model on the 5 separate train/test datasets and output max mean Youden Index.
# Each dataset is on a fixed 70/30 train/test split with prevalence preserved between train and test datasets.
#################

for (i in gamma){
  for (j in cost){
  try({
    
    #################
    # train 1
    #################
    
    train1 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "1_TRAIN.csv", sep=""));
    val1 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "1_TEST.csv", sep=""));
    
    taxa_train1 <- as.factor(train1[[species_name]])
    taxa_val1 <- as.factor(val1[[species_name]])
    
    class_instances1<-table(taxa_train1)
    absences1 = class_instances1[1]
    presences1 = class_instances1[2]
    number_presences1 = as.numeric(class_instances1[2])
    total1 = absences1 + presences1
    
    presences_fraction1 = presences1/total1
    presences_weight1 = 1 - presences_fraction1
    absences_fraction1 = absences1/total1
    absences_weight1 = 1 - absences_fraction1
    
    # Train the SVM.
    set.seed(422)
    SVM_model1 = svm(as.factor(taxa_train1) ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM, 
                    gamma=i,
                    cost=j,
                    class.weights=c(absences_weight1, presences_weight1),
                    data=train1)
    
    # This code is used to predict on the validation dataset and determine final model performance.
    SVM_model.prediction1 = predict(SVM_model1, val1, type = "response")

    # Validation confusion matrix.
    confusion_matrix1=table(taxa_val1, SVM_model.prediction1)
    
    numeric_SVM_model.prediction1 <- as.numeric(as.character(SVM_model.prediction1))
    numeric_taxa_val1 <- as.numeric(as.character(taxa_val1))
    
    # Performance metrics.
    metric1 <- class.sum(confusion_matrix1)
    Youden1 <- metric1[4,2]
    Youden1 <- as.numeric(Youden1)
    AUC1 <- roc.area(numeric_taxa_val1, numeric_SVM_model.prediction1)$A
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
    
    taxa_train2 <- as.factor(train2[[species_name]])
    taxa_val2 <- as.factor(val2[[species_name]])
    
    class_instances2<-table(taxa_train2)
    absences2 = class_instances2[1]
    presences2 = class_instances2[2]
    number_presences2 = as.numeric(class_instances2[2])
    total2 = absences2 + presences2
    
    presences_fraction2 = presences2/total2
    presences_weight2 = 1 - presences_fraction2
    absences_fraction2 = absences2/total2
    absences_weight2 = 1 - absences_fraction2
    
    set.seed(422)
    SVM_model2 = svm(as.factor(taxa_train2) ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM, 
                    gamma=i,
                    cost=j,
                    class.weights=c(absences_weight2, presences_weight2),
                    data=train2)
    
    # This code is used to predict on the validation dataset and determine final model performance.
    SVM_model.prediction2 = predict(SVM_model2, val2, type = "response")
    
    # Validation confusion matrix.
    confusion_matrix2=table(taxa_val2, SVM_model.prediction2)
    
    numeric_SVM_model.prediction2 <- as.numeric(as.character(SVM_model.prediction2))
    numeric_taxa_val2 <- as.numeric(as.character(taxa_val2))
    
    # Performance metrics.
    metric2 <- class.sum(confusion_matrix2)
    Youden2 <- metric2[4,2]
    Youden2 <- as.numeric(Youden2)
    AUC2 <- roc.area(numeric_taxa_val2, numeric_SVM_model.prediction2)$A
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
    
    taxa_train3 <- as.factor(train3[[species_name]])
    taxa_val3 <- as.factor(val3[[species_name]])
    
    class_instances3<-table(taxa_train3)
    absences3 = class_instances3[1]
    presences3 = class_instances3[2]
    number_presences3 = as.numeric(class_instances3[2])
    total3 = absences3 + presences3
    
    presences_fraction3 = presences3/total3
    presences_weight3 = 1 - presences_fraction3
    absences_fraction3 = absences3/total3
    absences_weight3 = 1 - absences_fraction3

    set.seed(422)
    SVM_model3 = svm(as.factor(taxa_train3) ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM, 
                    gamma=i,
                    cost=j,
                    class.weights=c(absences_weight3, presences_weight3),
                    data=train3)
    
    # This code is used to predict on the validation dataset and determine final model performance.
    SVM_model.prediction3 = predict(SVM_model3, val3, type = "response")
    
    # Validation confusion matrix.
    confusion_matrix3=table(taxa_val3, SVM_model.prediction3)
    
    numeric_SVM_model.prediction3 <- as.numeric(as.character(SVM_model.prediction3))
    numeric_taxa_val3 <- as.numeric(as.character(taxa_val3))
    
    # Performance metrics.
    metric3 <- class.sum(confusion_matrix3)
    Youden3 <- metric3[4,2]
    Youden3 <- as.numeric(Youden3)
    AUC3 <- roc.area(numeric_taxa_val3, numeric_SVM_model.prediction3)$A
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
    
    taxa_train4 <- as.factor(train4[[species_name]])
    taxa_val4 <- as.factor(val4[[species_name]])
    
    class_instances4<-table(taxa_train4)
    absences4 = class_instances4[1]
    presences4 = class_instances4[2]
    number_presences4 = as.numeric(class_instances4[2])
    total4 = absences4 + presences4
    
    presences_fraction4 = presences4/total4
    presences_weight4 = 1 - presences_fraction4
    absences_fraction4 = absences4/total4
    absences_weight4 = 1 - absences_fraction4
    
    set.seed(422)
    SVM_model4 = svm(as.factor(taxa_train4) ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM, 
                    gamma=i,
                    cost=j,
                    class.weights=c(absences_weight4, presences_weight4),
                    data=train4)
    
    # This code is used to predict on the validation dataset and determine final model performance.
    SVM_model.prediction4 = predict(SVM_model4, val4, type = "response")
    
    # Validation confusion matrix.
    confusion_matrix4=table(taxa_val4, SVM_model.prediction4)
    
    numeric_SVM_model.prediction4 <- as.numeric(as.character(SVM_model.prediction4))
    numeric_taxa_val4 <- as.numeric(as.character(taxa_val4))
    
    # Performance metrics.
    metric4 <- class.sum(confusion_matrix4)
    Youden4 <- metric4[4,2]
    Youden4 <- as.numeric(Youden4)
    AUC4 <- roc.area(numeric_taxa_val4, numeric_SVM_model.prediction4)$A
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
    
    taxa_train5 <- as.factor(train5[[species_name]])
    taxa_val5 <- as.factor(val5[[species_name]])
    
    class_instances5<-table(taxa_train5)
    absences5 = class_instances5[1]
    presences5 = class_instances5[2]
    number_presences5 = as.numeric(class_instances5[2])
    total5 = absences5 + presences5
    
    presences_fraction5 = presences5/total5
    presences_weight5 = 1 - presences_fraction5
    absences_fraction5 = absences5/total5
    absences_weight5 = 1 - absences_fraction5
    
    set.seed(422)
    SVM_model5 = svm(as.factor(taxa_train5) ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM, 
                    gamma=i,
                    cost=j,
                    class.weights=c(absences_weight5, presences_weight5),
                    data=train5)
    
    # This code is used to predict on the validation dataset and determine final model performance.
    SVM_model.prediction5 = predict(SVM_model5, val5, type = "response")
    
    # Validation confusion matrix.
    confusion_matrix5=table(taxa_val5, SVM_model.prediction5)
    
    numeric_SVM_model.prediction5 <- as.numeric(as.character(SVM_model.prediction5))
    numeric_taxa_val5 <- as.numeric(as.character(taxa_val5))
    
    # Performance metrics.
    metric5 <- class.sum(confusion_matrix5)
    Youden5 <- metric5[4,2]
    Youden5 <- as.numeric(Youden5)
    AUC5 <- roc.area(numeric_taxa_val5, numeric_SVM_model.prediction5)$A
    AUC5 <-as.numeric(AUC5)
    Sens5 <-metric5[3,2]
    Sens5 <- as.numeric(Sens5)
    PCC5 <- metric5[1,2]
    PCC5 <- as.numeric(PCC5)
    Kappa5<-kappa(confusion_matrix5)
    
    #################
    # Mean Youden Index and other performance metrics per loop.
    #################
    
    Mean_Youden = (Youden1+Youden2+Youden3+Youden4+Youden5)/5
    mean_Youden_list[k]=Mean_Youden
    k = k + 1
    
    Mean_Kappa = (Kappa1+Kappa2+Kappa3+Kappa4+Kappa5)/5
    Mean_AUC = (AUC1+AUC2+AUC3+AUC4+AUC5)/5
    Mean_Sens = (Sens1+Sens2+Sens3+Sens4+Sens5)/5
    Mean_PCC = (PCC1+PCC2+PCC3+PCC4+PCC5)/5
    
    print(cbind(c("gamma","cost"),c(i,j)))
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
    }, silent = TRUE)
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
