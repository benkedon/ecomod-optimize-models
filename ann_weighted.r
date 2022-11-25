# Print filename
print(c("File used: ", commandArgs(trailingOnly=FALSE)[4]));
# Pull the species_name from the command line argument

cmd_args = commandArgs(trailingOnly=TRUE);
species_name = paste(cmd_args[1], sep = "");

print(c("SPECIES NAME: ", species_name)); 


# The core neural network package.
library(nnet)
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

# Size is the number of neurons in hidden layer.
size = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100)

# maxit is the number of iterations during training.
maxit = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780, 790, 800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900, 910, 920, 930, 940, 950, 960, 970, 980, 990, 1000)

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

for (i in size){
  for (j in maxit){
    try({
    #################
    # train 1
    #################
    
    train1 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "1_TRAIN.csv", sep=""));
    val1 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "1_TEST.csv", sep=""));
    
    taxa_train1 <- as.factor(train1[[species_name]])
    taxa_val1 <- as.factor(val1[[species_name]])
    
    # determining class weights.
    class_instances1<-table(taxa_train1)
    absences1 = class_instances1[1]
    presences1 = class_instances1[2]
    total1 = absences1 + presences1
    
    presences_fraction1 = presences1/total1
    presences_weight1 = 1 - presences_fraction1
    absences_fraction1 = absences1/total1
    absences_weight1 = 1 - absences_fraction1
    
    myWeight1 <- rep(NA,length(taxa_train1))
    myWeight1[taxa_train1==1] <- presences_weight1
    myWeight1[taxa_train1==0] <- absences_weight1
    
    set.seed(422)
    ANN_model1 = nnet(taxa_train1 ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM,
                      trace=FALSE,
                      data=train1,
                      MaxNWts = 2000,
                      weights = myWeight1,
                      size = i,
                      maxit = j)
    
    # This code is used to predict on the validation dataset and determine final model performance.
    ANN_model.prediction1 = predict(ANN_model1, val1, type = "class")
    ANN_model.probability1 = predict(ANN_model1, val1, type = "raw")
    # we must declare the levels 1 and 0. Otherwise models without any predicted absences will not be included as a level and the confusion matrix will be incomplete.
    final_ANN_model.prediction1<-factor(ANN_model.prediction1, levels = c(0, 1))
    
    # Validation confusion matrix.
    confusion_matrix1=table(taxa_val1, final_ANN_model.prediction1)
    
    # roc.area function requires validation taxa classes (0 or 1) to be numeric.
    numeric_taxa_val1 <- as.numeric(as.character(taxa_val1))
    
    # Performance metrics.
    metric1 <- class.sum(confusion_matrix1)
    Youden1 <- metric1[4,2]
    Youden1 <- as.numeric(Youden1)
    AUC1 <- roc.area(numeric_taxa_val1, ANN_model.probability1)$A
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
    total2 = absences2 + presences2
    
    presences_fraction2 = presences2/total2
    presences_weight2 = 1 - presences_fraction2
    absences_fraction2 = absences2/total2
    absences_weight2 = 1 - absences_fraction2
    
    myWeight2 <- rep(NA,length(taxa_train2))
    myWeight2[taxa_train2==1] <- presences_weight2
    myWeight2[taxa_train2==0] <- absences_weight2
    
    set.seed(422)
    ANN_model2 = nnet(taxa_train2 ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM,
                      trace=FALSE,
                      data=train2,
                      MaxNWts = 2000,
                      weights = myWeight2,
                      size = i,
                      maxit = j)
    
    # This code is used to predict on the validation dataset and determine final model performance.
    ANN_model.prediction2 = predict(ANN_model2, val2, type = "class")
    ANN_model.probability2 = predict(ANN_model2, val2, type = "raw")
    # we must declare the levels 1 and 0. Otherwise models without any predicted absences will not be included as a level and the confusion matrix will be incomplete.
    final_ANN_model.prediction2<-factor(ANN_model.prediction2, levels = c(0, 1))
    
    # Validation confusion matrix.
    confusion_matrix2=table(taxa_val2, final_ANN_model.prediction2)
    
    # roc.area function requires validation taxa classes (0 or 1) to be numeric.
    numeric_taxa_val2 <- as.numeric(as.character(taxa_val2))
    
    # Performance metrics.
    metric2 <- class.sum(confusion_matrix2)
    Youden2 <- metric2[4,2]
    Youden2 <- as.numeric(Youden2)
    AUC2 <- roc.area(numeric_taxa_val2, ANN_model.probability2)$A
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
    total3 = absences3 + presences3
    
    presences_fraction3 = presences3/total3
    presences_weight3 = 1 - presences_fraction3
    absences_fraction3 = absences3/total3
    absences_weight3 = 1 - absences_fraction3
    
    myWeight3 <- rep(NA,length(taxa_train3))
    myWeight3[taxa_train3==1] <- presences_weight3
    myWeight3[taxa_train3==0] <- absences_weight3

    set.seed(422)
    ANN_model3 = nnet(taxa_train3 ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM,
                      trace=FALSE,
                      data=train3,
                      MaxNWts = 2000,
                      weights = myWeight3,
                      size = i,
                      maxit = j)
    
    # This code is used to predict on the validation dataset and determine final model performance.
    ANN_model.prediction3 = predict(ANN_model3, val3, type = "class")
    ANN_model.probability3 = predict(ANN_model3, val3, type = "raw")
    # we must declare the levels 1 and 0. Otherwise models without any predicted absences will not be included as a level and the confusion matrix will be incomplete.
    final_ANN_model.prediction3<-factor(ANN_model.prediction3, levels = c(0, 1))
    
    # Validation confusion matrix.
    confusion_matrix3=table(taxa_val3, final_ANN_model.prediction3)
    
    # roc.area function requires validation taxa classes (0 or 1) to be numeric.
    numeric_taxa_val3 <- as.numeric(as.character(taxa_val3))
    
    # Performance metrics.
    metric3 <- class.sum(confusion_matrix3)
    Youden3 <- metric3[4,2]
    Youden3 <- as.numeric(Youden3)
    AUC3 <- roc.area(numeric_taxa_val3, ANN_model.probability3)$A
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
    total4 = absences4 + presences4
    
    presences_fraction4 = presences4/total4
    presences_weight4 = 1 - presences_fraction4
    absences_fraction4 = absences4/total4
    absences_weight4 = 1 - absences_fraction4
    
    myWeight4 <- rep(NA,length(taxa_train4))
    myWeight4[taxa_train4==1] <- presences_weight4
    myWeight4[taxa_train4==0] <- absences_weight4
    
    set.seed(422)
    ANN_model4 = nnet(taxa_train4 ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM,
                      trace=FALSE,
                      data=train4,
                      MaxNWts = 2000,
                      weights = myWeight4,
                      size = i,
                      maxit = j)
    
    # This code is used to predict on the validation dataset and determine final model performance.
    ANN_model.prediction4 = predict(ANN_model4, val4, type = "class")
    ANN_model.probability4 = predict(ANN_model4, val4, type = "raw")
    # we must declare the levels 1 and 0. Otherwise models without any predicted absences will not be included as a level and the confusion matrix will be incomplete.
    final_ANN_model.prediction4<-factor(ANN_model.prediction4, levels = c(0, 1))
    
    # Validation confusion matrix.
    confusion_matrix4=table(taxa_val4, final_ANN_model.prediction4)
    
    # roc.area function requires validation taxa classes (0 or 1) to be numeric.
    numeric_taxa_val4 <- as.numeric(as.character(taxa_val4))
    
    # Performance metrics.
    metric4 <- class.sum(confusion_matrix4)
    Youden4 <- metric4[4,2]
    Youden4 <- as.numeric(Youden4)
    AUC4 <- roc.area(numeric_taxa_val4, ANN_model.probability4)$A
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
    total5 = absences5 + presences5
    
    presences_fraction5 = presences5/total5
    presences_weight5 = 1 - presences_fraction5
    absences_fraction5 = absences5/total5
    absences_weight5 = 1 - absences_fraction5
    
    myWeight5 <- rep(NA,length(taxa_train5))
    myWeight5[taxa_train5==1] <- presences_weight5
    myWeight5[taxa_train5==0] <- absences_weight5
    
    set.seed(422)
    ANN_model5 = nnet(taxa_train5 ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM,
                      trace=FALSE,
                      data=train5,
                      MaxNWts = 2000,
                      weights = myWeight5,
                      size = i,
                      maxit = j)
    
    # This code is used to predict on the validation dataset and determine final model performance.
    ANN_model.prediction5 = predict(ANN_model5, val5, type = "class")
    ANN_model.probability5 = predict(ANN_model5, val5, type = "raw")
    # we must declare the levels 1 and 0. Otherwise models without any predicted absences will not be included as a level and the confusion matrix will be incomplete.
    final_ANN_model.prediction5<-factor(ANN_model.prediction5, levels = c(0, 1))
    
    # Validation confusion matrix.
    confusion_matrix5=table(taxa_val5,final_ANN_model.prediction5)
    
    # roc.area function requires validation taxa classes (0 or 1) to be numeric.
    numeric_taxa_val5 <- as.numeric(as.character(taxa_val5))
    
    # Performance metrics.
    metric5 <- class.sum(confusion_matrix5)
    Youden5 <- metric5[4,2]
    Youden5 <- as.numeric(Youden5)
    AUC5 <- roc.area(numeric_taxa_val5, ANN_model.probability5)$A
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
    
    print(cbind(c("size","maxits"),c(i,j)))
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
