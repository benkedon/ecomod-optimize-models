# Print filename
print(c("File used: ", commandArgs(trailingOnly=FALSE)[4]));
# Pull the species_name from the command line argument

cmd_args = commandArgs(trailingOnly=TRUE);
species_name = paste(cmd_args[1], sep = "");

print(c("SPECIES NAME: ", species_name)); 


# The core gradient boosting package.
library(gbm)
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

interaction.depth = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
n.trees = c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800)
shrinkage = c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1 , 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2 , 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3 , 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4 , 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5 , 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6 , 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7 , 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8 , 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9 , 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0)

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

for (i in interaction.depth){
  for (j in n.trees){
      for (m in shrinkage){
      try({
      #################
      # train 1
      #################
      
      train1 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "1_TRAIN.csv", sep=""));
      val1 = read.csv(paste("datasets/", species_name, "_datasets/", species_name, "1_TEST.csv", sep=""));
      
      # with as.factor the the predict function returns all NaNs. This is unique to GBM.
      taxa_train1 <- as.numeric(as.character(train1[[species_name]]))
      taxa_val1 <- as.numeric(as.character(val1[[species_name]]))
      
      class_instances1<-table(taxa_train1)
      absences1 = class_instances1[1]
      presences1 = class_instances1[2]
      number_presences1 = as.numeric(class_instances1[2])
      total1 = absences1 + presences1
      
      presences_fraction1 = presences1/total1
      presences_weight1 = 1 - presences_fraction1
      absences_fraction1 = absences1/total1
      absences_weight1 = 1 - absences_fraction1
      
      set.seed(422)
      GBM_model1 = gbm(taxa_train1 ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM,
                      distribution="bernoulli",
                      interaction.depth=i,
                      n.trees=j, 
                      shrinkage=m,
                      data=train1)
      
      predicted_presence_probability1 = predict(GBM_model1, n.trees=j, newdata=val1, type="response")
      predicted_class1 <- ifelse(predicted_presence_probability1>=presences_fraction1, 1, 0)
      
      # we must declare the levels 1 and 0. Otherwise models without any predicted absences will not be included as a level and the confusion matrix will be incomplete.
      final_predicted_class1<-factor(predicted_class1, levels = c(0, 1))
      
      # confusion matrix.
      confusion_matrix1 = table(taxa_val1, final_predicted_class1)
      
      # Performance metrics.
      metric1 <- class.sum(confusion_matrix1)
      Youden1 <- metric1[4,2]
      Youden1 <- as.numeric(Youden1)
      AUC1 <- roc.area(taxa_val1, predicted_presence_probability1)$A
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
      
      # with as.factor the the predict function returns all NaNs. This is unique to GBM.
      taxa_train2 <- as.numeric(as.character(train2[[species_name]]))
      taxa_val2 <- as.numeric(as.character(val2[[species_name]]))
      
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
      GBM_model2 = gbm(taxa_train2 ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM,
                      distribution="bernoulli",
                      interaction.depth=i,
                      n.trees=j, 
                      shrinkage=m,
                      data=train2)
      
      predicted_presence_probability2 = predict(GBM_model2, n.trees=j, newdata=val2, type="response")
      predicted_class2 <- ifelse(predicted_presence_probability2>=presences_fraction2, 1, 0)
      
      # we must declare the levels 1 and 0. Otherwise models without any predicted absences will not be included as a level and the confusion matrix will be incomplete.
      final_predicted_class2<-factor(predicted_class2, levels = c(0, 1))
      
      # confusion matrix.
      confusion_matrix2 = table(taxa_val2, final_predicted_class2)
      
      # Performance metrics.
      metric2 <- class.sum(confusion_matrix2)
      Youden2 <- metric2[4,2]
      Youden2 <- as.numeric(Youden2)
      AUC2 <- roc.area(taxa_val2, predicted_presence_probability2)$A
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
      
      # with as.factor the the predict function returns all NaNs. This is unique to GBM.
      taxa_train3 <- as.numeric(as.character(train3[[species_name]]))
      taxa_val3 <- as.numeric(as.character(val3[[species_name]]))
  
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
      GBM_model3 = gbm(taxa_train3 ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM,
                      distribution="bernoulli",
                      interaction.depth=i,
                      n.trees=j, 
                      shrinkage=m,
                      data=train3)
      
      predicted_presence_probability3 = predict(GBM_model3, n.trees=j, newdata=val3, type="response")
      predicted_class3 <- ifelse(predicted_presence_probability3>=presences_fraction3, 1, 0)
      
      # we must declare the levels 1 and 0. Otherwise models without any predicted absences will not be included as a level and the confusion matrix will be incomplete.
      final_predicted_class3<-factor(predicted_class3, levels = c(0, 1))
      
      # confusion matrix.
      confusion_matrix3 = table(taxa_val3, final_predicted_class3)
      
      # Performance metrics.
      metric3 <- class.sum(confusion_matrix3)
      Youden3 <- metric3[4,2]
      Youden3 <- as.numeric(Youden3)
      AUC3 <- roc.area(taxa_val3, predicted_presence_probability3)$A
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
      
      # with as.factor the the predict function returns all NaNs. This is unique to GBM.
      taxa_train4 <- as.numeric(as.character(train4[[species_name]]))
      taxa_val4 <- as.numeric(as.character(val4[[species_name]]))
      
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
      GBM_model4 = gbm(taxa_train4 ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM,
                      distribution="bernoulli",
                      interaction.depth=i,
                      n.trees=j, 
                      shrinkage=m,
                      data=train4)
      
      predicted_presence_probability4 = predict(GBM_model4, n.trees=j, newdata=val4, type="response")
      predicted_class4 <- ifelse(predicted_presence_probability4>=presences_fraction4, 1, 0)
      
      # we must declare the levels 1 and 0. Otherwise models without any predicted absences will not be included as a level and the confusion matrix will be incomplete.
      final_predicted_class4<-factor(predicted_class4, levels = c(0, 1))
      
      # confusion matrix.
      confusion_matrix4 = table(taxa_val4, final_predicted_class4)
      
      # Performance metrics.
      metric4 <- class.sum(confusion_matrix4)
      Youden4 <- metric4[4,2]
      Youden4 <- as.numeric(Youden4)
      AUC4 <- roc.area(taxa_val4, predicted_presence_probability4)$A
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
      
      # with as.factor the the predict function returns all NaNs. This is unique to GBM.
      taxa_train5 <- as.numeric(as.character(train5[[species_name]]))
      taxa_val5 <- as.numeric(as.character(val5[[species_name]]))
      
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
      GBM_model5 = gbm(taxa_train5 ~ MSST+ANC+CA+CL+COND+K+MG+SODIUM+SO4+DOY+LSUB_DMM,
                      distribution="bernoulli",
                      interaction.depth=i,
                      n.trees=j, 
                      shrinkage=m,
                      data=train5)
   
      predicted_presence_probability5 = predict(GBM_model5, n.trees=j, newdata=val5, type="response")
      predicted_class5 <- ifelse(predicted_presence_probability5>=presences_fraction5, 1, 0)
      
      # we must declare the levels 1 and 0. Otherwise models without any predicted absences will not be included as a level and the confusion matrix will be incomplete.
      final_predicted_class5<-factor(predicted_class5, levels = c(0, 1))
      
      # confusion matrix.
      confusion_matrix5 = table(taxa_val5, final_predicted_class5)
      
      # Performance metrics.
      metric5 <- class.sum(confusion_matrix5)
      Youden5 <- metric5[4,2]
      Youden5 <- as.numeric(Youden5)
      AUC5 <- roc.area(taxa_val5, predicted_presence_probability5)$A
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
      
      print(cbind(c("interaction depth","n.trees","shrinkage"),c(i,j,m)))
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
}  

# Prints the mean Youden Index value for each set of hyperparameters (i.e. models).
# Max Youden Index value in this list corresponds to optimal model and set of hyperparameters. 
print("mean_Youden_list");
print(mean_Youden_list)

# Returns the highest Youden Index value from the 5-fold cross validation.
# Then query the raw output text file to find corresponding hyperparameters if desired (or scroll through R-console output).
print("max(unlist(mean_Youden_list))");
max(unlist(mean_Youden_list))
