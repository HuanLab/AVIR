## Load the required library
library(e1071) 
library(caret)

## Set a working directory (the folder postion) in your computer, this step specifies the location of your files
Working_directory <- "C:/Users/User/Desktop/Avir.Demo_2.0" 
setwd(Working_directory)

#################################################################
## The following code is an example of how I calculated the value of the Avir feature
## Read the table of metabolic feature using peak area and peak height to represent the intensity respectively
df_PA = read.csv('PeakArea_Demo.csv')
df_PH = read.csv('PeakHeight_Demo.csv')
df_Label = read.csv("Label.csv")

## Create a dataframe to store the prediction results of SVM model
df_result <- as.data.frame(df_PA[,1:3])
df_result$Prediction <- NA

## Create a dataframe to store the values of SVM model's machine learning feature
df_predict <- as.data.frame(df_result[, 1])
colnames(df_predict)[1] <- "Alignment ID"
df_predict$Spearman_Cor <- NA
df_predict$Pearson_Cor <- NA
df_predict$RSD_PAPH <- NA
df_predict$norm_diff_PA_PH_median <- NA

## Here is the example code for calculation of features of SVM, you may need to modify by case. 
## If using my given format then you can just use this.
## The purpose of each step is shown below.

for (i in 1:nrow(df_PA) ) {
  
  ## In the For loop, we first create a data frame to store the peak area and peak height of every metabolic feature.
  ## You may need to change this line of code if you use different formats as input.
  
  df_feature <- rbind(as.numeric(df_PA[i, 4:ncol(df_PA) ]), as.numeric(df_PH[i, 4:ncol(df_PA) ]))
  
  ##############################################################################  
  ## You don't need to change following code, it is just formatting
  df_feature  <- t(df_feature)
  colnames(df_feature ) = c("PA", "PH")
  df_feature[,1] = as.numeric(df_feature[,1])
  df_feature[,2] = as.numeric(df_feature[,2])
  df_feature  = as.data.frame(df_feature )
  
  #########################################################
  
    ## Calculate the Spearman correlation and Pearson correlation
    spearman_cor1 = cor(y = df_feature[,1], x = df_feature[,2], method = "spearman")
    pearson_cor1 = cor(y = df_feature[,1], x = df_feature[,2], method = "pearson")
    
    df_predict$Spearman_Cor[i] = spearman_cor1
    df_predict$Pearson_Cor[i] = pearson_cor1
    
    # Calculate the ratio of peak area to peak height first, then calculate the RSD of PA/PH and range-median ratio of PA/PH
    Ratio1 = df_feature[,1] / df_feature[,2]
    mean1 = mean(Ratio1)
    sd1= sd(Ratio1)
    rsd1 = sd1/mean1
    
    df_predict$RSD_PAPH[i] = rsd1
    df_predict$norm_diff_PA_PH_median[i] = ( max(Ratio1) - min(Ratio1) ) / median(Ratio1)
    
 
}

Trainingdata <- cbind(df_predict[,-1], df_Label[,2])
colnames(Trainingdata)[5] <- "Label"

## Train the SVM model
Classifier = svm(formula = Label ~ .,
           data = Trainingdata,
           type = 'C-classification',
           kernel = 'linear',
           probability=TRUE)

saveRDS(Classifier, file = "SVM.rds")

