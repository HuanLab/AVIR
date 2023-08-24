#2023-04
# Purpose: this code is a demonstration of training AVIR model for predicting the quality of peak integration
## given the peak area and peak height of a metabolic feature.

# Created: 2023-04
# Edited : 2023-05  
# Created by: Zixuan Zhang
# Edited by : Zixuan Zhang
# Version: 0.0.1

####################################################################################
# We utilize the correlation between peak area and peak height to detect inconsistent peak integration pattern
# therefore we can capture the false peak integration (you can consider there exists outlier) in a metabolic feature
## that causes large computational variation in the quantitative analysis.

# We provide the training data, and an example data for running the prediction in R.
# Given the training data, this model can also be generated in other programming language (Python). 

## I also have attached the raw files of peak area and peak height that are used to generate the testing data.
## It is a demo format for calculating the value of machine learning feature for the model prediction.
## The code that generate the testing data is also shown below.

## The detailed instruction be be found in our github, named as "AVIR User_Instruction.pdf"
## The design and discussion of the model can also be found in the user manual.



## First we need to install the required library to run our machine learning model.
## This code can either be uncommented or run in another script to install the necessary packages

## install.packages("e1071")
## install.packages("caret")


## Load the required library
library(e1071) 
library(caret)

## Set a working directory (the folder position) in your computer, this step specifies the location of your files
## Remember to use forward slashes, not backslashes if you copy the directory path
Working_directory <- "C:/User/Users/AVIR_Demo" 
setwd(Working_directory)

## Read the AVIR model
Avir = readRDS("AVIR.rds")

#################################################################
## The following code is an example of how I calculated the value of the Avir feature
## Read the table of metabolic feature using peak area and peak height to represent the intensity respectively
## When running your own data, replace the names of PeakArea and PeakHeight with the names of your area and height files
df_PA = read.csv('PeakArea_Demo.csv')
df_PH = read.csv('PeakHeight_Demo.csv')

## Create a dataframe to store the prediction results of Avir model
df_result <- as.data.frame(df_PA[,1:3])
df_result$Prediction <- NA

## Create a dataframe to store the values of Avir model's machine learning feature
df_predict<- as.data.frame(df_result[,1])

df_predict$Spearman_Cor <- NA
df_predict$Pearson_Cor <- NA
df_predict$RSD_PAPH <- NA
df_predict$norm_diff_PA_PH_median <- NA

## Here is the example code for calculation of features of Avir, you may need to modify by case. 
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

#############################
## Here is the filtering process, we set a threshold to remove those intensity (peak height) close to noise
## For the Bruker Impact II QTOF in the Huan lab, I used 500 as noise cutoff. Intensity above 1000 is an acceptable signal
## Modify this section's cutoffs (for both signal and percent of samples needed) as needed for your mass spectrometer
  
  ## n1 represent intensity filter
  n1 = sum(df_feature[,2] > 1000)
  
  ## Add one more noise filter, filter the feature that contains noise, low-intensity feature
  n2 = sum(df_feature[,2] <= 1000)
  
## Here is how I specifically filter the low-quality metabolic feature by noise.
## For a feature with an intensity of 1000 in over 20% of samples and no samples with noise-like intensity (< 500)
## You should set a custom cutoff for your noise level based on your mass spectrometer
  if( n1 > nrow(df_feature)*0.2 & n2 == 0  ) {
    
    df_feature= df_feature[apply(df_feature!=0, 1, all),]
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
  
  else if( n1 > nrow(df_feature)*0.2 & n2 > 0  ){

## For samples which have features which have noise-like intensity in some of the samples, 
## but more than 20% of samples show intensity over the acceptable threshold for the feature,
## extreme coefficients and values are assigned to the features of Avir so that you can recognize these metabolic features easily
    df_predict$Spearman_Cor[i] = -1
    df_predict$Pearson_Cor[i] = -1
    df_predict$RSD_PAPH[i] = 200
    df_predict$norm_diff_PA_PH_median[i] = 200
    
  }
  
  else{

## For those metabolic features of low-reproducibility (intensity < 500 in 80% of samples),
## I give different extreme values for this kind (suggests variation in the feature due to misalignment or misintegration -- computational variation)  
    df_predict$Spearman_Cor[i] = 0
    df_predict$Pearson_Cor[i] = 0
    df_predict$RSD_PAPH[i] = 100
    df_predict$norm_diff_PA_PH_median[i] = 100
    
  }
  
}

## Then we predict the integration quality for metabolic features.
Prediction <- predict(Avir, newdata = df_predict[-1])

# Add the prediction result to the dataframe that we stored the results, 1 represent TRUE, 0 represent FALSE.
# i.e., if assigned 1, it is a feature with low variation among the samples (proper alignment / integration)
# if assigned 0, it is either misaligned or misintegrated; there is computational variation present in the feature
df_result$Prediction <- Prediction
  
# Add the input value of Avir to the dataframe, so that we can see the category of the quality of the metabolic features.
# Shows the Spearman, Pearson coefficients as well as the area-to-height ratio
df_result <- cbind(df_result, df_predict[,2:5])

## Output the prediction result in the working directory.
write.csv(df_result, file = "AVIR_PredictionResult.csv", row.names = FALSE)
