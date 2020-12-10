# TEEBoT
Pipeline for tissue-specific gene expression prediction from the whole blood transcriptome.
For each tissue a linear model is build using GTEx data. 

## Data download
Download data files from GTEx Portal.

## Main script
The code for prediction of gene expression of a target tissue consists of the following steps:

##Model building 
For each gene the top 5 PCs of Whole blood transcriptome is used to build the model in training set of the data and then the prediction of the gene is performed in testing set. Five-fold cross validation prediction for all the genes are performed. 
Model is build using lasso regression method using cv.glmnet function from glmnet R-package.

##Measuring prediction accuracy 
The prediction accuracy for each of the genes are measure using pearson correlation coefficient. Also likelihood ratio test is performed for each gene  to assess which gene's expression are predictable from blood transcriptome above and beyond the confounders (age, race and sex).
 
Prediction accuracy for each tissue is reported based on the genes which pass the likelihood ratio test (FDR<=0.05) and also beyond some thresholf pearson correlation coefficient values. 



Folder : code : Includes all the code

Folder : input : Include all the input files neede to run the codes

Folder : output : Includes the model files for M1 and M2 
                  Also the output of the codes automatically gets saved in the output folder. 
