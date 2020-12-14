# TEEBoT 
## (Tissue Expression Estimation using Blood Transcriptome)

TEEBot is a tool to predict tissue-specific expression (TSGE) from an individual’s blood transcriptome and demographic 
information (age, gender and race), with clinical implications. If his/her genotype information is avalibale, it could boost the performance of TEEBot. 
We trained TEEBot on GTEx version 6 and evaluated its performance in a cross-validation manner. For each gene in each target tissue, we first evaluate its predictability based on LLR test and then fit a lasso regression model to estimate its TSGE. Our models require Whole blood gene expression (WBGE), Whole Blood splicing (WBSp) information, and three demographic ‘confounding’ factors (Age, Race, and Sex), with genetype information as one additonal option. 

## Data download
Download data files from GTEx Portal and dbGaP.

## Main script
The code for prediction of gene expression of a target tissue consists of the following steps:

### Model building 
For each gene the top 5 PCs of Whole blood transcriptome is used to build the model in training set of the data 
and then the model is used to predict its expression in testing set. Five-fold cross validation prediction for 
all the genes are performed. Model is build using lasso regression method using cv.glmnet() function from glmnet R-package.

### Measuring prediction accuracy 
The prediction accuracy for each of the genes are measure using pearson correlation coefficient. Also likelihood ratio test 
is performed for each gene  to assess which gene's expression are predictable from blood transcriptome above and beyond the 
confounders (age, race and sex). Thus for each gene we have its prediatability score (in terms of Pearson correlation coefficient)
and FDR of the LLR p-value. Prediction accuracy for each tissue is reported based on the genes which pass the likelihood ratio 
test (FDR<=0.05) and also have pearson correlation coefficient values beyond a predefined threshold. 


## Code running instructions
### For prediction of target tissue
Code: regression_model_articleoutput.r
#### Input
Input are the gene expression, phenotype and genotype file downloaded from GTEx Portal and dbGaP.
The input parameters to run the code are index, run, fold, ENS. 
 We also have a input file "tissue_pair_work_v6.RData" which contains a table with first column 
 as whole blood, second column target tissue and third colum are the number of common samples. 

1) "index" are the tissue ids, the ids for the 32 target tissues are as follows 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,32,34,37,39. 
2) "run" indicated which model to compute, M1 or M2. To compute model M1 run=100, and for M2 run=110
3) fold : how many fold cross validation we want.
4) ENS : how many independent predictions we want to do.

#### Output
Each Target tissue prediction file for example "Whole Blood_xx_lasso_regress_G_Spl_10_20.Rdata", 
        where xx will be replaced by the target tissue name. 
This Rdata contanis two objects, gene.regress.glmnet and gene.select.LLR.
1) gene.regress.glmnet : is a list, gene.regress.glmnet[["Mean-first"]] contains our desired correlation coefficient 
that is the predictibility scores.
2) gene.select.LLR : is a list with LLR ratio and pvalue

note: Download the code file (regression_model_articleoutput.r) and the function.r in the same folder. Example path to store the
codes is "/home/Desktop/Imputation". Within Imputation make a folder named "input" and "output". In the "output" folder the output of the codes automatically gets saved.
