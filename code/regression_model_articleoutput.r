#code to find predicted expression in target tissue

dir="/home/Desktop" #example of dir that needs to be given 
workdir=paste0(dir,"/imputation"); 

setwd(workdir)


library(data.table)
args<-commandArgs(TRUE);


index=as.numeric(args[1]) #target tissue
run=toString(args[2]) #category
fold=as.numeric(args[3]); #fold
ENS=as.numeric(args[4]); #ENS

PC1=as.numeric(args[3]); #GE
PC2=as.numeric(args[4]); #Sp
PC3=as.numeric(args[5]); #Snp
#PC1=10; PC2=20; 

#workdir="/Volumes/5TBbackup/UMD2019/project2_scratch_mbasu/gtex_v6/blood_cross_talk/gtex-v6"



#loading functions
#source('/Volumes/5TBbackup/UMD2019/project2_scratch_mbasu/gtex_v6/blood_cross_talk/gtex-v6/function_regression.r')
source(paste0(workdir,'/code/function.r'));

#load gtex data 
format_gtexdatav6(); 
#In this function we directly read gtex v6 phenotype data and gene expression data and merge them 
#into format convenient to be used in the code and create a gtexdata_v6.RData file in workdir/input/ folder.
load(paste0(workdir,"/input/gtexdata_v6.RData"));

#Common samples between blood and target tissue
#load("/Volumes/5TBbackup/UMD2019/project2_scratch_mbasu/gtex_v6/blood_cross_talk/tissue_pair_work_v6.RData")
load(paste0(workdir,"/input/tissue_pair_work_v6.RData"));

#finding patient ids that have snps information
load(paste0(workdir,"/input/genotype.v6.chr1_22.Rdata") #loading genotype data obatined from dbGaP
pat.snps=sapply(strsplit(colnames(genotype.v6[["chr2"]]), split="-"), '[[', 2)


gene=gtex.pc$Name
tss1="Whole Blood"
tss2=as.character(unlist(strsplit(tss_pair[index,],"Whole Blood"))[2]); #tss2: target tissue
allgeneid=as.numeric(seq(length(gene)))

comm.patients=patients[which(patients %in% pat.snps)]

#for common individuals seperate the gene expression for tss1 and tss2 
j1=c(); j2=c();
pat.com=c(); #common individuals between Blood (tss1) and target tissue (tss2)
for (ipat in comm.patients){
	i1=which(expc.nt$patient==ipat & expc.nt$SMTSD==tss1);
	i2=which(expc.nt$patient==ipat & expc.nt$SMTSD==tss2);
	if (length(i1)>0 & length(i2)>0){ j1=append(j1,i1[1]); j2=append(j2,i2[1]); pat.com=append(pat.com,ipat);}}
nn1=length(j1); nn2=length(j2); it1=matrix(0, 1, nn1); it2=matrix(0, 1, nn2);
it1[1,]=j1; it2[1,]=j2;
gtex.tss1=gtex.pc[,c(it1[1,]), with=F]; 
gtex.tss2=gtex.pc[,c(it2[1,]), with=F]; 


#PCA for tss1, whole Blood gtex
mdsk=100;
d = dist(t(gtex.tss1), method="euclidean"); mdsk=dim(t(gtex.tss1))[1]-1;
pcfit = cmdscale(d, eig=TRUE, k=mdsk);
gtex.tss1.pc=t(pcfit$points); #row:PCs; col: samples
colnames(gtex.tss1.pc)=pat.com;

#loading splicing PCs
#load('/Volumes/5TBbackup/UMD2019/project2_scratch_mbasu/gtex_v6/blood_cross_talk/gtex-v6/PC_splicing/PCcmd_splicing_blood.Rdata' )
load(paste0(workdir,"/input/PCcmd_splicing_blood.Rdata"))
splice.mat.pc=splicing.pc[[tss2]]

#confounders
x=sapply(1:length(pat.com),function(i){which(pheno.subj %in% pat.com[i])})
age=as.numeric(pheno.dt$AGE[x]);
race=as.numeric(pheno.dt$RACE[x]); race[race==99]=0;
gender=as.numeric(pheno.dt$GENDER[x]);
confund=cbind(age,race,gender)


 #row:PCs; col: samples

if(run=="110"){
#Gene + Splicing --------------------------------------------
ptm <- proc.time()

#outdir=paste0(workdir,"/output")
#outdir="/cbcb/project2-scratch/mbasu/gtex_v6/blood_cross_talk/gtex-v6/regress_result/lasso_geneexpr_snps_splicing_confound/data/lasso_cv_lambda/G_Sp_CF"
outdir=paste0(workdir,"/output/G_Sp_CF")
if(!file.exists(outdir))dir.create(outdir)
str=sprintf("%s/%s_%s_lasso_regress_G_Spl_%d_%d.Rdata",outdir,tss1,tss2,PC1,PC2)

gene.regress.glmnet=regression_glmnet_gene_splicing(fold,ENS,PC1,PC2,gtex.tss1.pc,splice.mat.pc,gtex.tss2,allgeneid,pat.com,"lasso",confund)

gene.select.LLR=regression_lm_loglik_selectgene_2(PC1,PC2,gtex.tss1.pc,splice.mat.pc,gtex.tss2,"llr",confund)
gene.select.LLR.Sp=regression_lm_loglik_selectgene_2Sp(PC1,PC2,gtex.tss1.pc,splice.mat.pc,gtex.tss2,"llr",confund)
gene.select.LLR.GE=regression_lm_loglik_selectgene_2GE(PC1,PC2,gtex.tss1.pc,splice.mat.pc,gtex.tss2,"llr",confund)

save(gene.regress.glmnet,gene.select.LLR,gene.select.LLR.Sp,gene.select.LLR.GE,file=str)
print(proc.time() - ptm)
#------------------------------------------------------------
}


if(run=="100"){
#Gene  ------------------------------------------------------
#outdir="/cbcb/project2-scratch/mbasu/gtex_v6/blood_cross_talk/gtex-v6/regress_result/lasso_geneexpr_snps_splicing_confound/data/lasso_cv_lambda/G_CF"
outdir=paste0(workdir,"/output/G_CF")
if(!file.exists(outdir))dir.create(outdir)
str=sprintf("%s/%s_%s_lasso_regress_G_%d_%dfold_%dENS.Rdata",myfolder,tss1,tss2,PC1,fold,ENS)

ptm <- proc.time()
gene.regress.glmnet=regression_glmnet_gene(fold,ENS,PC1,gtex.tss1.pc,gtex.tss2,allgeneid,pat.com,"lasso",confund)

gene.select.LLR=regression_lm_loglik_selectgene_1(PC1,gtex.tss1.pc,gtex.tss2,"llr",confund)
save(gene.regress.glmnet,gene.select.LLR,file=str)
print(proc.time() - ptm)
#------------------------------------------------------------
}








