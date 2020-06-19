#code to find predicted expression in target tissue


regression_glmnet_gene_splicing_M2<-function(nPC1,nPC2,temp.tss1,temp.snps,temp.tss2,genes,method,confund.x){

	#nPC1=PC1; nPC2=PC2;temp.tss1=gtex.tss1.pc; temp.snps=splice.mat.pc; temp.tss2=gtex.tss2;genes=gene[1:10];method="lasso";confund.x=confund;

require(foreach)
require(glmnet)
#library(doMC)
#registerDoMC(cores=25)
if(method=="lasso"){index=1}
if(method=="ridge"){index=0}

			gtex.train.tss2=as.matrix(temp.tss2); 
			gtex.train.tss1.pc=as.matrix(rbind(temp.tss1[1:nPC1,],temp.snps[1:nPC2,],t(confund.x)));

			eachgene=foreach(itr=1:length(genes),.inorder=T,.combine='c') %do% {
				geneid=itr;
				fit.cv=cv.glmnet(t(gtex.train.tss1.pc),gtex.train.tss2[geneid,], alpha=index,nfolds=5)
				lam=fit.cv$lambda.min
				fit1 = glmnet(t(gtex.train.tss1.pc),gtex.train.tss2[geneid,], alpha=index, lambda=lam)
				#pred=predict(fit1,newx=t(gtex.test.tss1.pc))
				
				temp=list(fit1)
				return(temp)
				} #gene loop ends
				names(eachgene)=genes
return(eachgene)

}

regression_glmnet_gene_M1<-function(nPC,temp.tss1,temp.tss2,genes,method,confund.x){
require(foreach)
require(glmnet)
#library(doMC)
#registerDoMC(cores=25)
if(method=="lasso"){index=1}
if(method=="ridge"){index=0}

			gtex.train.tss2=as.matrix(temp.tss2); 
			gtex.train.tss1.pc=as.matrix(rbind(temp.tss1[1:nPC,],t(confund.x)));

			eachgene=foreach(itr=1:length(genes),.inorder=T,.combine='c') %do% {
				geneid=itr;
				fit.cv=cv.glmnet(t(gtex.train.tss1.pc),gtex.train.tss2[geneid,], alpha=index,nfolds=5)
				lam=fit.cv$lambda.min
				fit1 = glmnet(t(gtex.train.tss1.pc),gtex.train.tss2[geneid,], alpha=index, lambda=lam)
				#pred=predict(fit1,newx=t(gtex.test.tss1.pc))
				
				temp=list(fit1)
				return(temp)
				} #gene loop ends
				names(eachgene)=genes
return(eachgene)

}





dir="/home/Desktop" #example of dir that needs to be given 
workdir=paste0(dir,"/imputation"); 

setwd(workdir)


library(data.table)
args<-commandArgs(TRUE);


index=as.numeric(args[1]) #target tissue
run=toString(args[2]) #category
#fold=as.numeric(args[3]); #fold
#ENS=as.numeric(args[4]); #ENS

#PC1=as.numeric(args[3]); #GE
#PC2=as.numeric(args[4]); #Sp
#PC3=as.numeric(args[5]); #Snp
PC1=10; PC2=20; 

#workdir="/Volumes/5TBbackup/UMD2019/project2_scratch_mbasu/gtex_v6/blood_cross_talk/gtex-v6"



#loading functions
#source('/Volumes/5TBbackup/UMD2019/project2_scratch_mbasu/gtex_v6/blood_cross_talk/gtex-v6/function_regression.r')
#source(paste0(workdir,'/code/function.r'));

#load gtex data 
#load("/Volumes/5TBbackup/UMD2019/project2_scratch_mbasu/gtex_v6/data/gtexdata_v6.RData")
load(paste0(workdir,"/input/gtexdata_v6.RData"));

#Common samples between blood and target tissue
#load("/Volumes/5TBbackup/UMD2019/project2_scratch_mbasu/gtex_v6/blood_cross_talk/tissue_pair_work_v6.RData")
load(paste0(workdir,"/input/tissue_pair_work_v6.RData"));

#loading patient ids that have snps information
#load("/Volumes/5TBbackup/UMD2019/project2_scratch_mbasu/gtex_v6/blood_cross_talk/gtex-v6/PC_SNPs/snps_sample_v6.RData")
load(paste0(workdir,"/input/snps_sample_v6.RData"));


gene=gtex.pc$Name
tss1="Whole Blood"
for(index in c(1,2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 32, 34, 37, 39)){
tss2=as.character(unlist(strsplit(tss_pair[index,],"Whole Blood"))[2]); #tss2: target tissue
allgeneid=as.numeric(seq(length(gene)))
print(tss2);
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
outdir=paste0(workdir,"/output/Model/M2")
if(!file.exists(outdir))dir.create(outdir, recursive = TRUE)
str=sprintf("%s/%s_lasso_regress_G_Spl_%d_%d.Rdata",outdir,tss2,PC1,PC2)
ptm <- proc.time()

modelM2=regression_glmnet_gene_splicing_M2(PC1,PC2,gtex.tss1.pc,splice.mat.pc,gtex.tss2,gene,"lasso",confund)

save(modelM2,file=str)
print(proc.time() - ptm)
#------------------------------------------------------------
}


if(run=="100"){
#Gene  ------------------------------------------------------
outdir=paste0(workdir,"/output/Model/M1")
if(!file.exists(outdir))dir.create(outdir, recursive = TRUE)
str=sprintf("%s/%s_lasso_regress_G_%d.Rdata",outdir,tss2,PC1)

ptm <- proc.time()
modelM1=regression_glmnet_gene_M1(PC1,gtex.tss1.pc,gtex.tss2,gene,"lasso",confund)
save(modelM1,file=str)

print(proc.time() - ptm)
#------------------------------------------------------------
}
}









