format_gtexdatav6<-function(dirx){
library(data.table)
dt=fread("phs000424.v6.pht002742.v6.p1.c1.GTEx_Subject_Phenotypes.GRU.txt") #phenotype data 
pheno.dt=dt[-1,]
setnames(pheno.dt,as.character(dt[1,]))
#subj=unlist(lapply(phe.dt$SUBJID, function(tt) {substr(tt,6,9)}))
pheno.subj=sapply(1:nrow(pheno.dt), function(x){ y=unlist(strsplit(pheno.dt$SUBJID[x], split="-")); l=length(y); y[l] })

#gene expression
rpkm=fread("All_Tissue_Site_Details_Analysis.combined.rpkm.gct") #gene expression data
rpkm$Name = sapply(strsplit(rpkm$Name, '[.]'), '[[', 1)

rpkm1=unique(rpkm)

rpkm.col = colnames(rpkm); #gtex.row = rownames(gtex);
patients = unique(sapply(strsplit(rpkm.col[-(1:2)], split="-"), '[[', 2))


annt=fread("GTEx_Data_V6_Annotations_SampleAttributesDS.txt")
tmp = data.table(SAMPID=rpkm.col[-(1:2)], myid = 1:(length(rpkm.col) -2)); 
setkey(tmp, SAMPID); 
setkey(annt, SAMPID);
annt = merge(x=tmp, y=annt, by="SAMPID"); 
annt$patient = sapply(strsplit(annt$SAMPID, split="-"), '[[', 2)
annt1=annt[order(annt$myid),]
expc.nt=annt1
#========================================================================

#=======================================================================
#### gene information to identify protein coding genes
require(pracma)
gene=fread("Homo_sapiens.gene_info") 
gene.symbol=gene$V3

gene.IDs=gene$V6 
ensmbl=array("", c(length(gene.IDs),1))
for (i in seq(length(gene.IDs))){
	gIDs=gene.IDs[i]
	tt=strfind(gIDs, 'Ensembl:', overlap = TRUE)
	if (!is.null(tt)){
		ensmbl[i]=substr(gIDs, tt+8, tt+22)
	}
}

gene[, c("Name") := ensmbl]

tmp = data.table(Name=gene$Name, myid = 1:dim(gene)[1]); setkey(tmp, Name); 
newrpkm = merge(x=tmp, y=rpkm, by="Name") 

tmp = data.table(Name=rpkm$Name, myid = 1:dim(rpkm)[1]); setkey(tmp, Name);
newgene = merge(x=tmp, y=gene, by="Name") #contains the newgene with only those which are in gtex
  
newrpkm$myid = newgene$myid; 
#the order of genes in newrpkm and newgene are same, checked using which(newrpkm$myid!=newgene$myid)
newrpkm$geneSym = newgene$V3 #gene name
newrpkm$geneFun = newgene$V9; #gene description
newrpkm$geneKnd = newgene$V10; #gene type
newrpkm$chromosome = newgene$V7; #chromosome

newrpkm.pc = newrpkm[!is.na(newrpkm$geneKnd) & newrpkm$geneKnd=="protein-coding" ] #only protein coding genes
g=newrpkm.pc$Name
x=as.data.table(table(g))
xg=x$g[which(x$N>1)]
remove=unlist(sapply(1:length(xg),function(i){r=which(g %in% xg[i]); r[2:length(r)] }))
newrpkm.pc1=newrpkm.pc[-c(remove),];

newrpkm.pc2=newrpkm.pc
newrpkm.pc=newrpkm.pc1
l=dim(newrpkm.pc)[2]; 
gtex.pc1=newrpkm.pc[,c(-1,-2,-3,-(l-3),-(l-2),-(l-1),-l),with=F]
zero.l=sum(gtex.pc1==0)
extr.gtex.pc=newrpkm.pc[,c(1,2,3,(l-3),(l-2),(l-1),l),with=F]
gtex.pc=cbind(gtex.pc1,extr.gtex.pc)
save(gtex.pc,expc.nt,pheno.dt,pheno.subj,patients,file=paste0(dirx,"/input/gtexdata_v6.RData"))
}



regression_glmnet_gene<-function(nfold,ENS,nPC,temp.tss1,temp.tss2,genes,pat.com,method,confund.x){
require(foreach)
require(glmnet)
if(method=="lasso"){index=1}
if(method=="ridge"){index=0}

#library(doMC)
#registerDoMC(cores=25)

	k=round(length(pat.com)/nfold-.5);
	folds=rep(c(1:nfold),c(rep(k,(nfold-1)),(length(pat.com)-(nfold-1)*k)))
	avg.pred=matrix(0,length(pat.com),length(genes))
	corr.perENS=NULL; pval.perENS=NULL;
	sample.matrix=matrix(0,ENS+1,length(pat.com))
	colnames(sample.matrix)=pat.com; 
	sample.matrix[ENS+1,]=folds
	lambda.matrix=matrix(0,length(genes),ENS*nfold)
	jj=0;
for(ee in seq(ENS)){
	sam=sample(1:length(pat.com),length(pat.com),replace=F); sample.matrix[ee,]=sam;
	temp.pred=matrix(0,length(pat.com),length(genes))	
	print(ee)

for(ifold in seq(nfold)){ print(paste0("e",ee,"f",ifold)); jj=jj+1;
	test=sam[which(folds==ifold)]
	train=sam[which(folds!=ifold)]
	gtex.train.tss2=as.matrix(temp.tss2[,train,with=F]); gtex.test.tss2=as.matrix(temp.tss2[,test,with=F]);
	gtex.train.tss1.pc=as.matrix(rbind(temp.tss1[1:nPC,train],t(confund.x[train,]))); 
        gtex.test.tss1.pc =as.matrix(rbind(temp.tss1[1:nPC,test],t(confund.x[test,])));


#eachgene=foreach(itr=1:length(genes),.inorder=T,.combine='cbind') %dopar% {
	eachgene=foreach(itr=1:length(genes),.inorder=T,.combine='cbind') %do% {

		geneid=genes[itr];
		fit.cv=cv.glmnet(t(gtex.train.tss1.pc),gtex.train.tss2[geneid,], alpha=index,nfolds=5); 
		lam=fit.cv$lambda.min; 
		fit1 = glmnet(t(gtex.train.tss1.pc),gtex.train.tss2[geneid,], alpha=index, lambda=lam)
		pred=predict(fit1,newx=t(gtex.test.tss1.pc))
		return(c(pred,lam))
		} #gene loop ends
temp.pred[test,]=eachgene[1:length(test),]; lambda.matrix[,jj]=eachgene[length(test)+1,];
} #fold loop ends

avg.pred=avg.pred+temp.pred; 

temp.tss22=temp.tss2[genes,]
#temp=foreach(itr=1:nrow(temp.tss22),.inorder=T,.combine='rbind') %dopar% {
	temp=foreach(itr=1:nrow(temp.tss22),.inorder=T,.combine='rbind') %do% {

x=cor.test(as.numeric(temp.tss22[itr,]),temp.pred[,itr])
c(x$estimate,x$p.value) }

corr.perENS=cbind(corr.perENS,temp[,1])
pval.perENS=cbind(pval.perENS,temp[,2])

} #ENS  loop ends
avg.pred=avg.pred/ENS

temp.tss22=temp.tss2[genes,]
#corr=foreach(itr=1:nrow(temp.tss22),.inorder=T,.combine='rbind') %dopar% {
	corr=foreach(itr=1:nrow(temp.tss22),.inorder=T,.combine='rbind') %do% {

x=cor.test(as.numeric(temp.tss22[itr,]),avg.pred[,itr])
c(x$estimate,x$p.value)
}

rownames(avg.pred)=pat.com
colnames(corr)=c("cor","pval");

result=list()
result[["Mean-first"]]=corr
result[["Mean-end"]][["cor"]]=corr.perENS
result[["Mean-end"]][["pval"]]=corr.perENS
result[["prediction"]]=avg.pred
result[["sampling"]]=sample.matrix
result[["lambda"]]=lambda.matrix
return(result)

}


regression_glmnet_gene_splicing<-function(nfold,ENS,nPC1,nPC2,temp.tss1,temp.snps,temp.tss2,genes,pat.com,method,confund.x){
require(foreach)
require(glmnet)
library(doMC)
registerDoMC(cores=25)
if(method=="lasso"){index=1}
if(method=="ridge"){index=0}


	k=round(length(pat.com)/nfold-.5);
	folds=rep(c(1:nfold),c(rep(k,(nfold-1)),(length(pat.com)-(nfold-1)*k)))
	avg.pred=matrix(0,length(pat.com),length(genes))
	corr.perENS=NULL; pval.perENS=NULL;
	sample.matrix=matrix(0,ENS+1,length(pat.com))
	colnames(sample.matrix)=pat.com; 
	sample.matrix[ENS+1,]=folds
	lambda.matrix=matrix(0,length(genes),ENS*nfold)
	jj=0;

for(ee in seq(ENS)){
	sam=sample(1:length(pat.com),length(pat.com),replace=F);sample.matrix[ee,]=sam;
	temp.pred=matrix(0,length(pat.com),length(genes))	
	#print(ee)

		for(ifold in seq(nfold)){ print(paste0("e",ee,"f",ifold));
	    	jj=jj+1;
			test=sam[which(folds==ifold)];
			train=sam[which(folds!=ifold)];
			gtex.train.tss2=as.matrix(temp.tss2[,train,with=F]); gtex.test.tss2=as.matrix(temp.tss2[,test,with=F]);
			gtex.train.tss1.pc=as.matrix(rbind(temp.tss1[1:nPC1,train],temp.snps[1:nPC2,train],t(confund.x[train,]))); 
	        gtex.test.tss1.pc =as.matrix(rbind(temp.tss1[1:nPC1,test],temp.snps[1:nPC2,test],t(confund.x[test,])));
			#print("step1")

			eachgene=foreach(itr=1:length(genes),.inorder=T,.combine='cbind') %dopar% {
				geneid=genes[itr];
				fit.cv=cv.glmnet(t(gtex.train.tss1.pc),gtex.train.tss2[geneid,], alpha=index,nfolds=5)
				lam=fit.cv$lambda.min
				fit1 = glmnet(t(gtex.train.tss1.pc),gtex.train.tss2[geneid,], alpha=index, lambda=lam)
				pred=predict(fit1,newx=t(gtex.test.tss1.pc))
				return(c(pred,lam))
				} #gene loop ends
		temp.pred[test,]=eachgene[1:length(test),]; lambda.matrix[,jj]=eachgene[length(test)+1,];
		} #fold loop ends

	avg.pred=avg.pred+temp.pred; 
	temp.tss22=temp.tss2[genes,];
	temp=foreach(itr=1:nrow(temp.tss22),.inorder=T,.combine='rbind') %dopar% {
	x=cor.test(as.numeric(temp.tss22[itr,]),temp.pred[,itr]);
	c(x$estimate,x$p.value) }

	corr.perENS=cbind(corr.perENS,temp[,1])
	pval.perENS=cbind(pval.perENS,temp[,2])

} #ENS  loop ends

avg.pred=avg.pred/ENS
temp.tss22=temp.tss2[genes,]
corr=foreach(itr=1:nrow(temp.tss22),.inorder=T,.combine='rbind') %dopar% {
x=cor.test(as.numeric(temp.tss22[itr,]),avg.pred[,itr])
c(x$estimate,x$p.value)
}

rownames(avg.pred)=pat.com
colnames(corr)=c("cor","pval");

result=list()
result[["Mean-first"]]=corr
result[["Mean-end"]][["cor"]]=corr.perENS
result[["Mean-end"]][["pval"]]=corr.perENS
result[["prediction"]]=avg.pred
result[["sampling"]]=sample.matrix
result[["lambda"]]=lambda.matrix
return(result)

} #function ends


#using loglikelihood of lm test  (GE+Sp+CF~CF): contribution from GE+Sp
regression_lm_loglik_selectgene_2<-function(nPC1,nPC2,temp.tss1,temp.snps,temp.tss2,method,confund.x)
{
require('lmtest')
xH1=as.matrix(cbind(t(temp.tss1[1:nPC1,]),t(temp.snps[1:nPC2,]),confund.x)); xH0=as.matrix(confund.x);
y=as.matrix(t(temp.tss2));

LR=list();pval=c();

for(itr in seq(ncol(y))){
	fit1=lm(y[,itr]~xH1); fit0=lm(y[,itr]~xH0); LR[[itr]]=lrtest(fit1,fit0);
	pval=append(pval,LR[[itr]]$"Pr(>Chisq)"[2]);
	}

loglike.Sgenes=list();
loglike.Sgenes[["loglikelihood-ratio"]]=LR
loglike.Sgenes[["pval"]]=pval
return(loglike.Sgenes)
}

#using loglikelihood of lm test contribution from splicing (GE+Sp+CF~GE+CF): contribution from Sp
regression_lm_loglik_selectgene_2Sp<-function(nPC1,nPC2,temp.tss1,temp.snps,temp.tss2,method,confund.x)
{
require('lmtest')
xH1=as.matrix(cbind(t(temp.tss1[1:nPC1,]),t(temp.snps[1:nPC2,]),confund.x)); 
xH0=as.matrix(cbind(t(temp.tss1[1:nPC1,]),confund.x));
y=as.matrix(t(temp.tss2));

LR=list();pval=c();

for(itr in seq(ncol(y))){
	fit1=lm(y[,itr]~xH1); fit0=lm(y[,itr]~xH0); LR[[itr]]=lrtest(fit1,fit0);
	pval=append(pval,LR[[itr]]$"Pr(>Chisq)"[2]);
	}

loglike.Sgenes=list();
loglike.Sgenes[["loglikelihood-ratio"]]=LR
loglike.Sgenes[["pval"]]=pval
return(loglike.Sgenes)
}


#using loglikelihood of lm test contribution from gene (GE+Sp+CF~Sp+CF): contribution from GE
regression_lm_loglik_selectgene_2GE<-function(nPC1,nPC2,temp.tss1,temp.snps,temp.tss2,method,confund.x)
{
require('lmtest')
xH1=as.matrix(cbind(t(temp.tss1[1:nPC1,]),t(temp.snps[1:nPC2,]),confund.x)); 
xH0=as.matrix(cbind(t(temp.snps[1:nPC2,]),confund.x));
y=as.matrix(t(temp.tss2));

LR=list();pval=c();
  
for(itr in seq(ncol(y))){
	fit1=lm(y[,itr]~xH1); fit0=lm(y[,itr]~xH0); LR[[itr]]=lrtest(fit1,fit0);
	pval=append(pval,LR[[itr]]$"Pr(>Chisq)"[2]);
	}

loglike.Sgenes=list();
loglike.Sgenes[["loglikelihood-ratio"]]=LR
loglike.Sgenes[["pval"]]=pval
return(loglike.Sgenes)
}


#using loglikelihood of lm test  (GE/Sp/SNP+CF~CF): contribution from GE/Sp/SNP
regression_lm_loglik_selectgene_1<-function(nPC,temp.tss1,temp.tss2,method,confund.x)
{
require('lmtest')
xH1=as.matrix(cbind(t(temp.tss1[1:nPC,]),confund.x)); xH0=as.matrix(confund.x);
y=as.matrix(t(temp.tss2));

LR=list();pval=c();

for(itr in seq(ncol(y))){
	fit1=lm(y[,itr]~xH1); fit0=lm(y[,itr]~xH0); LR[[itr]]=lrtest(fit1,fit0);
	#g=append(g,((LR[[itr]]$"Pr(>Chisq)"[2]<=0.05)*1));
	pval=append(pval,LR[[itr]]$"Pr(>Chisq)"[2]);
	}

loglike.Sgenes=list();
loglike.Sgenes[["loglikelihood-ratio"]]=LR
loglike.Sgenes[["pval"]]=pval
return(loglike.Sgenes)
}



prediction_dis_CVfeature<-function(expr.pre,expr.org,blood.expr,istat.x,nfold,ENS,geneid.x){

library(data.table) #cntrl.ox,cntrl.op,cntrl.ob,
require(ROCR)
require(glmnet)
require(foreach)
require('lmtest')

		cvlasso<-function(xinput.train.x,xinput.test.x,yinput.train.x,yinput.test.x){
			fit.cv=cv.glmnet(xinput.train.x,yinput.train.x, alpha=1, nfolds=4);
			lam=fit.cv$lambda.min;
			fit = glmnet(xinput.train.x,yinput.train.x, alpha=1, lambda=lam);
			pre=predict(fit,newx=xinput.test.x);
			pred<-prediction(pre[,1],yinput.test.x);aucval<-performance(pred,"auc");
		return(aucval@y.values[[1]])
		}


		DEgene_disease_wilcox<-function(expr.x,istat.xx,confound.x,genes){

		library(data.table)
		case=which(istat.xx==1);
		cntrl=which(istat.xx==-1);
		pval=c();
		for(g in seq(length(genes))){
			geneid.x=genes[g];
			pval=append(pval,wilcox.test(as.numeric(expr.x[geneid.x,case]),as.numeric(expr.x[geneid.x,cntrl]))$"p.value");}
		return(pval)
		}


		DEgene_disease<-function(expr.x,istat.xx,confound.x,genes){

		library(data.table)

		confound.x=t(confound.x)

		pval=foreach(itr=1:length(genes),.inorder=T,.combine='cbind') %do% { #par
				geneid.x=genes[itr];
				yinput=istat.xx
				xinput=t(as.matrix(rbind(expr.x[geneid.x,],confound.x))); xinput0=t(as.matrix(confound.x)); 
				colnames(xinput)=c("expr",colnames(xinput)[-c(1)]);
				H1=lm(yinput~xinput); H0=lm(yinput~xinput0); LR=lrtest(H1,H0); 
				#pval=append(pval,LR[[itr]]$"Pr(>Chisq)"[2]);
				return(LR$"Pr(>Chisq)"[2])
				}
		return(pval)
		}



result=list();

expr.pre=as.matrix(expr.pre); expr.org=as.matrix(expr.org); blood.expr=as.matrix(blood.expr);

pat.com.x=colnames(expr.pre)

ind3_d=which(istat.x==1);ind3_n=which(istat.x==-1);

x.pred=c(); x.raw=c();x.blood=c(); 

	#eachens.x=foreach(itr=1:ENS,.inorder=T,.combine='cbind') %dopar% {
	eachens.x=foreach(itr=1:ENS,.inorder=T) %do% {
	tag=0
	itag=0
	while(!tag){
	folds_d=sample(1:nfold,length(ind3_d),replace=T); folds_n=sample(1:nfold,length(ind3_n),replace=T);
	tnn=data.frame(table(folds_n)); td=data.frame(table(folds_d));
	u=(length(unique(folds_d))==nfold)*1; v=(length(unique(folds_n))==nfold)*1;
	tag=u*v*((sum((tnn$Freq>4)*1)==nfold)*1)*((sum((td$Freq>4)*1)==nfold)*1); 
	itag=itag+1; print(c(itag,tag)); }

foldsamples=c(folds_d,folds_n);

for(ifold in 1:nfold) {

	val1_d=which(folds_d==ifold); 	val1_n=which(folds_n==ifold);
	trn1_d=which(folds_d!=ifold); 	trn1_n=which(folds_n!=ifold);

	train=c(ind3_d[trn1_d],ind3_n[trn1_n])
	test=c(ind3_d[val1_d],ind3_n[val1_n]) # val has indices of istat

	istat_trn=istat.x[train]; istat_val=istat.x[test]; 

	expr.TEST.im=expr.pre[,test]
	expr.TRAIN.im=expr.pre[,train]

	expr.TEST.original=expr.org[,test]
	expr.TRAIN.original=expr.org[,train]	

	expr.TRAIN.blood=blood.expr[,train]
	expr.TEST.blood=blood.expr[,test]
	
		x.pred=append(x.pred,cvlasso(t(expr.TRAIN.im),t(expr.TEST.im),istat_trn,istat_val))
		x.raw=append(x.raw,cvlasso(t(expr.TRAIN.original),t(expr.TEST.original),istat_trn,istat_val))
		x.blood=append(x.blood,cvlasso(t(expr.TRAIN.blood),t(expr.TEST.blood),istat_trn,istat_val))
		

} #nfold loop ends
		
print("step2")
asd=rbind(x.pred,x.raw,x.blood);
return(list(asd,foldsamples))
} #---ENS loop

eachens=NULL;foldsamples.x=NULL;
for(i in seq(ENS)){eachens=cbind(eachens,eachens.x[[i]][[1]]);
foldsamples.x=cbind(foldsamples.x,eachens.x[[i]][[2]]);}


result[["predicted-AUC"]]=eachens[1,];
result[["original-AUC"]]=eachens[2,];
result[["Blood-AUC"]]=eachens[3,];
result[["ensfold"]]=foldsamples.x;
return(result)
}



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



