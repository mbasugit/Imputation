
library(data.table)
args<-commandArgs(TRUE);

nfold=as.numeric(args[1]); #fold
ENS=as.numeric(args[2]); #ENS

#nfold=5; ENS=50; #one can change this parameters, nfold is partision for CV and ENS is number of independent runs

dir="/home/Desktop" #example of dir that needs to be given 
workdir=paste0(dir,"/imputation"); 
setwd(workdir)
#files for predicted expression 
folder.GSp=paste0(workdir,"/output/given"); #if user do not have predicted expression of target tissue use the one provided here

source(paste0(workdir,'/code/function.r'));



#load gtex data
load(paste0(workdir,"/input/gtexdata_v6.RData"));
gene=gtex.pc$Name

#load disease tissue and sample count
load(paste0(workdir,"/input/disease_blood_32tissue_common_sample.Rdata"))

disease=names(comm_sample); dis.tissue=list();samplecut=25;

for(dis in disease){
	tmp=comm_sample[[dis]]; print(dis);
	k=which(as.numeric(tmp[,2])>=samplecut & as.numeric(tmp[,3])>=samplecut)
	dis.tissue[[dis]]=tmp[k,1];
}
disease.study=names(dis.tissue)[lapply(dis.tissue,length)>0]


tss1="Whole Blood"
files=list.files(path=folder.GSp,pattern=".Rdata");

preddisllr1=list();

#for(dis in disease.study){
#tissue=dis.tissue[[dis]];
#for(tss in tissue){


#For all disease and tissue disease prediction uncomment above 3 lines, and tissue disease loop end line
#and comment the following line
dis="MHHTN"; tissue=c("Artery - Tibial","Adipose - Subcutaneous");

fl=files[grep(tss,files,fixed=TRUE)]
tss2=tss; print(tss2);

#load predicted expression
load(paste0(folder.GSp,"/",fl))
mat.pred=gene.regress.glmnet[["prediction"]]; 

if(ncol(mat.pred)<length(gene)){   
nn=as.numeric(colnames(mat.pred)); colnames(mat.pred)=gene[nn]; patients=rownames(mat.pred); geneid=nn;
} else { colnames(mat.pred)=gene; patients=rownames(mat.pred); geneid=seq(length(gene))}

j1=c(); j2=c(); pat.com=c();
	for (ipat in patients){
		i1=which(expc.nt$patient==ipat & expc.nt$SMTSD==tss1);
		i2=which(expc.nt$patient==ipat & expc.nt$SMTSD==tss2);
		if (length(i1)>0 & length(i2)>0){ j1=append(j1,i1[1]); j2=append(j2,i2[1]); pat.com=append(pat.com,ipat);}
				 }

nn1=length(j1); nn2=length(j2); it1=matrix(0, 1, nn1); it2=matrix(0, 1, nn2);
it1[1,]=j1; it2[1,]=j2;
gtex.tss1=gtex.pc[,c(it1[1,]), with=F]; #expr of blood
gtex.tss2=gtex.pc[,c(it2[1,]), with=F]; #expr of target tissue

pat=pat.com; mhh2=pheno.dt[,dis,with=F]; subj=pheno.subj
istat=matrix(-1, 1, length(pat))

for (i in seq(length(pat))){  ii=which(pat[i]==subj);  istat[i]=mhh2[ii];  }
istat[istat==0]=-1; istat[istat==99]=0; istat=as.numeric(istat);

k=which(istat==1|istat==-1)
temp.p=(t(mat.pred))
temp.o=as.matrix(gtex.tss2); 
temp.p=temp.p[,k]; temp.o=temp.o[,k]; istat=istat[k]; blood.oo=as.matrix(gtex.tss1)[,k]; 

#LLR gene
llrfdr=0.05
fdr=p.adjust(gene.select.LLR[['pval']],method="BH")
x2=gene.regress.glmnet[["Mean-end"]]; predPCC=rowMeans(x2[["cor"]]);
llrgene=which(fdr<=llrfdr & predPCC>0.3);

#print(c(dis,tss,length(llrgene)))

#perform prediciton if tere is minimum of 20 llr genes
if(length(llrgene)>20){ 
	temp.p1=temp.p[llrgene,]; temp.o1=temp.o[llrgene,]; blood.oo1=blood.oo[llrgene,];

	ptm <- proc.time()
	z=prediction_dis_CVfeature(temp.p1,temp.o1,blood.oo1,istat,nfold,ENS,geneid) #using age,race,gender
	print(proc.time() - ptm)
	preddisllr1[[dis]][[tss2]][[toString(paste0(llrfdr))]]=z;
	} else { 
  	preddisllr1[[dis]][[tss2]][[toString(paste0(llrfdr))]]=NA;
  	}

save(preddisllr1,file=paste0(workdir,"/output/given/dispred_ens50.Rdata"))

#} #tissue loop ends
#} #disease loop ends




