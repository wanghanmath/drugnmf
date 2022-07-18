#integrate by SNF#
library("SNFtool")
cor_faers<-read.csv("/home/r9user10/Documents/clinical_trial/COR_FAERS_casimerge.csv",header=T,stringsAsFactors=F)
cor_assay<-read.csv("/home/r9user10/Documents/clinical_trial/casi_bam/MDS/COR_ASSAY_update.csv",header=T,stringsAsFactors=F)
cor_struc<-read.csv("/home/r9user10/Documents/clinical_trial/casi_bam/MDS/COR_STRUCTURE_update.csv",header=T,stringsAsFactors=F)
cor_expre<-read.csv("/home/r9user10/Documents/clinical_trial/casi_bam/MDS/COR_LINCS_update.csv",header=T,stringsAsFactors=F)
cor_faers<-cor_faers[,2:dim(cor_faers)[2]]
cor_assay<-cor_assay[,2:dim(cor_assay)[2]]
cor_struc<-cor_struc[,2:dim(cor_struc)[2]]
cor_expre<-cor_expre[,2:dim(cor_expre)[2]]
rownames(cor_faers)<-colnames(cor_faers)
rownames(cor_assay)<-colnames(cor_assay)
rownames(cor_struc)<-colnames(cor_struc)
rownames(cor_expre)<-colnames(cor_expre)
for (i in 2:dim(cor_expre)[1]){
	for (j in 1:(i-1)){
		cor_expre[i,j]<-cor_expre[j,i]
	}
	cor_expre[i,i]<-1
}
cor_expre[1,1]<-1
cor_expre<-as.matrix(cor_expre)
cor_expre[which(cor_expre<0)]<-0
cor_expre<-cor_expre/max(cor_expre)
faers_name<-colnames(cor_faers)
assay_name<-toupper(colnames(cor_assay))
struc_name<-toupper(colnames(cor_struc))
expre_name<-toupper(colnames(cor_expre))
int1<-intersect(faers_name,assay_name)
int2<-intersect(int1,struc_name)
int3<-intersect(int2,expre_name)
faers_loc<-assay_loc<-struc_loc<-expre_loc<-rep(0,length(int3))
for (i in 1:length(int3)){
	faers_loc[i]<-which(faers_name==int3[i])
	assay_loc[i]<-which(assay_name==int3[i])
	struc_loc[i]<-which(struc_name==int3[i])	
    expre_loc[i]<-which(expre_name==int3[i])
}
cor_faersnew<-cor_faers[faers_loc,faers_loc]
cor_assaynew<-cor_assay[assay_loc,assay_loc]
cor_strucnew<-cor_struc[struc_loc,struc_loc]
cor_exprenew<-cor_expre[expre_loc,expre_loc]
rownames(cor_assaynew)<-colnames(cor_assaynew)<-c(int3)
rownames(cor_strucnew)<-colnames(cor_strucnew)<-c(int3)
rownames(cor_exprenew)<-colnames(cor_exprenew)<-c(int3)
cor_inte<-matrix(ncol=dim(cor_faersnew)[1],nrow=dim(cor_faersnew)[1])
for (i in 1:dim(cor_inte)[1]){
	for (j in 1:dim(cor_inte)[1]){
		cor_inte[i,j]<-cor_faersnew[i,j]*cor_assaynew[i,j]*cor_strucnew[i,j]*cor_exprenew[i,j]
	}
}
rownames(cor_inte)<-c(int3)
colnames(cor_inte)<-c(int3)
cor_faers1<-as.matrix(1-cor_faersnew)
cor_assay1<-as.matrix(1-cor_assaynew)
cor_struc1<-as.matrix(1-cor_strucnew)
cor_expre1<-as.matrix(1-cor_exprenew)
K = 10 
alpha = 0.5 
T = 20
W_faers<-affinityMatrix(cor_faers1, K, alpha)
W_assay<-affinityMatrix(cor_assay1, K, alpha)
W_struc<-affinityMatrix(cor_struc1, K, alpha)
W_expre<-affinityMatrix(cor_expre1, K, alpha)
W = SNF(list(W_faers,W_assay,W_struc,W_expre), K, T)
SNFname<-colnames(cor_faersnew)[order(W[94,],decreasing=T)[1:20]]
rownames(W)<-rownames(cor_inte)
colnames(W)<-rownames(cor_inte)
