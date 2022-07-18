#imputation by MDS#
cor_faers<-read.csv("COR_FAERS_casimerge.csv",header=T,stringsAsFactors=F)
cor_faers<-cor_faers[,2:dim(cor_faers)[2]]
rownames(cor_faers)<-colnames(cor_faers)
datasets<-c("COR_ASSAYnew.csv","COR_STRUCTURE.csv","COR_LINCS_overall.csv")
datasetsname<-c("COR_ASSAY","COR_STRUCTURE","COR_LINCS")
assayK_data<-c(39,44,49)
for (data in 1:3){
	cor_assay<-read.csv("./first-stage/",datasets[data],header=T,stringsAsFactors=F)
	cor_assay<-cor_assay[,2:dim(cor_assay)[2]]
    rownames(cor_assay)<-colnames(cor_assay)
    colnames(cor_assay)<-toupper(colnames(cor_assay))
    rownames(cor_assay)<-toupper(rownames(cor_assay))
	if (data==3){
		for (i in 2:dim(cor_assay)[2]){
			for (j in 1:(i-1)){
				cor_assay[i,j]<-cor_assay[j,i]
			}
			cor_assay[i,i]<-1
		}
        cor_assay[1,1]<-1
	}
	cor_faers<-as.matrix(cor_faers)
	cor_faers[which(cor_faers<0)]<-0
	cor_faers_dist<-1-cor_faers
	cor_assay<-as.matrix(cor_assay)
	cor_assay[which(cor_assay<0)]<-0
	cor_assay<-cor_assay/max(cor_assay)
	cor_assay_dist<-1-cor_assay
	cor_assay2<-cor_assay_dist*cor_assay_dist
	cor_faers2<-cor_faers_dist*cor_faers_dist
	n<-dim(cor_assay)[2]
	n_faers<-dim(cor_faers)[2]
	H1<-1/n*as.matrix(rep(1,n))%*%t(as.matrix(rep(1,n)))
	H<-diag(n)-H1	
	cor_tomds_assay=-1/2*H%*%cor_assay_dist%*%H
	H1_faers<-1/n_faers*as.matrix(rep(1,n_faers))%*%t(as.matrix(rep(1,n_faers)))
	H_faers<-diag(n_faers)-H1_faers
	cor_tomds_faers=-1/2*H_faers%*%cor_faers_dist%*%H_faers
	eigen_faers<-eigen(cor_tomds_faers)
	eigen_assay<-eigen(cor_tomds_assay,symmetric=T)
	faersK<-5
	assayK<-assayK_data[data]
faers_reduc<-matrix(ncol=faersK,nrow=dim(cor_faers)[1])
assay_reduc<-matrix(ncol=assayK,nrow=dim(cor_assay)[1])
for (i in 1:dim(faers_reduc)[2]){
	faers_reduc[,i]<-sqrt(eigen_faers$values[i])*eigen_faers$vectors[,i]
}
for (i in 1:dim(assay_reduc)[2]){
	assay_reduc[,i]<-sqrt(eigen_assay$values[i])*eigen_assay$vectors[,i]
}
rownames(faers_reduc)<-rownames(cor_faers)
rownames(assay_reduc)<-toupper(rownames(cor_assay))
drug_int<-intersect(rownames(cor_faers),rownames(assay_reduc))
loc_faers<-loc_assay<-rep(0,length(drug_int))
for (i in 1:length(loc_faers)){
	loc_faers[i]<-which(rownames(faers_reduc)==drug_int[i])
	loc_assay[i]<-which(rownames(assay_reduc)==drug_int[i])
}
faers_reducnew<-faers_reduc[loc_faers,]
assay_reducnew<-assay_reduc[loc_assay,]
faers_assay<-cbind(assay_reducnew,faers_reducnew)
colnames(faers_assay)<-c(paste0("assay",1:assayK),paste0("faers",1:faersK))
coefficient<-matrix(ncol=faersK+1,nrow=assayK)
for (i in 1:assayK){
	relation<-lm(faers_assay[,i]~faers_assay[,(assayK+1)]+faers_assay[,(assayK+2)]+faers_assay[,(assayK+3)]+faers_assay[,(assayK+4)]+faers_assay[,(assayK+5)])
	coefficient[i,]<-c(relation$coefficients)
}
predict_casi<-matrix(ncol=assayK,nrow=2)
drugname<-c("CASIRIVIMAB.AND.IMDEVIMAB","BAMLANIVIMAB")
for (kk in 1:length(drugname)){
	drug<-drugname[kk]
	drug_loc<-which(rownames(faers_reduc)==drug)
	for (i in 1:dim(predict_casi)[2]){
		predict_casi[kk,i]<-coefficient[i,1]+coefficient[i,2]*faers_reduc[drug_loc,1]+coefficient[i,3]*faers_reduc[drug_loc,2]+coefficient[i,4]*faers_reduc[drug_loc,3]+coefficient[i,5]*faers_reduc[drug_loc,4]+coefficient[i,6]*faers_reduc[drug_loc,5]
	}
    cor_casi<-rep(0,dim(assay_reduc)[1])
    for (i in 1:dim(assay_reduc)[1]){
    	bb<-c(predict_casi[kk,])-c(assay_reduc[i,])
    	bb_norm<-norm(as.matrix(bb),"F")^2
    	cor_casi[i]<-1-bb_norm
    }	
    cor_casi<-as.matrix(cor_casi)
    rownames(cor_casi)<-rownames(assay_reduc)
    if (kk==1){
    	cor_casinew<-cor_casi       
    }
    if (kk==2){
    	cor_bamnew<-cor_casi
    }
    write.csv(cor_casi,paste0("cor_",drug,"_",datasetsname[data],".csv"))
}
casi_bam_matrix<-cbind(cor_casinew,cor_bamnew)
colnames(casi_bam_matrix)<-c(drugname)
casi_bam<-predict_casi[1,]%*%predict_casi[2,]
casi_bam_mat<-matrix(ncol=2,nrow=2)
casi_bam_mat[1,1]<-casi_bam_mat[2,2]<-1
casi_bam_mat[1,2]<-casi_bam_mat[2,1]<-casi_bam
rownames(casi_bam_mat)<-colnames(casi_bam_mat)<-c(drugname)
casi_bam_matrixnew<-rbind(casi_bam_matrix,casi_bam_mat)
cor_assay1<-cbind(cor_assay,casi_bam_matrix)
cor_assaynew<-rbind(cor_assay1,t(casi_bam_matrixnew))
write.csv(cor_assaynew,paste0("./MDS/",datasetsname[data],"_update.csv"))
}
