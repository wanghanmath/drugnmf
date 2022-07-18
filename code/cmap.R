cmap<-function(query,ref){
	upset<-as.matrix(query[which(query>0),])
	upset_order<-order(upset,decreasing=F)
	upsetnew<-as.matrix(upset[upset_order,])
	upset<-upsetnew
	downset<-as.matrix(query[which(query<0),])
	downset_order<-order(downset,decreasing=F)
	downsetnew<-as.matrix(downset[downset_order,])
	downset<-downsetnew
	ref_order<-order(ref,decreasing=F)
    refnew<-as.matrix(ref[ref_order,])
    upset_aset<-upset_bset<-rep(0,length(upset))
    for (i in 1:length(upset)){
    	loc<-which(rownames(refnew)==rownames(upset)[i])
    	upset_aset[i]<-(i/dim(upset)[1])-(loc/dim(refnew)[1])
    	upset_bset[i]<-(loc/dim(refnew)[1])-((i-1)/dim(upset)[1])
    }
    downset_aset<-downset_bset<-rep(0,length(downset))
    for (i in 1:length(downset)){
    	loc<-which(rownames(refnew)==rownames(downset)[i])
    	downset_aset[i]<-(i/dim(downset)[1])-(loc/dim(refnew)[1])
    	downset_bset[i]<-(loc/dim(refnew)[1])-((i-1)/dim(downset)[1])
    }
    upset_a<-max(upset_aset)
    upset_b<-max(upset_bset)
    downset_a<-max(downset_aset)
    downset_b<-max(downset_bset)
    if (upset_a>upset_b){
    	ks_up<-upset_a
    }
    if (upset_a<upset_b){
    	ks_up<-(-upset_b)
    }
    if (downset_a>downset_b){
    	ks_down<-downset_a
    }
    if (downset_a<downset_b){
    	ks_down<-(-downset_b)
    }
    if (sign(ks_up)==sign(ks_down)){
    	Con_score<-abs(ks_down)-abs(ks_up)
    }
    if (sign(ks_up)!=sign(ks_down)){
    	Con_score<-ks_down-ks_up
    }
    return(Con_score)
}    
