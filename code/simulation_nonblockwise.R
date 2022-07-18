#simulation with non-block-wise pattern#
library("DirichletReg")
library("NMF")
library("combinat")
pattern<-matrix(ncol=50,nrow=5)
prob<-matrix(ncol=10,nrow=5)
prob[1,]<-c(0.02736860,0.12541887,0.05641219,0.07464338,0.09908623,0.06283853,0.16954483,0.08033740,0.12642745,0.17792252)
prob[2,]<-c(0.05518038,0.03670391,0.20534650,0.01779860,0.19465400,0.15284080,0.12575085,0.02593506,0.11578989,0.07000000)
prob[3,]<-c(1.713488e-01,5.411259e-02,1.476365e-03,2.709419e-03,1.494267e-01,6.800434e-02,1.616755e-13,3.024545e-01,1.602027e-13,2.504672e-01)
prob[4,]<-c(0.03419887,0.02268002,0.05456208,0.03440245,0.19438291,0.21710751,0.02188385,0.13304830,0.16673734,0.12099668)
prob[5,]<-c(0.02736860,0.12541887,0.05641219,0.07464338,0.09908623,0.06283853,0.16954483,0.08033740,0.12642745,0.17792252)
pattern[1,]<-c(c(prob[1,]/2),c(prob[2,]/2),rep(0,30))
pattern[2,]<-c(c(prob[2,]/2),c(prob[1,]/2),rep(0,30))
pattern[3,]<-c(rep(0,20),c(prob[3,]/2),c(prob[4,]/2),rep(0,10))
pattern[4,]<-c(rep(0,20),c(prob[4,]/2),c(prob[3,]/2),rep(0,10))
pattern[5,]<-c(rep(0,30),c(prob[3,]/2),c(prob[5,]/2))
pattern_array<-array(dim=c(5,50,101))
pattern_array[,,1]<-pattern
drug_array<-array(dim=c(100,10,100))
drug_prop<-matrix(rep(0,500),ncol=5,nrow=100)
drug_num<-rep(0,100)
w_simplex<-h_simplex<-rep(0,100)
w_rank<-h_rank<-rep(0,100)
pattern_rank<-matrix(ncol=20,nrow=5)
pattern_rank[1,]<-c(1:20)
pattern_rank[2,]<-c(1:20)
pattern_rank[3,]<-c(21:40)
pattern_rank[4,]<-c(21:40)
pattern_rank[5,]<-c(31:50)
for (time in 1:100){
	poisson_num<-10
	drug_prop[1:20,1:2]<-rdirichlet(20,alpha=c(0.5,0.5))
    drug_prop[21:40,3:4]<-rdirichlet(20,alpha=c(0.5,0.5))
    drug_prop[41:60,c(1,3,5)]<-rdirichlet(20,alpha=c(0.2,0.2,0.6))
    drug_prop[61:80,c(2,4,5)]<-rdirichlet(20,alpha=c(0.3,0.4,0.3))
    drug_prop[81:100,1:5]<-rdirichlet(20,alpha=c(0.1,0.2,0.3,0.3,0.1))
    drug_array[1:100,1:5,time]<-drug_prop
	drug_num[1:20]<-floor(runif(20,1500,2000))
	drug_num[21:40]<-floor(runif(20,1000,1500))
	drug_num[41:60]<-floor(runif(20,1500,2000))
	drug_num[61:80]<-floor(runif(20,2000,2500))
	drug_num[81:100]<-floor(runif(20,2500,3000))
	drug_prop_num<-matrix(ncol=50,nrow=100)
	for (i in 1:100){
		proportion<-rep(0,50)
		for (j in 1:5){
			proportion<-proportion+c(drug_prop[i,j]*pattern[j,])
		}
		drug_prop_num[i,]<-c(drug_num[i]*proportion)
		poisson_pos<-sample(1:50)[1:10]
		drug_prop_num[i,poisson_pos]<-drug_prop_num[i,poisson_pos]+rpois(10,5)
	}
	norm<-rep(0,50)
	w_array<-array(dim=c(50,5,50))
	h_array<-array(dim=c(5,100,50))
	for (iter in 1:50){
		res<-nmf(t(drug_prop_num),5)
		w<-basis(res)
		h<-coef(res)
		aaa<-w%*%h
		norm1<-norm(t(drug_prop_num)-aaa,"F")
		norm[iter]<-norm1
		w_array[,,iter]<-w
		h_array[,,iter]<-h
	}
	iter_max<-which.max(norm)
	wnew<-w_array[,,iter_max]
	hnew<-h_array[,,iter_max]
    w1<-wnew
    h1<-hnew
    for (i in 1:dim(w1)[2]){
    	w1[,i]<-w1[,i]/sum(w1[,i])
    }
    for (i in 1:dim(h1)[2]){
    	h1[,i]<-h1[,i]/sum(h1[,i])
    }
    permut<-permn(c(1:5))
    w_norm<-h_norm<-rep(0,length(permut))
    for (permut_time in 1:length(permut)){
    	w1_permut<-w1[,c(permut[[permut_time]])]
    	h1_permut<-h1[c(permut[[permut_time]]),]
    	w_norm[permut_time]<-norm(t(pattern)-w1_permut,"F")
    	h_norm[permut_time]<-norm(t(drug_prop)-h1_permut,"F")
    }
    permut_max<-which.min(w_norm)
    w1_final<-w1[,c(permut[[permut_max]])]
    h1_final<-h1[c(permut[[permut_max]]),]
    pattern_array[,,(time+1)]<-t(w1_final)
    drug_array[,6:10,time]<-t(h1_final)
    w_rank_sub<-rep(0,5)
    h_rank_sub<-rep(0,100)
    for (i in 1:5){
    	w1_final_rank<-order(w1_final[,i],decreasing=T)[1:20]
    	w_rank_sub[i]<-length(intersect(w1_final_rank,c(pattern_rank[i,])))
    }
    for (i in 1:100){
    	if (i<=40){
    		h1_final_rank<-order(h1_final[,i],decreasing=T)[1:2]
    		h1_origi_rank<-order(drug_prop[i,],decreasing=T)[1:2]
    		h_rank_sub[i]<-length(intersect(h1_final_rank,h1_origi_rank))/2
    	}
    	if (i>40&&i<=80){
    		h1_final_rank<-order(h1_final[,i],decreasing=T)[1:3]
    		h1_origi_rank<-order(drug_prop[i,],decreasing=T)[1:3]
    		h_rank_sub[i]<-length(intersect(h1_final_rank,h1_origi_rank))/3
    	}
    	if (i>80&&i<=100){
    		h1_final_rank<-order(h1_final[,i],decreasing=T)[1:5]
    		h1_origi_rank<-order(drug_prop[i,],decreasing=T)[1:5]
    		h_rank_sub[i]<-length(intersect(h1_final_rank,h1_origi_rank))/5
    	}
    }
    w_rank[time]<-mean(w_rank_sub)
    h_rank[time]<-mean(h_rank_sub)
    w_simplex[time]<-norm(t(pattern)-w1_final,"F")/5
    h_simplex[time]<-norm(t(drug_prop)-h1_final,"F")/100
}
