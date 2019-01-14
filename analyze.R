LF=nw.prior$Latent.Var
LF2=nw.prior$Confidence
 xbar=mean(c(LF))
 netMat=matrix(0, ncol=length(testGenes), nrow=length(testGenes))
 for(i in 1:length(testGenes)){
 	for(j in 1:length(testGenes)){
 		z=LF[i,j]/LF2[i,j]
 		netMat[i,j]=1-pnorm(z)
 	}
 }
diag(netMat)=1

net=threshold(netMat, 0.99, high=F)
TP = sum(net[cand.net.mat==1]==1,na.rm=TRUE)
TN = sum(net[cand.net.mat==0]==0,na.rm=TRUE)
FP = sum(net[cand.net.mat==0]==1,na.rm=TRUE)
FN = sum(net[cand.net.mat==1]==0,na.rm=TRUE)
TN/(TN+FP)

TP/(TP+FN)
as.numeric(auc(roc(c(net), c(cand.net.mat))))
net=threshold(netMat, 0.60, high=F)
plot(roc(c(net), c(cand.net.mat)))

##############################################################
 ROCauc = list()
 TPR = list()
 SP = list()
 C=list()
 th=c(0.99, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50,0.40)
 for(i in 1:length(th)){
 net=threshold(netMat, th[i], high=F)
 C[[i]]=cor(c(net),c(cand.net.mat))
 ROCauc[[i]]=as.numeric(auc(roc(c(net), c(cand.net.mat))))
 #ROCauc
 #net=threshold(netMat, 0.95, high=F)
 TP = sum(net[cand.net.mat==1]==1,na.rm=TRUE)
 TN = sum(net[cand.net.mat==0]==0,na.rm=TRUE)
 FP = sum(net[cand.net.mat==0]==1,na.rm=TRUE)
 FN = sum(net[cand.net.mat==1]==0,na.rm=TRUE)
 TPR[[i]]=TP/(TP+FN)
 SP[[i]]=TN/(TN+FP)
}
TPR=unlist(TPR)
SP=unlist(SP)
AUC=unlist(ROCauc)
unlist(C)
PV.cutoff=1-th
###############################################################

diag(NOP)=0
ROCauc1 = list()
TPR1 = list()
SP1= list()
PCut=c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99, 1)
for(i in 1:length(PCut)){
net=threshold(NOP, PCut[i], high=T)
ROCauc1[[i]]=as.numeric(auc(roc(c(net), c(cand.net.mat))))
#ROCauc
#net=threshold(netMat, 0.95, high=F)
TP1 = sum(net[cand.net.mat==1]==1,na.rm=TRUE)
TN1= sum(net[cand.net.mat==0]==0,na.rm=TRUE)
FP1 = sum(net[cand.net.mat==0]==1,na.rm=TRUE)
FN1 = sum(net[cand.net.mat==1]==0,na.rm=TRUE)
TPR1[[i]]=TP1/(TP1+FN1)
SP1[[i]]=TN1/(TN1+FP1)
}
unlist(TPR1)
unlist(SP1)
unlist(ROCauc1)
#prob.cutoff=PCut
################################
Indep.prior = function(data, testGenes){
	data2=data
	IPr = matrix(1, nrow = length(testGenes), ncol = length(testGenes))
	for(i in 1:ncol(data)){
		mat = matrix(data[,i], nrow = length(testGenes), ncol = length(testGenes))
		IPr=IPr*(mat)
	}
	return(IPr/10^(ncol(data)))
}

##
avg=function(data, testGenes){
	data2=data
	Av = matrix(0, nrow = length(testGenes), ncol = length(testGenes))
	for(i in 1:ncol(data)){
		mat = matrix(data[,i], nrow = length(testGenes), ncol = length(testGenes))
		Av=IPr+(mat)
	}
	return(Av/(ncol(data)))
}

#####
diag(IPr)=0
ROCauc2 = list()
TPR2 = list()
SP2= list()
PCut=c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99, 1)
for(i in 1:length(PCut)){
net=threshold(IPr, PCut[i], high=T)
ROCauc2[[i]]=as.numeric(auc(roc(c(net), c(cand.net.mat))))
#ROCauc
#net=threshold(netMat, 0.95, high=F)
TP2 = sum(net[cand.net.mat==1]==1,na.rm=TRUE)
TN2= sum(net[cand.net.mat==0]==0,na.rm=TRUE)
FP2 = sum(net[cand.net.mat==0]==1,na.rm=TRUE)
FN2 = sum(net[cand.net.mat==1]==0,na.rm=TRUE)
TPR2[[i]]=TP2/(TP2+FN2)
SP2[[i]]=TN1/(TN2+FP2)
print(c(PCut[i], TPR2[[i]], SP2[[i]], ROCauc2[[i]]))
}
TPR2=unlist(TPR2)
SP2=unlist(SP2)
AUC2=unlist(ROCauc2)

###
unlist(TPR1)
unlist(SP1)
unlist(ROCauc1)
################################
unlist(TPR)
unlist(SP)
unlist(ROCauc)
1-th

################################################################

> unlist(TPR1)
 [1] 0.3362832 0.3362832 0.3362832 0.3362832 0.3362832 0.3362832 0.3362832
 [8] 0.3362832 0.3362832 0.3362832 0.3362832 0.3362832
> unlist(SP1)
 [1] 0.359337 0.359337 0.359337 0.359337 0.359337 0.359337 0.359337 0.359337
 [9] 0.359337 0.359337 0.359337 0.359337
> unlist(ROCauc)
 [1] 0.5041232 0.5007745 0.5011194 0.5047883 0.5073414 0.5107908 0.5107120
 [8] 0.5092123 0.5123348 0.5158596 0.5158485 0.5158374

