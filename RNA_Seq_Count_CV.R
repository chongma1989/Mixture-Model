##===============================================##
## Purpose: RNA-Seq count study for the microarray
##          project using cross-validation. 
## Author : Chong Ma <chongm@email.sc.edu>
## Updated_Date_On: January 1st, 2018
##===============================================##

## read the arguments
args=commandArgs(TRUE)

## load relevant libraries
.libPaths("~/Rlibs")
# source("https://www.bioconductor.org/biocLite.R")
# biocLite("tweeDEseqCountData")
library("tweeDEseqCountData")
library("edgeR")
library("statmod")
library("sda")
library("maxLik")
library("foreach")
library("doParallel")
library("doSNOW")
source("Headfile_lfdr1.R")

#==============================================#
# obtain and filter the RNA-seq count data     #
#==============================================#
data(pickrell1)
Gender=pickrell1.eset$gender
y=readRDS("pickrell_trim.rds")

#==============================================#
# approximately splits samples as training and 
# verfication data sets 
#==============================================#
Tsample=sort(c(sample(which(Gender=="female"),20),
               sample(which(Gender=="male"),15)))
Vsample=(1:69)[-Tsample]

# Training samples
Ty=y[,Tsample]
# estimating the dispersion
Tdesign=model.matrix(~Gender[Tsample])
Ty=estimateDisp(Ty,Tdesign,robust=TRUE)
Tfit=glmQLFit(Ty,Tdesign,robust = TRUE)
Tqlf=glmQLFTest(Tfit)
quants=Tqlf$table$PValue

# Verfication samples
Vy=y[,Vsample]
# estimating the dispersion
Vdesign=model.matrix(~Gender[Vsample])
Vy=estimateDisp(Vy,Vdesign,robust=TRUE)
Vfit=glmQLFit(Vy,Vdesign,robust = TRUE)
Vqlf=glmQLFTest(Vfit)
quantv=Vqlf$table$PValue

# op<-par(mfrow=c(1,2))
# hist(quants,breaks=100,freq=FALSE,main="")
# hist(quantv,breaks=100,freq=FALSE,main="")
# par(op)
#==============================================#
# Use proposed methods by cross-validation     #
#==============================================#
##used empiri.null in the headfile
emnull.Train=empiri.null(quants,c(1,2),c(.99,.01))
# get the raw weights and parameters for the 
# uniform-beta mixture model
parm=emnull.Train$parm 
tau=emnull.Train$tau

#==============================================#
A11=integrate(f11,0,1,stop.on.error = FALSE)$value  ##reset the tolerance at 1e-8
new.tau=c(1-tau[2]*A11,tau[2]*A11)   ##makes results more precise

##get the likelihood for the verification data using 
##the adjusted empirical null model?? Need clarify this...
f0.verf=new.f0(quantv)
f1.verf=new.f1(quantv)
fw.verf=new.fw(quantv)
Tverf=mean(h(quantv))

##calculate the tilted distribution
theta=tryCatch(uniroot(TiltGmean,interval=c(-10,10))$root,
               error=function(e) 0)
tilted.tau=tryCatch(tilt.tau(theta),error=function(e) tau)

# cl<-makeCluster(4)
# registerDoParallel(cl)
str=Sys.time()
##ROC curve for verified ptvalues
verf.ROC1<-ROC1(quantv,rownames(Vqlf))

##ROC curve for verified ptvalues
verf.ROC2<-ROC2(quantv,rownames(Vqlf))
Sys.time()-str
# stopCluster(cl)

#====================#
#    save results
cv.parm=parm
cv.tau1=tau
cv.tau2=new.tau
cv.tau3=tilted.tau
RNA_Seq.result=list(verf.ROC1,verf.ROC2,cv.parm,cv.tau1,cv.tau2,cv.tau3,Vsample)
saveRDS(RNA_Seq.result,file=paste0("./RNA_Seq_Cluster_lfdr/RNA_Seq_CV",args,".rds")) 


