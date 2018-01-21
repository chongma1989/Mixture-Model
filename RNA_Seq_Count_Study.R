##===============================================##
## Purpose: RNA-Seq count study for the microarray
##          project using cross-validation. 
## Author : Chong Ma <chongm@email.sc.edu>
## Updated_Date_On: January 1st, 2018
##===============================================##

## load relevant libraries
.libPaths("~/Rlibs")
source("https://www.bioconductor.org/biocLite.R")
biocLite("tweeDEseqCountData")
library("edgeR")
library("statmod")
library("sda")
library("maxLik")
library("foreach")
library("doParallel")
library("doSNOW")
library("plyr")
library("igraph")
library("MASS")
library("plot3D")
library("xtable")
library("nleqslv")
source("Headfile_lfdr1.R")

#==============================================#
# obtain and filter the RNA-seq count data     #
#==============================================#
data(pickrell1);data("annotEnsembl63")
Counts=exprs(pickrell1.eset)
Gender=pickrell1.eset$gender
annot=annotEnsembl63[,c("Symbol","Chr")]
y=DGEList(counts=Counts, genes=annot[rownames(Counts),])

# keep genes with at least 1 count-per-million 
# (cpm) in at least 20 samples
isexpr=rowSums(cpm(y)>1)>=20
# keep only genes with defined annotation
hasannot= rowSums(is.na(y$genes))==0
y=y[isexpr & hasannot, , keep.lib.sizes=FALSE]
y=calcNormFactors(y) ## calculate normal factors
saveRDS(y,"pickrell_trim.rds")

#==============================================# 
data(pickrell1)
Gender=pickrell1.eset$gender
y=readRDS("pickrell_trim.rds")
barplot(y$samples$lib.size*1e-6,names=1:69,
        ylab="Library size (million)")

# estimating the dispersion
design=model.matrix(~Gender)
y=estimateDisp(y,design,robust=TRUE)
fit=glmQLFit(y,design,robust = TRUE)
qlf=glmQLFTest(fit)

plotBCV(y) 
plotQLDisp(fit)

#============================================#
# Use full data set to fit the mixture model 
#============================================#
quants=qlf$table$PValue

#=============================#
# Method of moment estimation #
#=============================#
mu=c(mean(quants),mean(quants^2),mean(quants^3))
MOM(quants,start=c(0.99,1,3))
#=============================#
#      EM algorithm           #
#=============================#
emnull.Train=empiri.null(quants,c(1,1),c(.7,.3))
##get the raw weights and parameters for the uniform-beta mixture model
parm=emnull.Train$parm
tau=emnull.Train$tau

# get FDR using the whole genes 
A11=integrate(f11,0,1,stop.on.error = FALSE)$value  ##reset the tolerance at 1e-8
new.tau=c(1-tau[2]*A11,tau[2]*A11)   ##makes results more precise
cl<-makeCluster(4)
registerDoParallel(cl)
str=Sys.time()
whole.ROC1<-data.frame(ROC1(quants,rownames(qlf$table)))
stopCluster(cl)

qlf_table=data.frame(qlf$table,
                     FDR1=p.adjust(qlf$table$PValue,method="BH"),
                     FDR2=whole.ROC1$FDR)
# qlf_table=qlf_table[order(qlf_table$FDR1),]

#============================================#
#  graphical displays comparisons
#============================================#
pdf("whole_pvalue1.pdf")
op<-par(mar=c(4,4.1,4.1,1.1),tck=-.01,cex=1)
hist(qlf$table$PValue,breaks = seq(0,1,0.01),col = "grey70",
     freq=FALSE,xlab="",ylab="",main="",xaxs="i",yaxs="i",
     xaxt="n",yaxt="n",ylim=c(0,2))
axis(1,at=seq(0,1,0.1),line=0,padj=-1.25)
axis(2,at=seq(0,2,by=0.5),line=0,hadj=0.5,las=2)
mtext(side=1,text="x (p-values from lrt-statistic)",
      line=1.5,font = 1)
mtext(side=2,text="Density",line=1.75,font = 1)
box()
curve(fw(x),0,1,col="red",lwd=3,add=T)
curve(dunif(x),0,1,col="blue",lwd=3,lty=4,add=T)
curve(dbeta(x,parm[1],parm[2]),0,1,col="blue",lwd=3,lty=2,add=T)

legend(.05,2,legend=c(expression(paste(f[scriptscriptstyle(0)]^{"*"})),
                      expression(paste(f[scriptscriptstyle(1)]^{"*"})),
                      expression(f)),
       col=c("blue","blue","red"),lty=c(4,2,1),lwd=2,cex = 1.25)
par(op)
dev.off()

pdf("whole_pvalue2.pdf")
op<-par(mar=c(4,4.1,4.1,1.1),tck=-.01,cex=1)
hist(qlf$table$PValue,breaks = seq(0,1,0.01),col = "grey70",
     freq=FALSE,xlab="",ylab="",main="",xaxs="i",yaxs="i",
     xaxt="n",yaxt="n",ylim=c(0,2))
axis(1,at=seq(0,1,0.1),line=0,padj=-1.25)
axis(2,at=seq(0,2,by=0.5),line=0,hadj=0.5,las=2)
mtext(side=1,text="x (p-values from lrt-statistic)",
      line=1.5,font = 1)
mtext(side=2,text="Density",line=1.75,font = 1)
box()
curve(fw(x),0,1,col="red",lwd=3,add=T)
curve(new.f0(x),0,1,col="blue",lwd=3,lty=4,add=T)
curve(new.f1(x),0,1,col="blue",lwd=3,lty=2,add=T)

legend(.05,2,legend=c(expression(paste(f[scriptscriptstyle(0)])),
                        expression(paste(f[scriptscriptstyle(1)])),
                        expression(f)),
       col=c("blue","blue","red"),lty=c(4,2,1),lwd=2,cex = 1.25)
par(op)
dev.off()

pdf("FDR_cmp1.pdf")
op<-par(mar=c(4,4.1,4.1,1.1),tck=-.01,cex=1)
curve(ecdf(qlf_table$FDR1)(x),0,1,lwd=3,lty=1,
      xlab="",ylab="",xaxs="i",yaxs="i",xaxt="n",yaxt="n")
curve(ecdf(qlf_table$FDR2)(x),0,1,lwd=3,lty=2,add=TRUE)
axis(1,at=seq(0,1,0.1),line=0,padj=-1.25)
axis(2,at=seq(0,1,by=0.2),line=0,hadj=0.5,las=2)
mtext(side=1,text="FDR",line=1.5,font = 1)
mtext(side=2,text="CDF",line=1.75,font = 1)
box()
par(op)
dev.off()

pdf("FDR_cmp2.pdf")
op<-par(mar=c(4,4.1,4.1,1.1),tck=-.01,cex=1)
curve(ecdf(qlf_table$FDR1)(x),0,0.5,lwd=3,lty=1,ylim=c(0,0.02),
      xlab="",ylab="",xaxs="i",yaxs="i",xaxt="n",yaxt="n")
curve(ecdf(qlf_table$FDR2)(x),0,0.5,lwd=3,lty=2,add=TRUE)
axis(1,at=seq(0,1,0.1),line=0,padj=-1.25)
axis(2,at=seq(0,0.02,by=0.004),line=0,hadj=0.75,las=2)
mtext(side=1,text="FDR",line=1.5,font = 1)
mtext(side=2,text="CDF",line=2,font = 1)
box()
par(op)
dev.off()

#============================#
#  cross-validation method   #
#============================#
files=list.files("./RNA_Seq_Cluster_lfdr",pattern=".rds",full.names=TRUE)
files=files[order(nchar(files))][1:100]
# cv.result=vector("list",100)
cv.result=lapply(files,function(x) readRDS(x))
FDR=lapply(cv.result[!sapply(cv.result,is.null)],
           function(x) as.data.frame(x[[1]][,1:2]))

##create significant genes with ordered FDRs
tmp.gene=lapply(FDR,function(x) x[x$FDR<=.01,])
##get which cv has genes with FDR <.01
a=lapply(tmp.gene,nrow)>0
sigcv=which(a==T)
tmp.gene1=tmp.gene[lapply(tmp.gene,nrow)>0]
tmp.gene20=lapply(tmp.gene1,function(x) na.omit(x[order(x$FDR),]))
tmp.gene2=tmp.gene20[lapply(tmp.gene20,nrow)>0]
sigcv=sigcv[unlist(lapply(tmp.gene20,nrow))>0]

## get the frequency of significant genes and the path in each 
## cross-validation significant genes which FDR <0.01 and 100 "paths"
tmp.genes2=count(unlist(lapply(tmp.gene2,rownames)))
sigenes2=tmp.genes2[order(tmp.genes2$freq,decreasing=T),]
rownames(sigenes2)=NULL
sigenes2.path=lapply(tmp.gene2,rownames)


#=============================================================#
#           construct the table of significant genes          #
#=============================================================#
##construct gene_nodes data frame
##define colors for nodes
pal1=rainbow(100,alpha=.8)
pal2=heat.colors(100,alpha=.8)
pal3=c("grey","darkgrey",rep(heat.colors(10,alpha=1)[2],9))

gene_nodes=sigenes2
gene_nodes$symbol=y$genes$Symbol[match(as.character(sigenes2$x),rownames(y$genes))]
gene_nodes$size=8+sigenes2$freq/100*10
gene_nodes$col1=pal1[factor(sigenes2$freq)]
gene_nodes$col2=pal2[factor(sigenes2$freq)]
colnames(gene_nodes)=c("ID","freq","symbol","size","col1","col2")

##output table
gene_table=cbind(y$genes[as.character(gene_nodes$ID),],
                 qlf_table[as.character(sigenes2$x),],
                 freq=gene_nodes$freq)
gene_table=gene_table[order(gene_table$PValue),]
print(xtable(gene_table,digits=c(0,0,0,2,2,2,-2,-2,-2,2)),
      include.rownames=T,sanitize.colnames.function=identity)
saveRDS(gene_table,"gene_table.rds")

gene_index=match(rownames(gene_table),rownames(y$genes))
Verfsam=lapply(cv.result[!sapply(cv.result,is.null)],function(x) x[[7]])
Verfsam1=Verfsam[sigcv];Verfsam2=Verfsam[-sigcv]
ng = nrow(sigenes2)  ## number of significant genes
n=length(Verfsam)    ## number of cv's
n1=length(Verfsam1)  ## number of significant cv's
n2=length(Verfsam2)  ## number of non-significant cv's

#==================================================#
# obtain the summary table including the gene_ID,
# gene_name, frequency, median FDR, mean, median,
# and sd of p-values from varification sets
#==================================================#
cl=makeCluster(4)
registerDoParallel(cl)
quantv=foreach(i=1:100,.combine = rbind,.multicombine = TRUE)%dopar%{
  Vsample=cv.result[[i]][[7]]
  # Verfication samples
  Vy=y[,Vsample]
  # estimating the dispersion
  Vdesign=model.matrix(~Gender[Vsample])
  Vy=estimateDisp(Vy,Vdesign,robust=TRUE)
  Vfit=glmQLFit(Vy,Vdesign,robust = TRUE)
  Vqlf=glmQLFTest(Vfit)
  
  c(Vqlf$table$PValue[gene_index])
}
rownames(quantv)=NULL
stopCluster(cl)

quantv_table=t(apply(quantv,2,function(t) c(median(t),mean(t),sd(t))))
FDR_table=t(sapply(FDR, function(t) t$FDR[gene_index]))
med_FDR=apply(FDR_table,2, function(t) median(t[sigcv]))

sum_table=cbind(gene_table[,c(1,2,9)],med.FDR=med_FDR,quantv_table)
colnames(sum_table)[5:7]=c("med(x)","avg(x)","sd(x)")
sum_table=sum_table[order(sum_table$freq,decreasing = T),]
print(xtable(sum_table,digits=c(0,0,0,2,-2,-2,-2,-2)),
      include.rownames=T,sanitize.colnames.function=identity)

##================================================##
##  create a sparse network of singificant genes  ##
##================================================##
sparse_network<-function(gene_nodes,gene_path,min_freq=5,
                         savePDF=NULL,...){
  
  ##create edges 
  create_edges<-function(x){
    if(length(x)>1){
      return(c(combn(x,2)))
    }else{
      return(rep(x,2))
    }
  }  
  
  ##remove the significant genes merely appear less than min_freq
  minor_gene=as.character(with(gene_nodes,ID[freq<=min_freq]))
  sparse_gene_path=lapply(gene_path, function(x) x[is.na(match(x,minor_gene))])
  sparse_gene_edges<-unlist(lapply(sparse_gene_path,create_edges))
  
  ## create network graph object
  sparse_gene_net<-graph(sparse_gene_edges,directed=F)
  E(sparse_gene_net)$weight=rep(1,length(E(sparse_gene_net)))
  
  ## simplify the network graph 
  sparse_gene_nets=simplify(sparse_gene_net,remove.multiple=T,
                            remove.loops=T,edge.attr.comb=c(weight="sum"))
  V(sparse_gene_nets)$name=as.character(with(gene_nodes,
                           symbol[match(V(sparse_gene_nets)$name,as.character(ID))]))
  V(sparse_gene_nets)$size=gene_nodes$size[match(V(sparse_gene_nets)$name,gene_nodes$symbol)]
  V(sparse_gene_nets)$col="grey70"
  ew=ifelse(E(sparse_gene_nets)$weight<10,0.5,
            ifelse(E(sparse_gene_nets)$weight<20,1,
                   ifelse(E(sparse_gene_nets)$weight<30,1.5,
                          ifelse(E(sparse_gene_nets)$weight<40,2,2.5))))
  
  pdf(savePDF,width = 8,height = 8)
  op<- par(mar=c(0,0,0,0))
  plot(sparse_gene_nets,
       # layout=layout.fruchterman.reingold,
       layout=layout.auto,
       vertex.color=V(sparse_gene_nets)$col,vertex.size=V(sparse_gene_nets)$size,
       vertex.frame.color="black",vertex.label.color="black",
       vertex.label.cex=.6,vertex.label.dist=0,
       edge.width=ew,edge.color="grey")
  # colkey(col=colorspace::diverge_hsv(100,alpha=0.8),clim=c(0,1),clab="",add=TRUE,
  #        length = 0.75,width=0.5,tcl=-0.25,hadj=0.5)
  par(op)
  dev.off()
}
sparse_network(gene_nodes,sigenes2.path,min_freq = 2,savePDF = "gene_network_lfdr1.pdf")
sparse_network(gene_nodes,sigenes2.path,min_freq = 50,savePDF = "gene_network_lfdr2.pdf")
sparse_network(gene_nodes,sigenes2.path,min_freq = 70,savePDF = "gene_network_lfdr3.pdf")
sparse_network(gene_nodes,sigenes2.path,min_freq = 80,savePDF = "gene_network_lfdr4.pdf")
sparse_network(gene_nodes,sigenes2.path,min_freq = 90,savePDF = "gene_network_lfdr5.pdf")

##======================================##
##  stability selection study(tilted)   ##
##======================================##
t=seq(0,1,0.01)
FDR_path=do.call(rbind, lapply(FDR, function(x) 
  ifelse(abs(x[,1])>1,1,abs(x[,1]))))
FDR_cdf=apply(FDR_path,2,function(x) ecdf(x)(t))
FDR_auc=apply(FDR_path,2, function(x) integrate(ecdf(x),0,1,subdivisions = 500)$value)

pdf("FDR_AUC_lfdr1.pdf")
hist(FDR_auc,breaks=100,freq = F, xlab="FDR AUC", ylab = "Density",main="",
     col = "gray80",xaxs="i",yaxs="i",xaxt="n")
axis(1,at=seq(0,1,.1),labels = seq(0,1,.1))
box()
dev.off()

pdf("FDR_CDF_lfdr2.pdf")
op<-par(mar=c(4,4.1,4.1,1.1),tck=-.01,cex=1)
matplot(t,FDR_cdf[,gene_index],type="l",
        lty=rep(1,83),lwd=rep(2,83),
        col=ifelse(!is.na(match(gene_index,which(FDR_auc>=0.6))),"grey30","grey80"),
        # col=c(rep("grey30",50),rep("grey80",33)),
        xlab="",ylab="",xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,at=seq(0,1,0.1),line=0,padj=-1.25)
axis(2,at=seq(0,1,by=0.2),line=0,hadj=0.5,las=2)
mtext(side=1,text="FDR",line=1.5,font = 1)
mtext(side=2,text="CDF",line=1.75,font = 1)
par(op)
dev.off()




