#prepare phenotype file and covariates file
load('ROSMAP_rawcounts_38704genes.RData')
load('rosmap_s1112_deconvolution.RData')
setwd('rosmap_ieqtl')
mrna<-bulk[,match(rownames(mrmeta),colnames(bulk))]
# 1. normalization 
library(edgeR)
logcpm = cpm(calcNormFactors(DGEList(counts = mrna),method = 'TMM'),log=T)

# 2. Remove lowly expressed
# cutoff=0.1, 25% subjects
keep <- (rowSums(logcpm > .1)) > 0.25*dim(logcpm)[2]
print(table(keep))
logcpm.fgene <- logcpm[keep,]
counts.fgene<-mrna[keep,]

#3. Remove outliers
library(WGCNA)
normadj <- adjacency(logcpm.fgene,type = 'signed',corFnc = 'bicor')   #Calculate network adjacency
netsummary <- fundamentalNetworkConcepts(normadj)
C <- netsummary$Connectivity   #Extract connectivity of each sample
Z.C <- (C-mean(C))/sqrt(var(C))   #Covert to Z-score

pdf('bulk_outlier_plot.pdf')
par(mfrow=c(1,1))
plot(1:length(Z.C),Z.C,main="Outlier Plot",xlab = "Samples",ylab="Connectivity Z Score")
abline(h=(-3), col="red")
outliers <- (Z.C < -3)###can be adjusted 45
print(which(outliers))
dev.off()

logcpm.fgene.fsample <- logcpm.fgene[,!outliers]
counts.fgene.fsample <- counts.fgene[,!outliers]
mrmeta<-mrmeta[!outliers,]#

##4. quantile normalization
library(preprocessCore)
logcpm.fgene.fsample.qn<-normalize.quantiles(as.matrix(logcpm.fgene.fsample),copy=T)  # Quantile normalization across columns
rownames(logcpm.fgene.fsample.qn)<-rownames(logcpm.fgene.fsample)
colnames(logcpm.fgene.fsample.qn)<-colnames(logcpm.fgene.fsample)


mrmeta$age_death<-sub("90[+]",95,mrmeta$age_death)
mrmeta$age_death<-as.numeric(mrmeta$age_death)

PC = prcomp(t(logcpm.fgene.fsample.qn))$x[,1:5] 
##correlation between variables
library(corrplot)
design<-mrmeta[,c("RIN","sequencingBatch","msex","educ","race","age_death","pmi","ceradsc","cogdx")]
design$sequencingBatch<-as.numeric(factor(design$sequencingBatch))

metatmp<-data.matrix(cbind(PC,design))
mrmeta$RIN<-as.numeric(mrmeta$RIN)
mrmeta$ceradsc<-as.factor(mrmeta$ceradsc)
mrmeta$sequencingBatch<-as.factor(mrmeta$sequencingBatch)
mrmeta$race<-as.factor(mrmeta$race)
mrmeta$pmi<-as.numeric(mrmeta$race)
mrmeta$educ<-as.numeric(mrmeta$educ)
library(sva)
data.batch<-mrmeta$sequencingBatch
exprMat <- as.matrix(logcpm.fgene.fsample.qn)

combat_expr <- ComBat(dat = exprMat, batch = data.batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)

genosample<-read.table("ROSMAP_wgsid.csv",header=T,sep=',')
overlap<-intersect(genosample[,1],mrmeta[,1])
mrmeta<-mrmeta[na.omit(match(overlap,mrmeta[,1])),]
mrmeta$wgsid<-genosample[match(overlap,genosample[,1]),2]


info<-read.table("gene.bed",sep=' ',header=F)
info2<-read.table("gene.id",sep='\t',header=F)
#annotate gene location and change gene name to Ensembl id
info[,4]<-info2[,1]
info[,1]<-gsub('^chr','',info[,1])
info<-info[info[,1] %in% c(1:22),]
gene<-intersect(info[,4],rownames(combat_expr))
combat_expr<-combat_expr[match(gene,rownames(combat_expr)),]
combat_expr<-combat_expr[,match(rownames(mrmeta),colnames(combat_expr))]
colnames(combat_expr)<-mrmeta[match(colnames(combat_expr),rownames(mrmeta)),]$wgsid

#covarites
library(peer)
cov<-function(quant,nFactor=0){
if(nFactor==0){
if(ncol(quant) < 150) nFactor = 15 else if (ncol(quant) < 250)
     nFactor = 30 else if (ncol(quant) < 350)
       nFactor = 45 else nFactor = 60}
message("Setting nFactor depend on sample size (GTEx suggestion): ", nFactor)
  expr = t(as.matrix(quant))
  model = PEER()  # create the model object
  invisible(PEER_setPhenoMean(model,expr))# set the observed data
  invisible(PEER_setNk(model,nFactor)) # gradient number of factors
  invisible(PEER_getNk(model))
  invisible(PEER_setAdd_mean(model, TRUE))  # include an additional factor (covariate) to account for the mean expression
  #PEER_setNmax_iterations(model, 100)  # If the model is not converged after the default 1,000 iterations, and the variance of the residuals keeps decreasing, choose a higher value of iterations, e.g., 10,000.
  time = system.time(PEER_update(model))  # perform the inference
  factors = data.frame(PEER_getX(model));  # inferred confounders samples x PEER factors
  factors=factors[-1];
  rownames(factors) = rownames(expr);
  colnames(factors) = paste0("peer",1:ncol(factors))
  weights = PEER_getW(model)  # their weights
  precision = PEER_getAlpha(model)     # precision (inverse variance) of the weights
  residuals = PEER_getResiduals(model) # the residual dataset
  #plot(precision)
  #PEER_plotModel(model)
  #Variance_factor_plot
  Variance = apply(factors,2,var);
  Variance<-sort(Variance, decreasing=TRUE); 
  Variance<-100*Variance/sum(Variance)
  deltavariance = Variance[-length(Variance)] - Variance[-1]
  n<-c()#拐点
  for(i in 2:(length(deltavariance)-1)){
    point<-deltavariance[i]
    if(point>deltavariance[i-1] & point>deltavariance[i+1] & point>(10/nFactor)){
        n<-c(n,i)
    }
  }

  if(length(n)==0){
    n<-c()
    for(i in 2:(length(deltavariance)-1)){
    point<-deltavariance[i]
    if(point>deltavariance[i-1] & point>deltavariance[i+1]& point>0.005){
        n<-c(n,i)
    }} 
  }
  if(length(n)==0){n<-1}
   nFactor<-n[length(n)]+1
   print(nFactor)
   factors<-t(factors)
   data<-factors[1:nFactor,]
   data<-cbind(rownames(data),data)
   colnames(data)[1]<-"Factor"
   write.table(data,"cov.txt", sep="\t", row.names=F, quote=F,col.names=T)
}
cov(combat_expr)

#expression
info[match(gene,info[,4]),]->temp
tmp<-cbind(temp[,c(1:2)],temp[,2]+1,temp[,4])
colnames(tmp)<-c("#chr","start","end","phenotype_id")
data<-as.data.frame(combat_expr)
data<-na.omit(cbind(tmp,data))
print(dim(data))
write.table(data,"pheno.txt",row.names=F,col.names=T,quote=F,sep='\t') 

interaction<-prop[match(rownames(mrmeta),rownames(prop)),]
interaction<-cbind(mrmeta$wgsid,interaction)
colnames(interaction)[1]<-"sample"
write.table(interaction,'interaction.txt',col.names=T,row.names=F,quote=F,sep='\t')
cellname<-colnames(interaction)[-1]
for(i in 1:length(cellname)){
	tmp<-interaction[,c(1,i+1)]
	write.table(tmp,paste0("interaction/",cellname[i],'.txt'),col.names=T,row.names=F,quote=F,sep='\t')
}


info[match(gene,info[,4]),]->temp
tmp<-cbind(temp[,c(1:2)],temp[,2]+1,temp[,4],temp[,4],temp[,5])
colnames(tmp)<-c("#chr","start","end","gene","gene","strand")
data<-as.data.frame(combat_expr)
data<-na.omit(cbind(tmp,data))
print(dim(data))
write.table(data,"../bulk/pheno.txt",row.names=F,col.names=T,quote=F,sep='\t') 


