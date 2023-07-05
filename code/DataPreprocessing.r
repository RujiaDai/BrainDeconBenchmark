##Bulk rnaseq data propcessing 
##ROSMAP 
setwd('F:/deconbench/rosmap/raw/bulk')
load('ROSMAP_rawcounts_38704genes.RData')##raw count matrix and metadata
tmp<-read.table('rosmap_sn_ihc_sample.txt',header=T,sep='\t')



mrna<-bulk[,bulkmeta$projid%in%tmp$projID]
mrmeta<-bulkmeta[bulkmeta$projid%in%tmp$projID,]

mrna<-mrna[,mrmeta$region=='DLPFC']#select DLPFC samples
mrmeta<-subset(mrmeta,region=='DLPFC')

#samples having IHC or snRNAseq data were used in method evaluation
#all DLPFC samples were used in decon-eQTL mapping

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

# optional
#remove samples with missing information 
# logcpm.fgene.fsample<-logcpm.fgene.fsample[,!is.na(mrmeta$pmi)]
# mrmeta<-mrmeta[!is.na(mrmeta$pmi),]


##4. quantile normalization
library(preprocessCore)
logcpm.fgene.fsample.qn<-normalize.quantiles(as.matrix(logcpm.fgene.fsample),copy=T)  # Quantile normalization across columns
rownames(logcpm.fgene.fsample.qn)<-rownames(logcpm.fgene.fsample)
colnames(logcpm.fgene.fsample.qn)<-colnames(logcpm.fgene.fsample)

save(logcpm.fgene.fsample,counts.fgene.fsample,mrmeta,logcpm.fgene.fsample.qn,file='rosmap_qc_bulk.RData')

mrmeta$age_death<-sub("90[+]",95,mrmeta$age_death)
mrmeta$age_death<-as.numeric(mrmeta$age_death)

PC = prcomp(t(logcpm.fgene.fsample.qn))$x[,1:5] 
##correlation between variables
library(corrplot)
design<-mrmeta[,c("RIN","sequencingBatch","msex","educ","race","age_death","pmi","ceradsc","cogdx")]
design$sequencingBatch<-as.numeric(factor(design$sequencingBatch))

metatmp<-data.matrix(cbind(PC,design))

metacor<-cor(metatmp,method='spearman')
write.csv(metacor,'rosmap_bulk_correlation_pc_meta.csv')
pdf('rosmap_bulk_correlation_pc_meta.pdf',width=8,height=8)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(metacor, method="color", col=col(200),  
         type="upper", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
         )
dev.off()


##pca plot before covariate correction
mrmeta$RIN<-as.numeric(mrmeta$RIN)
mrmeta$ceradsc<-as.factor(mrmeta$ceradsc)
mrmeta$sequencingBatch<-as.factor(mrmeta$sequencingBatch)
mrmeta$race<-as.factor(mrmeta$race)
mrmeta$pmi<-as.numeric(mrmeta$pmi)
mrmeta$educ<-as.numeric(mrmeta$educ)

library(ggfortify)
library(gridExtra)
pc<-prcomp(t(logcpm.fgene.fsample.qn))
g1<-autoplot(pc, data = mrmeta, colour = 'ceradsc')
g2<-autoplot(pc, data = mrmeta, colour = 'sequencingBatch')
g3<-autoplot(pc, data = mrmeta, colour = 'age_death')
g4<-autoplot(pc, data = mrmeta, colour = 'msex')
g5<-autoplot(pc, data = mrmeta, colour = 'educ')
g6<-autoplot(pc, data = mrmeta, colour = 'race')
g7<-autoplot(pc, data = mrmeta, colour = 'RIN')
g8<-autoplot(pc, data = mrmeta, colour = 'pmi')
g10<-autoplot(pc, data = mrmeta, colour = 'cogdx')
pdf("rosmap.bulk.qn.pca.biological.pdf",width=10,height=6)
grid.arrange(g1,g10,g3,g4,g5,g6,nrow=2)
dev.off()


pdf("rosmap.bulk.qn.pca.technical.pdf",width=8,height=6)
grid.arrange(g2,g7,g8,nrow=2)
dev.off()



# 5.remove batch effect with ComBat function
library(sva)
data.batch<-mrmeta$sequencingBatch
exprMat <- as.matrix(logcpm.fgene.fsample.qn)
combat_expr <- ComBat(dat = exprMat, batch = data.batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)


# 6.other covariates correction (linear regression)
mrmeta$ceradsc<-as.numeric(mrmeta$ceradsc)
mrmeta$pmi<-as.numeric(mrmeta$pmi)
mrmeta$educ<-as.numeric(mrmeta$educ)

mod1 = model.matrix(~ceradsc+msex+age_death+educ+RIN+pmi,data=mrmeta)
Y = as.matrix(combat_expr)
X = as.matrix(mod1)
beta = as.matrix((solve(t(X)%*%X)%*%t(X))%*%t(Y))

bulk.techRegressed = Y - t(X[,c(6:ncol(X))] %*% beta[c(6:nrow(beta)),])#regress technical covariates
bulk.AllRegressed = Y - t(X[,c(2:ncol(X))] %*% beta[c(2:nrow(beta)),])#regress all covariates

save(combat_expr,bulk.techRegressed,bulk.AllRegressed,mrmeta,file='rosmap_covcorrected_bulk.RData')
#7. post-processing PCA
library(ggfortify)
library(gridExtra)
pc<-prcomp(t(bulk.techRegressed))
g1<-autoplot(pc, data = mrmeta, colour = 'ceradsc')
g2<-autoplot(pc, data = mrmeta, colour = 'sequencingBatch')
g3<-autoplot(pc, data = mrmeta, colour = 'age_death')
g4<-autoplot(pc, data = mrmeta, colour = 'msex')
g5<-autoplot(pc, data = mrmeta, colour = 'educ')
g6<-autoplot(pc, data = mrmeta, colour = 'race')
g7<-autoplot(pc, data = mrmeta, colour = 'RIN')
g8<-autoplot(pc, data = mrmeta, colour = 'pmi')
g10<-autoplot(pc, data = mrmeta, colour = 'cogdx')
pdf("rosmap.bulk.regtech.pca.biological.pdf",width=10,height=6)
grid.arrange(g1,g10,g3,g4,g5,g6,nrow=2)
dev.off()


pdf("rosmap.bulk.regtech.pca.technical.pdf",width=8,height=6)
grid.arrange(g2,g7,g8,nrow=2)
dev.off()


