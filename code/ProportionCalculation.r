setwd('F:/deconbench/rosmap')
load('raw/bulk/rosmap_covcorrected_bulk.RData')#bulk-tissue expression
load('new/rosmap_meanref.RData')#reference
mg<-read.csv('new/mg_rosmap.csv',header=T,sep=',',row.names=1)#mg: edgerqlf, wilcox


#intersection of genes in bulk-tissue expression and reference
geneid<-intersect(rownames(meanexp),rownames(combat_expr))
length(geneid)#13488

bulk<-2^(combat_expr[match(geneid,rownames(combat_expr)),])
ref<-2^meanexp[match(geneid,rownames(meanexp)),]

#retain marker genes for deconvolution
bulk<-bulk[rownames(bulk)%in%mg$gene,]
ref<-ref[rownames(ref)%in%mg$gene,]





####deconvolution methods
#dtangle
dtangle_decon<-function(bulk,ref){
library('dtangle')
mixture_samples = t(log2(bulk+1))
reference_samples = t(log2(ref+1))
dtangleprop = dtangle::dtangle(Y=mixture_samples, reference=reference_samples)$estimates
}

#dsa
dsa_decon<-function(bulk,ref,mg){
require(CellMix)
md = mg
ML = CellMix::MarkerList()
ML@.Data <- tapply(as.character(rownames(mg)),as.character(mg$maxid),list)
dsaprop = CellMix::ged(as.matrix(bulk), ML, method = "DSA", log = FALSE)@fit@H
dsaprop = apply(dsaprop,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
dsaprop = t(apply(dsaprop,2,function(x) x/sum(x))) #explicit SbulkO constraint
#write.csv(t(dsaprop),'rosmap_regtechcpm_dasprop.csv')
}

#ols
ols_decon<-function(bulk,ref){
olsprop = apply(bulk,2,function(x) lm(x ~ as.matrix(ref))$coefficients[-1])
olsprop = apply(olsprop,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
olsprop = apply(olsprop,2,function(x) x/sum(x)) #explicit SbulkO constraint
rownames(olsprop) <- unlist(lapply(strsplit(rownames(olsprop),")"),function(x) x[2]))
t(olsprop)
}

#cibersort
source('D:/CIBERSORT.R')

cbprop<-CIBERSORT(ref,bulk, perm=0, QN=FALSE, absolute=FALSE, abs_method='sig.score')[,1:ncol(ref)]
dsaprop<-dsa_decon(bulk,ref,mg)
dtprop<-dtangle_decon(bulk,ref)
olsprop<-ols_decon(bulk,ref)


#prepare data from deconvolution using single-cell reference
counts_bulk <- as.matrix(2^combat_expr)
scdata <- ros_raw@assays$RNA@counts#single-cell refernce
smeta <- ros_raw@meta.data

#remove genes with no variation across cell types
tmp<-apply(scdata,1,sd)
scdata<-scdata[tmp!=0,]

geneid<-intersect(rownames(counts_bulk),rownames(scdata))
bulk<-counts_bulk[match(geneid,rownames(counts_bulk)),]
scmatrix<-scdata[match(geneid,rownames(scdata)),]

 

#bisque and music
library(MuSiC)
library(BisqueRNA)
library(Biobase)
library(Matrix)

sample.ids=colnames(scmatrix)
sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=sample.ids,
                       SubjectName=smeta$projid,
                       cellType=smeta$broad.cell.type)

sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))

sc.pdata <- new("AnnotatedDataFrame",
                data=sc.pheno,
                varMetadata=sc.meta)

sc.eset <- Biobase::ExpressionSet(assayData=as.matrix(scmatrix),
                                  phenoData=sc.pdata)

bulk.eset <- Biobase::ExpressionSet(assayData = bulk)

 
mg4decon<-list()
mgx<-mg11
cellnames<-unique(mgx$maxid)
for(j in 1:length(cellnames)){mg4decon[[j]]<-subset(mgx,maxid==cellnames[[j]])$gene}
names(mg4decon)<-cellnames
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=mg4decon, use.overlap=FALSE)
bisqueprop<-t(res$bulk.props)
Est.prop= music_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, clusters = 'cellType',
                               samples = 'SubjectName', markers=mg4decon, verbose = F)
musicprop<-Est.prop$Est.prop.weighted[,cellnames]

topprop<-list(dtprop,dsaprop,olsprop,cbprop,bisqueprop,musicprop)
names(topprop)<-c('dtangle','dsa','ols','cibersort','bisque','music')

for(i in 1:length(topprop)){
topprop[[i]]<-topprop[[i]][,match(colnames(topprop[[1]]),colnames(topprop[[i]]))]
}
save(topprop,file='new/props_rosmap.RData')

