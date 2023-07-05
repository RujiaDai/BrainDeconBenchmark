#marker genes identification
load("ros_obj_rnaslot.RData")

##cell level: wilcox
load('ros_obj_rnaslot.RData')
Idents(ros_raw)<-ros_raw@meta.data$broad.cell.type
wilcox_ovs<-FindAllMarkers(ros_raw, test.use = "wilcox",min.pct = 0.1,logfc.threshold=0.25,only.pos=TRUE)

t1<-rep(colnames(meanexp),rep(ncol(meanexp)-1,ncol(meanexp)))
t2<-c()
for(x in 1:ncol(meanexp)){t2<-c(t2,colnames(meanexp)[-x])}



library(future)
plan('multiprocess',workers=8)
#options(future.globals.maxSize = 8000 * 1024^2)
wilcox_2nd<-list()
for(i in 1:length(t1)){
print(i)
wilcox_2nd[[i]]<-FindMarkers(ros_raw,ident.1=t1[[i]],ident.2=t2[[i]], test.use = "wilcox",min.pct = 0.1,logfc.threshold=0.25,only.pos=TRUE)
}




##pseudobulk level: deseq2

#construte count matrix of pseudobulk
pscounts<-list()
tmp<-unique(smeta$class)
for(i in 1:length(tmp)){
cs<-scdata[,smeta$class==tmp[[i]]]
csmeta<-subset(smeta,class==tmp[[i]])
write.csv(csmeta,'tmpfile.csv')
csmeta<-read.csv('tmpfile.csv',header=T,sep=',')
pscounts[[i]]<-t(apply(cs,1,function(x){tapply(x,csmeta$individualID,sum)}))
}
 names(pscounts)<-tmp
#load('rosmap_pscounts.RData')


psmatrix<-Reduce(cbind,pscounts)
celltype<-rep(names(pscounts),c(rep(47,6),36,36))
sample<-c(rep(colnames(pscounts[[1]]),6),rep(colnames(pscounts[[7]]),2))
colData<-data.frame(sample=sample,celltype=celltype)
rownames(colData)<-colnames(psmatrix)<-paste0(colData$sample,"_",colData$celltype)
colmeta<-colData

#identification of marker gene with DESeq2
deseq_ovs<-list()
cellid<-unique(wilcox_ovs$cluster)
for(i in 1:length(cellid)){
exp<-psmatrix+1
colData<-colmeta
colnames(colData)[[2]]<-'group'
colData$group<-ifelse(colData$group==cellid[[i]],cellid[[i]],'others')
colData$group<-factor(colData$group,levels=c('others',cellid[[i]]))

dds <- DESeqDataSetFromMatrix(countData = exp,
                              colData = colData,
                              design= ~ group)
dds <- DESeq(dds,test='LRT',reduced=~1)
resultsNames(dds) # lists the coefficients
res <- results(dds)
deseq_ovs[[i]]<-res
deseq_ovs[[i]]$cluster<-cellid[[i]]
}



deseq_2nd<-list()
for(i in 1:length(t1)){
exp<-(psmatrix+1)[,colmeta$celtype==t1[[i]]|colmeta$celltype==t2[[i]]]
colData<-subset(colmeta,celltype==t1[[i]]|celltype==t2[[i]])
colnames(colData)[[2]]<-'group'
colData$group<-factor(colData$group,levels=c(t1[[i]],t2[[i]]))

dds <- DESeqDataSetFromMatrix(countData = exp,
                              colData = colData,
                              design= ~ group)
dds <- DESeq(dds,test='LRT',reduced=~1)
resultsNames(dds) # lists the coefficients
res <- results(dds)
deseq_2nd[[i]]<-res
deseq_2nd[[i]]$cluster<-cellid[[i]]
}

for(i in 1:length(t1)){
wilcox_2nd[[i]]$maxid<-t1[[i]]
wilcox_2nd[[i]]$max2id<-t2[[i]]
deseq_2nd[[i]]$maxid<-t1[[i]]
deseq_2nd[[i]]$max2id<-t2[[i]]
}
wilcox_ovs<-subset(wilcox_ovs,p_val_adj<0.05)
wilcox_ovs<-wilcox_ovs[order(wilcox_ovs$avg_log2FC,decreasing=T),]
wilcox_ovs<-wilcox_ovs[!duplicated(wilcox_ovs$gene),]

for(i in 1:length(deseq_ovs)){
deseq_ovs[[i]]$gene<-rownames(deseq_ovs[[i]])
deseq_ovs[[i]]<-subset(deseq_ovs[[i]],padj<0.05&log2FoldChange>1.8)
}
deseq_ovs<-Reduce(rbind,deseq_ovs)
deseq_ovs<-deseq_ovs[order(deseq_ovs,decreasing=T),]
deseq_ovs<-deseq_ovs[!duplicated(deseq_ovs$gene),]

for(i in 1:length(wilcox_2nd)){
wilcox_2nd[[i]]<-subset(wilcox_2nd[[i]],p_val_adj<0.05)
wilcox_2nd[[i]]$gene<-rownames(wilcox_2nd[[i]])
}
wilcox_2nd<-Reduce(rbind,wilcox_2nd)
wilcox_2nd<-wilcox_2nd[order(wilcox_2nd$avg_log2FC,decreasing=T),]
wilcox_2nd<-wilcox_2nd[!duplicated(wilcox_2nd$gene),]


for(i in 1:length(deseq_2nd)){
deseq_2nd[[i]]<-subset(deseq_2nd[[i]],padj<0.05&log2FoldChange>1.8)
deseq_2nd[[i]]$gene<-rownames(deseq_2nd[[i]])
}
deseq_2nd<-Reduce(rbind,deseq_2nd)
deseq_2nd<-deseq_2nd[order(deseq_2nd$log2FoldChange,decreasing=T),]
deseq_2nd<-deseq_2nd[!duplicated(deseq_2nd$gene),]

save(wilcox_2nd,deseq_2nd,wilcox_ovs,deseq_ovs,file='rosmap_marker.RData')


mg1<-subset(deseq_2nd,log2FoldChange>2)
mg2<-subset(wilcox_2nd,avg_log2FC>1)
mg1$cid<-paste(mg1$gene,mg1$maxid)
mg2$cid<-paste(mg2$gene,mg2$maxid)

mg11<-mg1[match(intersect(mg1$cid,mg2$cid),mg1$cid),]
mg22<-mg2[match(intersect(mg1$cid,mg2$cid),mg2$cid),]


write.csv(mg11,'mg_rosmap.csv')