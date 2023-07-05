#Evaluation of deconvoluted expression by bMIND, swCAM and TCA	
load('tca_rosmap_withpropfile.RData')
load('swcam_rosmap_withpropfile.RData')
load('bmind_rosmap_withpropfile.RData')

load('ros_obj_rnaslot.RData')#load ground truth
scdata<-as.matrix(ros_raw@assays$RNA@counts)
smeta<-ros_raw@meta.data
scdata<-scdata[,smeta$projid%in%mrmeta$projid]
smeta<-smeta[smeta$projid%in%mrmeta$projid,]

library(edgeR)
truth<-c()
cellname<-colnames(meanexp)
for(i in 1:ncol(meanexp)){
metatmp<-subset(smeta,broad.cell.type==cellname[[i]])
sctmp<-scdata[,smeta$broad.cell.type==cellname[[i]]]
truth[[i]]<-cpm(t(apply(sctmp,1,function(x){tapply(x,metatmp$projid,mean)})),log=F)
}
names(truth)<-colnames(meanexp)

#Calculate correlation between deconvoluted expression and ground truth 
bmindcor_gene<-c()
bmindcor_sample<-c()
for(i in 1:length(bmind)){
sampletmp<-c()
genetmp<-c()
for(j in 1:dim(bmind[[i]]$A)[[2]]){
test1<-bmind[[i]]$A[,j,]
colnames(test1)<-mrmeta$projid
cid<-colnames(bmind[[i]]$A)[[j]]
test2<-truth[[ match(colnames(bmind[[i]]$A)[[j]],names(truth))]]
tmp1<-intersect(rownames(test1),rownames(test2))
tmp2<-intersect(colnames(test1),colnames(test2))
sampletmp[[j]]<-diag(cor(test1[match(tmp1,rownames(test1)),match(tmp2,colnames(test1))],log2(1+test2[match(tmp1,rownames(test2)),match(tmp2,colnames(test2))]),method='spearman'))
genetmp[[j]]<-diag(cor(t(test1[match(tmp1,rownames(test1)),match(tmp2,colnames(test1))]),t(log2(1+test2[match(tmp1,rownames(test2)),match(tmp2,colnames(test2))])),method='spearman'))
}
bmindcor_sample[[i]]<-Reduce(cbind,sampletmp)
bmindcor_gene[[i]]<-Reduce(cbind,genetmp)
colnames(bmindcor_sample[[i]])<-colnames(bmindcor_gene[[i]])<-colnames(bmind[[i]]$A[,,1])
}


swcamcor_sample<-c()
swcamcor_gene<-c()
for(z in 1:length(swcam)){
genetmp<-c()
sampletmp<-c()
cid<-colnames(rawprop2[[z]])
for(j in 1:dim(swcam[[z]])[[1]]){
test1<-swcam[[z]][j,,]
rownames(test1)<-rownames(bmind[[1]]$A)
colnames(test1)<-mrmeta$projid
test2<-truth[[match(cid[[j]],names(truth))]]
tmp1<-intersect(rownames(test1),rownames(test2))
tmp2<-intersect(colnames(test1),colnames(test2))
sampletmp[[j]]<-diag(cor(log2(1+test1[match(tmp1,rownames(test1)),match(tmp2,colnames(test1))]),log2(1+test2[match(tmp1,rownames(test2)),match(tmp2,colnames(test2))]),method='spearman'))
genetmp[[j]]<-diag(cor(t(log2(1+test1[match(tmp1,rownames(test1)),match(tmp2,colnames(test1))])),t(log2(1+test2[match(tmp1,rownames(test2)),match(tmp2,colnames(test2))])),method='spearman'))

}
swcamcor_sample[[z]]<-Reduce(cbind,sampletmp)
swcamcor_gene[[z]]<-Reduce(cbind,genetmp)
colnames(swcamcor_sample[[z]])<-colnames(swcamcor_gene[[z]])<-cid
}



tcacor_gene<-c()
tcacor_sample<-c()
for(z in 1:length(tca)){
cid<-colnames(rawprop2[[z]])
genetmp<-c()
sampletmp<-c()
for(j in 1:length(tca[[z]])){
test1<-tca[[z]][[j]] 
colnames(test1)<-mrmeta$projid
test2<-truth[[match(cid[[j]],names(truth))]]
tmp1<-intersect(rownames(test1),rownames(test2))
tmp2<-intersect(colnames(test1),colnames(test2))
sampletmp[[j]]<-diag(cor(test1[match(tmp1,rownames(test1)),match(tmp2,colnames(test1))],log2(1+test2[match(tmp1,rownames(test2)),match(tmp2,colnames(test2))]),method='spearman'))
genetmp[[j]]<-diag(cor(t(test1[match(tmp1,rownames(test1)),match(tmp2,colnames(test1))]),t(log2(1+test2[match(tmp1,rownames(test2)),match(tmp2,colnames(test2))])),method='spearman'))

}
tcacor_sample[[z]]<-Reduce(cbind,sampletmp)
tcacor_gene[[z]]<-Reduce(cbind,genetmp)
colnames(tcacor_sample[[z]])<-colnames(tcacor_gene[[z]])<-cid
}
names(bmindcor_sample)<-names(bmindcor_gene)<-names(swcamcor_gene)<-names(swcamcor_sample)<-names(tcacor_gene)<-names(tcacor_sample)<-names(rawprop2)

#reformate correlations to tables for plotting
list_to_table<-function(x,expname){
tmp<-lapply(x,melt)
for(i in 1:length(tmp)){
tmp[[i]]$prop<-names(x)[[i]]
tmp[[i]]$exp<-expname
}
Reduce(rbind,tmp)
}

bmind_sample<-list_to_table(bmindcor_sample,'bmind')
swcam_sample<-list_to_table(swcamcor_sample,'swcam')
tca_sample<-list_to_table(tcacor_sample,'tca')

bmind_gene<-list_to_table(bmindcor_gene,'bmind')
swcam_gene<-list_to_table(swcamcor_gene,'swcam')
tca_gene<-list_to_table(tcacor_gene,'tca')

genecor<-rbind(bmind_gene,swcam_gene,tca_gene)
samplecor<-rbind(bmind_sample,swcam_sample,tca_sample)

write.csv(genecor,'rosmap_exp_genecor_withcovariance.csv')
write.csv(samplecor,'rosmap_exp_samplecor_withcovariance.csv')

