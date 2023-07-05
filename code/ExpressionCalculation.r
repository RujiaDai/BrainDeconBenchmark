#Computing cell type expression with bMIND, swCAM and TCA

##load bulk-tissue data and cell-pooled reference 
load('rosmap_covcorrected_bulk.RData')
load('props_rosmap.RData')
load('rosmap_meanref.RData')


bulk<-as.matrix(combat_expr)
ref<-as.matrix(meanexp)
sdtmp<-apply(ref,1,sd)
ref<-ref[sdtmp!=0,]
tmp<-intersect(rownames(bulk),rownames(ref))
bulk<-bulk[match(tmp,rownames(bulk)),]
ref<-ref[match(tmp,rownames(ref)),]




#remove cell types with zero proportion
rawprop2<-list()
for(j in 1:length(topprop)){
tmp<-apply(topprop[[j]],2,mean)
rawprop2[[j]]<-topprop[[j]][,tmp!=0]
}
names(rawprop2)<-names(topprop)

load('ros_obj_rnaslot.RData')
scdata<-as.matrix(ros_raw@assays$RNA@counts)
smeta<-ros_raw@meta.data
##bMIND
library(MIND)
bmind<-list()
prior = get_prior(sc = as.matrix(scdata), sample = smeta$projid, cell_type = smeta$broad.cell.type)
profile<-prior$profile
covariance<-prior$covariance
gene1<-intersect(rownames(bulk),rownames(profile))

for(i in 1:length(topprop)){
bmind[[i]] = bMIND(bulk[match(gene1,rownames(bulk)),], frac = as.matrix(rawprop2[[i]]), profile=profile[match(gene1,rownames(profile)),match(colnames(rawprop2[[i]]),colnames(profile))],  covariance= covariance[match(gene1,rownames(covariance)),match(colnames(rawprop2[[i]]),colnames(profile)),match(colnames(rawprop2[[i]]),colnames(profile))],ncore = 16)
}
names(bmind)<-names(rawprop2)
save(bmind,file='bmind_rosmap_withpropfile.RData')



#swcam
swcam<-list()
source('/mnt/compute/Groups/LiuLab/User/Rdai/tools/swCAM/sCAMfastNonNeg.R')

for(j in 1:length(rawprop2)){
print(j)
Aest<-rawprop2[[j]]
Xn<-t(2^(bulk[match(gene1,rownames(bulk)),]))
Sest<-t(2^(profile[match(gene1,rownames(profile)),match(colnames(rawprop2[[j]]),colnames(profile))]))
eta <- 1000
iteradmm <- 1000
rsCAM <- sCAMfastNonNeg(as.matrix(Xn), as.matrix(Aest), as.matrix(Sest), eta = eta, iteradmm=iteradmm, silent = T)
swcam[[j]]<- rsCAM$S
}
names(swcam)<-names(rawprop2)
save(swcam,file='swcam_rosmap_withpropfile.RData')


##TCA
tca<-list() 
library(TCA)
for(j in 1:length(rawprop2)){
tmp <- tca(X = bulk[match(gene1,rownames(bulk)),],
                      W = rawprop2[[j]],
                      )
tca[[j]]<- tensor(bulk[match(gene1,rownames(bulk)),], tmp, verbose=FALSE)}
names(tca)<-names(rawprop2)

save(tca,file='tca_rosmap_withpropfile.RData')


