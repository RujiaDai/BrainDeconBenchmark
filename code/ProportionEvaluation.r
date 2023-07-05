##Evaluation of cell porprotion estimated deconvolution methods

##load ground truth
ihc<-read.csv('raw/ihc/ihcprop.csv',header=T,sep=',',row.names=1)
colnames(ihc)<-c('Ast','End','Mic','Neu','Oli')
scprop<-read.csv('raw/ihc/cell_counts_sc.csv',header=T,sep=',',row.names=1)
scprop<-t(apply(scprop,1,function(x){x/sum(x)}))
scprop<-scprop[,match(colnames(topprop[[1]]),colnames(scprop))]
scprop<-scprop[rownames(scprop)%in%mrmeta$projid,]


##calculate RMSE and Pearson correlation at sample level (using scprop as groud truth)
library(Metrics)
rmse_sample_scprop<-list()
cor_sample_scprop<-list()
for(j in 1:length(topprop)){
deconprop<-topprop[[j]][match(rownames(scprop),mrmeta$projid),]
rownames(deconprop)<-rownames(scprop)
tmp<-as.matrix(cbind(scprop,deconprop))
rmse_sample_scprop[[j]]<-apply(tmp,1,function(x){rmse(as.numeric(x[1:8]),as.numeric(x[9:16]))})
cor_sample_scprop[[j]]<-apply(tmp,1,function(x){cor(as.numeric(x[1:8]),as.numeric(x[9:16]),method='pearson')})
}
names(rmse_sample_scprop)<-names(topprop)
names(cor_sample_scprop)<-names(topprop)


cor_cell_scprop<-list()
rmse_cell_scprop<-list()
for(j in 1:length(topprop)){
deconprop<-topprop[[j]][match(rownames(scprop),mrmeta$projid),]
rownames(deconprop)<-rownames(scprop)
tmp<-as.matrix(rbind(scprop,deconprop))
rmse_cell_scprop[[j]]<-apply(tmp,2,function(x){rmse(as.numeric(x[1:35]),as.numeric(x[36:70]))})
cor_cell_scprop[[j]]<-apply(tmp,2,function(x){cor(as.numeric(x[1:35]),as.numeric(x[36:70]),method='pearson')})
}
names(rmse_cell_scprop)<-names(topprop)
names(cor_cell_scprop)<-names(topprop)



##calculate RMSE and Pearson correlation at sample level (using IHC proprotions as groud truth)
rmse_sample_ihc<-list()
cor_sample_ihc<-list()
for(j in 1:length(topprop)){
deconprop<-data.frame(topprop[[j]][match(rownames(ihc),mrmeta$projid),])
deconprop$Neu<-deconprop$Ex+deconprop$In
deconprop<-deconprop[,c('Ast','End','Mic','Neu','Oli')]
deconprop<-t(apply(deconprop,1,function(x){x/sum(x)}))
colnames(deconprop)<-colnames(ihc)
rownames(deconprop)<-rownames(ihc)
tmp<-as.matrix(cbind(ihc,deconprop))
rmse_sample_ihc[[j]]<-apply(tmp,1,function(x){rmse(as.numeric(x[1:5]),as.numeric(x[6:10]))})
cor_sample_ihc[[j]]<-apply(tmp,1,function(x){cor(as.numeric(x[1:5]),as.numeric(x[6:10]),method='pearson')})
}
names(rmse_sample_ihc)<-names(topprop)
names(cor_sample_ihc)<-names(topprop)



rmse_cell_ihc<-list()
cor_cell_ihc<-list()
for(j in 1:length(topprop)){
deconprop<-data.frame(topprop[[j]][match(rownames(ihc),mrmeta$projid),])
deconprop$Neu<-deconprop$Ex+deconprop$In
deconprop<-deconprop[,c('Ast','End','Mic','Neu','Oli')]
deconprop<-t(apply(deconprop,1,function(x){x/sum(x)}))
rownames(deconprop)<-rownames(ihc)
tmp<-as.matrix(rbind(ihc,deconprop))
cor_cell_ihc[[j]]<-apply(tmp,2,function(x){cor(as.numeric(x[1:35]),as.numeric(x[36:70]),method='pearson')})
rmse_cell_ihc[[j]]<-apply(tmp,2,function(x){rmse(as.numeric(x[1:35]),as.numeric(x[36:70]))})

}
names(rmse_cell_ihc)<-names(topprop)
names(cor_cell_ihc)<-names(topprop)




##Reformate calculated RMSE and Pearson correlation
rmse_cell_ihc<-Reduce(rbind,rmse_cell_ihc)
rmse_cell_scprop<-Reduce(rbind,rmse_cell_scprop)
rownames(rmse_cell_ihc)<-rownames(rmse_cell_scprop)<-names(topprop)
cor_cell_ihc<-Reduce(rbind,cor_cell_ihc)
cor_cell_scprop<-Reduce(rbind,cor_cell_scprop)
rownames(cor_cell_ihc)<-rownames(cor_cell_scprop)<-names(topprop)

rmse_sample_ihc<-Reduce(cbind,rmse_sample_ihc)
rmse_sample_scprop<-Reduce(cbind,rmse_sample_scprop)
colnames(rmse_sample_ihc)<-colnames(rmse_sample_scprop)<-names(topprop)
cor_sample_ihc<-Reduce(cbind,cor_sample_ihc)
cor_sample_scprop<-Reduce(cbind,cor_sample_scprop)
colnames(cor_sample_ihc)<-colnames(cor_sample_scprop)<-names(topprop)


##Plotting
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggsci)
rmseplot<-function(x,y,z){
if(z=='cell'){
meanihc<-apply(ihc,2,mean)
x<-t(apply(x,1,function(a){a/rank(meanihc)}))
scpropmean<-apply(scprop,2,mean)
y<-t(apply(y,1,function(a){a/rank(scpropmean)}))
test1<-melt(x)
test2<-melt(y)
mean1<-apply(x,1,mean)
mean2<-apply(y,1,mean)
test1$Var1<-factor(test1$Var1,levels=names(mean1)[order(mean1)])
test2$Var1<-factor(test2$Var1,levels=names(mean2)[order(mean2)])

mean3<-apply(ihc,2,mean)
mean4<-apply(scprop,2,mean)
test1$Var2<-factor(test1$Var2,levels=names(mean3)[order(mean3,decreasing=T)])
test2$Var2<-factor(test2$Var2,levels=names(mean4)[order(mean4,decreasing=T)])

g1<-ggplot(test1,aes(x=Var1,y=value,fill=Var2))+geom_histogram(stat='identity',position='dodge')+theme_light(base_size=15)+labs(x='',y='RMSE/prop.rank',title='IHC')+ scale_fill_npg()
g2<-ggplot(test2,aes(x=Var1,y=value,fill=Var2))+geom_histogram(stat='identity',position='dodge')+theme_light(base_size=15)+labs(x='',y='RMSE/prop.rank',title='snRNAseq') + scale_fill_lancet()
grid.arrange(g1,g2)}

if(z=='sample'){
test1<-melt(x)
test2<-melt(y)
mean1<-apply(x,2,mean)
mean2<-apply(y,2,mean)
test1$Var2<-factor(test1$Var2,levels=names(mean1)[order(mean1)])
test2$Var2<-factor(test2$Var2,levels=names(mean2)[order(mean2)])

g1<-ggplot(test1,aes(x=Var2,y=value))+geom_boxplot()+theme_light(base_size=15)+labs(x='',y='RMSE',title='IHC')+ylim(0,0.3)
g2<-ggplot(test2,aes(x=Var2,y=value))+geom_boxplot()+theme_light(base_size=15)+labs(x='',y='RMSE',title='snRNAseq') +ylim(0,0.3)
grid.arrange(g1,g2)}
}
rmseplot(rmse_cell_ihc,rmse_cell_scprop,'cell')
rmseplot(rmse_sample_ihc,rmse_sample_scprop,'sample')



corplot<-function(x,y,z){

if(z =='cell'){

test1<-melt(x)
test2<-melt(y)
mean1<-apply(x,1,mean)
mean2<-apply(y,1,mean)
test1$Var1<-factor(test1$Var1,levels=names(mean1)[order(mean1)])
test2$Var1<-factor(test2$Var1,levels=names(mean2)[order(mean2)])

mean3<-apply(ihc,2,mean)
mean4<-apply(scprop,2,mean)
test1$Var2<-factor(test1$Var2,levels=names(mean3)[order(mean3,decreasing=F)])
test2$Var2<-factor(test2$Var2,levels=names(mean4)[order(mean4,decreasing=F)])

g1<-ggplot(test1,aes(x=Var1,y=value,fill=Var2))+geom_histogram(stat='identity',position='dodge')+theme_light(base_size=15)+labs(x='',y='Correlation coefficient',title='IHC')+ scale_fill_npg()
g2<-ggplot(test2,aes(x=Var1,y=value,fill=Var2))+geom_histogram(stat='identity',position='dodge')+theme_light(base_size=15)+labs(x='',y='Correlation coefficient',title='snRNAseq') + scale_fill_lancet()
grid.arrange(g1,g2)}

if(z=='sample'){
test1<-melt(x)
test2<-melt(y)
mean1<-apply(x,2,mean)
mean2<-apply(y,2,mean)
test1$Var2<-factor(test1$Var2,levels=names(mean1)[order(mean1,decreasing=T)])
test2$Var2<-factor(test2$Var2,levels=names(mean2)[order(mean2,decreasing=T)])

g1<-ggplot(test1,aes(x=Var2,y=value))+geom_boxplot()+theme_light(base_size=15)+labs(x='',y='Correlation coefficient',title='IHC')
g2<-ggplot(test2,aes(x=Var2,y=value))+geom_boxplot()+theme_light(base_size=15)+labs(x='',y='Correlation coefficient',title='snRNAseq') 
grid.arrange(g1,g2)}


}

corplot(cor_cell_ihc,cor_cell_scprop,'cell')
corplot(cor_sample_ihc,cor_sample_scprop,'sample')


write.csv(rmse_cell_ihc,'new/rmse_cell_ihc.csv')
write.csv(rmse_cell_scprop,'new/rmse_cell_scprop.csv')
write.csv(rmse_sample_ihc,'new/rmse_sample_ihc.csv')
write.csv(rmse_sample_scprop,'new/rmse_sample_scprop.csv')

write.csv(cor_cell_ihc,'new/cor_cell_ihc.csv')
write.csv(cor_cell_scprop,'new/cor_cell_scprop.csv')
write.csv(cor_sample_ihc,'new/cor_sample_ihc.csv')
write.csv(cor_sample_scprop,'new/cor_sample_scprop.csv')


rmse_sample<-data.frame(rbind(rmse_sample_ihc,rmse_sample_scprop))
rmse_sample$groundtruth<-rep(c('IHC','snRNAseq'),c(nrow(rmse_sample_ihc),nrow(rmse_sample_scprop)))
rmse_sample<-melt(rmse_sample,id='groundtruth')
ggplot(rmse_sample,aes(x=variable,y=value,fill=groundtruth))+geom_boxplot()+labs(x='',y='RMSE')+theme_light(base_size=15)

cor_sample<-data.frame(rbind(cor_sample_ihc,cor_sample_scprop))
cor_sample$groundtruth<-rep(c('IHC','snRNAseq'),c(nrow(cor_sample_ihc),nrow(cor_sample_scprop)))
cor_sample<-melt(cor_sample,id='groundtruth')
ggplot(cor_sample,aes(x=variable,y=value,fill=groundtruth))+geom_boxplot()+labs(x='',y='cor')


