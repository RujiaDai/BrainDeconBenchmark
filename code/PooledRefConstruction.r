#load snRNAseq count matrix (ROSMAP)
library(Matrix)
genenames<-read.delim('filtered_gene_row_names.txt',header=F)
scdata<-readMM('filtered_count_matrix.mtx')
smeta<-read.table('filtered_column_metadata.txt',header=T,sep='\t',row.names=1)
scdata<-as.matrix(scdata)
rownames(scdata)<-genenames[,1]
colnames(scdata)<-rownames(smeta)

clmeta<-read.csv('snmeta_clinical.csv',header=T,sep=',')
rownames(clmeta)<-rownames(smeta)

smeta$class<-smeta$broad.cell.type
smeta$class<-gsub('Ex','Neuron',smeta$class)
smeta$class<-gsub('In','Neuron',smeta$class)

#clinical data, replace age 90+ with 95
clmeta$age_death<-gsub('90[+]',95,clmeta$age_death)
clmeta$age_death<-as.numeric(clmeta$age_death)
design<-data.matrix(clmeta[,c('msex','educ','race','age_death','ceradsc','cogdx','pmi')])
smeta<-cbind(smeta,design)

#remove missing data
scdata<-scdata[,!is.na(smeta$pmi)]
smeta<-smeta[!is.na(smeta$pmi),]


#create seurat object
#data is already qc-ed
library(Seurat)
ros_raw <- CreateSeuratObject(counts = scdata, meta=smeta, project = "ros_rawmap", min.cells = 0, min.features = 0)
ros_raw[["percent.mt"]] <- PercentageFeatureSet(ros_raw, pattern = "^MT-")
ros_raw <- NormalizeData(object = ros_raw, normalization.method = "LogNormalize", scale.factor = 10000)
ros_raw <- FindVariableFeatures(ros_raw, selection.method = "vst", nfeatures = 2000)
ros_raw <- ScaleData(ros_raw,features = rownames(ros_raw))
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ros_raw  <- CellCycleScoring(ros_raw , s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
ros_raw <- RunPCA(ros_raw, features = VariableFeatures(object = ros_raw))
ros_raw <- RunUMAP(ros_raw, dims = 1:20)

ros_reg_tech <- ScaleData(ros_raw ,vars.to.regress=c('nFeature_RNA','projid','pmi','S.Score', 'G2M.Score'),features = rownames(ros_raw))
ros_reg_tech <- RunPCA(ros_reg_tech, features = VariableFeatures(object = ros_reg_tech))
ros_reg_tech <- RunUMAP(ros_reg_tech, dims = 1:20)


ros_reg_all <- ScaleData(ros_raw, vars.to.regress=c('age_death','msex','ceradsc','race','educ','nFeature_RNA','projid','pmi','S.Score', 'G2M.Score'),features = rownames(ros_raw))
ros_reg_all <- RunPCA(ros_reg_all, features = VariableFeatures(object = ros_reg_all))
ros_reg_all <- RunUMAP(ros_reg_all, dims = 1:20)

save(ros_raw,ros_reg_all,ros_reg_tech,file='ros_obj_rnaslot.RData')


#pooled ref construction 
library(edgeR)
meancounts<-AverageExpression(ros_raw,slot='counts')[[1]]
meanexp<-log2(1+cpm(meancounts,log=F))
save(meancounts,meanexp,file='rosmap_meanref.RData')
