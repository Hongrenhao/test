library(data.table)
filedir="F:/BaiduNetdiskDownload/NBT_CD34.counts.HS_BM_3.csv"
NBT=fread(input=filedir,sep=",",header=TRUE,data.table=FALSE)
row.names(NBT)=NBT[,1]
NBT=NBT[,-1]    #免疫细胞
NBT_mtr=as.matrix(NBT)
NBT_sp=as(NBT_mtr,"dgCMatrix")
library(dplyr)
library(Seurat)
NBT111=t(NBT_sp)
NBT_seurat<-CreateSeuratObject(counts =NBT111, project = "NBT",min.cells=290,min.features = 4110)
NBT_seurat
write.csv(NBT_seurat@assays$RNA@counts,file="E:/myresource/NBT/data/NBT_data_test.csv",quote=F,row.names=T,col.names=T)
