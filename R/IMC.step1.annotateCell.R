#!/usr/bin/Rscript
pars <- commandArgs(trailingOnly = TRUE)
if(length(pars)!=4){
  print("Description: cluster dimension reduction and cell identity was defined by seurat R package")
  print("Rscript IMC.step1.annotateCell.R histocatfile panelfile outdire prefix")
  quit()
}

library(Seurat)
library(RColorBrewer)
library(ggbreak)
library(ggthemes)
library(ggplot2)
library(viridis)


histocatfile<-pars[1]
panelfile<-pars[2]
outdire<-pars[3]
prefix<-pars[4]

# histocatfile<-"example_data/P0_histocat_outfile.csv"
# panelfile<-"example_data/panel.csv"
# outdire<-"1Cellannotate/"
# prefix<-"P0"

mypanel<-read.table(panelfile,sep = ",",header = T)
Clean_Target<-mypanel$Clean_Target
MetalTag<-mypanel$Metal.Tag

splitlabel<-function(x){
  unlist(strsplit(x,split = "_"))[c(4)]
}

histocatfile<-read.table(histocatfile,sep = ",",header = T,stringsAsFactors = F)
mycounts<-as.matrix(t(histocatfile[,3:(length(Clean_Target))]))
count.rowname<-sapply(rownames(mycounts),splitlabel)
count.rowname[grep("PD_L1",names(count.rowname))]<-"PDL1"
count.rowname[grep("PD_1",names(count.rowname))]<-"PD1"
names(count.rowname)<-NULL
rownames(mycounts)<-count.rowname
colnames(mycounts)<-paste(prefix,c(1:dim(mycounts)[2]),sep = "_")
rownames(mycounts)<-gsub("_","",rownames(mycounts))
seurat.obj = CreateSeuratObject(counts = mycounts, assay="Spatial")
coord.df = data.frame(x=histocatfile$X_position, y=histocatfile$Y_position, stringsAsFactors=FALSE)
rownames(coord.df) = paste(prefix,c(1:dim(histocatfile)[1]),sep = "_")

seurat.obj@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = coord.df
  )


seurat.obj.st <- SCTransform(seurat.obj, assay = "Spatial", verbose = FALSE,ncells =5000)
merge.sc.pca <- RunPCA(seurat.obj.st, npcs = 30, features = rownames(seurat.obj.st))
merge.sc.nb <- FindNeighbors(merge.sc.pca, dims = 1:6)
merge.sc.cluster <- FindClusters(merge.sc.nb,resolution = 0.6, verbose = FALSE)
merge.sc <- RunUMAP(merge.sc.cluster, dims = 1:6)

#-----conventional drawing---------------#
Cluster_10_cols<-c(brewer.pal(9,"YlOrRd")[2:5],brewer.pal(9,"YlGnBu")[2:5],brewer.pal(9,"YlGn")[2:5],brewer.pal(9,"Reds")[2:5]
                   ,brewer.pal(9,"RdPu")[2:5],brewer.pal(9,"Purples")[2:5],brewer.pal(9,"PuRd")[2:5],brewer.pal(9,"Greens")[2:5],
                   brewer.pal(9,"GnBu")[2:5],brewer.pal(9,"BuPu")[2:5])

p1<-DimPlot(merge.sc,cols =Cluster_10_cols,label = T)+theme(axis.title= element_text( size=10,face="bold"),
                                                            axis.text  =  element_text( size=10,face="bold"))
ggsave(paste0(outdire,"/sc_cluster_umap.png"),plot = p1,width = 4,height = 3)


library(viridis)
DefaultAssay(merge.sc)<-"SCT"
#membrane gene
f.cd3<-FeaturePlot(merge.sc,features = "CD3")+scale_color_viridis()
f.cd20<-FeaturePlot(merge.sc,features = "CD20")+scale_color_viridis()
f.cd68<-FeaturePlot(merge.sc,features = "CD68")+scale_color_viridis()
f.cd4<-FeaturePlot(merge.sc,features = "CD4")+scale_color_viridis()
f.cd8<-FeaturePlot(merge.sc,features = "CD8a")+scale_color_viridis()
f.cd31<-FeaturePlot(merge.sc,features = "CD31")+scale_color_viridis()
f.panck<-FeaturePlot(merge.sc,features = "PanCK")+scale_color_viridis()
f.ecad<-FeaturePlot(merge.sc,features = "Ecad")+scale_color_viridis()
f.asma<-FeaturePlot(merge.sc,features = "aSMA")+scale_color_viridis()
f.vim<-FeaturePlot(merge.sc,features = "Vim")+scale_color_viridis()



ggsave(paste0(outdire,"/CD3_featureplot.png"),plot = f.cd3,width = 3,height = 2)
ggsave(paste0(outdire,"/CD20_featureplot.png"),plot = f.cd20,width = 3,height = 2)
ggsave(paste0(outdire,"/CD68_featureplot.png"),plot = f.cd68,width = 3,height = 2)
ggsave(paste0(outdire,"/CD4_featureplot.png"),plot = f.cd4,width = 3,height = 2)
ggsave(paste0(outdire,"/CD8_featureplot.png"),plot = f.cd8,width = 3,height = 2)
ggsave(paste0(outdire,"/CD31_featureplot.png"),plot = f.cd31,width = 3,height = 2)
ggsave(paste0(outdire,"/PanCK_featureplot.png"),plot = f.panck,width = 3,height = 2)
ggsave(paste0(outdire,"/ECad_featureplot.png"),plot = f.ecad,width = 3,height = 2)
ggsave(paste0(outdire,"/aSMA_featureplot.png"),plot = f.asma,width = 3,height = 2)
ggsave(paste0(outdire,"/Vim_featureplot.png"),plot = f.vim,width = 3,height = 2)

# merge.sc<-readRDS("1Cellannotate/sc.Rds")
#---function gene----#
function.gene<-c("PDL1","PD1","FOXP3","PD1","PDL1","GZMB","KI67","CD45RO")


for(imarker in function.gene){
  f.plot<-FeaturePlot(merge.sc,features = imarker)+scale_color_viridis()
  ggsave(paste0(outdire,"/",imarker,"_featureplot.png"),plot = f.plot,width = 3,height = 2)
}


exprs.sc<-as.matrix(merge.sc@assays$SCT@data)
#-------density positive marker-----------------
exprs.df<-as.data.frame(t(exprs.sc))
cd20.df<-as.data.frame(table(exprs.df$CD20))
colnames(cd20.df)<-c("exprs","freq")
cd20.df$exprs<-as.numeric(as.character(cd20.df$exprs))
p.cd20<-ggplot(cd20.df,aes(x=exprs,y=freq))+geom_line()+theme_bw()+geom_vline(xintercept =0.69,linetype='dashed')+
  theme_stata()+labs(x="",y="")+
  labs(x="Normalized mean intensity",y="cell numbers",title="CD20")+
  scale_x_continuous(breaks = c(0,0.69,1,2,3,4))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.3))
# p.cd20
ggsave(paste0(outdire,"/CD20.densityplot.png"),plot = p.cd20,width = 3,height = 4)

cd68.df<-as.data.frame(table(exprs.df$CD68))
colnames(cd68.df)<-c("exprs","freq")
cd68.df$exprs<-as.numeric(as.character(cd68.df$exprs))
P.CD68<-ggplot(cd68.df,aes(x=exprs,y=freq))+geom_line()+theme_bw()+geom_vline(xintercept =1.6 ,linetype='dashed')+
  theme_stata()+labs(x="Normalized mean intensity",y="cell numbers",title="CD68")+
  scale_x_continuous(breaks = c(0,1,1.6,2,3,4))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.6))
ggsave(paste0(outdire,"/CD68.densityplot.png"),plot = P.CD68,width = 3,height = 4)


aSMA.df<-as.data.frame(table(exprs.df$aSMA))
colnames(aSMA.df)<-c("exprs","freq")
aSMA.df$exprs<-as.numeric(as.character(aSMA.df$exprs))
p.asam<-ggplot(aSMA.df,aes(x=exprs,y=freq))+geom_line()+theme_bw()+geom_vline(xintercept =1.6,linetype='dashed')+
  theme_stata()+scale_x_continuous(breaks = c(0,1,1.4,2,3,4))+
  labs(x="Normalized mean intensity",y="cell numbers",title="aSMA")+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.6))
ggsave(paste0(outdire,"/aSAM.densityplot.png"),plot = p.asam,width = 3,height = 4)

Vim.df<-as.data.frame(table(exprs.df$Vim))
colnames(Vim.df)<-c("exprs","freq")
Vim.df$exprs<-as.numeric(as.character(Vim.df$exprs))
p.vim<-ggplot(Vim.df,aes(x=exprs,y=freq))+geom_line()+theme_bw()+geom_vline(xintercept =2,linetype='dashed')+
  theme_stata()+labs(x="Normalized mean intensity",y="cell numbers",title="Vimentin")+
  scale_x_continuous(breaks = c(0,1,2,3,4))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.6))
ggsave(paste0(outdire,"/VIM.densityplot.png"),plot = p.vim,width = 3,height = 4)

CD163.df<-as.data.frame(table(exprs.df$CD163))
colnames(CD163.df)<-c("exprs","freq")
CD163.df$exprs<-as.numeric(as.character(CD163.df$exprs))
p.CD163<-ggplot(CD163.df,aes(x=exprs,y=freq))+geom_line()+geom_vline(xintercept =1,linetype='dashed')+
  theme_stata()+
  labs(x="Normalized mean intensity",y="cell numbers",title="CD163")+
  scale_x_continuous(breaks = c(0,1,2,3,4))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.6))
ggsave(paste0(outdire,"/CD163.densityplot.png"),plot = p.CD163,width = 3,height = 4)


PanCK.df<-as.data.frame(table(exprs.df$PanCK))
colnames(PanCK.df)<-c("exprs","freq")
PanCK.df$exprs<-as.numeric(as.character(PanCK.df$exprs))
p.panck<-ggplot(PanCK.df,aes(x=exprs,y=freq))+geom_line()+geom_vline(xintercept =1,linetype='dashed')+
  theme_stata()+
  labs(x="Normalized mean intensity",y="cell numbers",title="Pan-keratin")+
  scale_x_continuous(breaks = c(0,0.69,1,2,3))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.3))
ggsave(paste0(outdire,"/panck.densityplot.png"),plot = p.panck,width =4,height = 4)


PDL1.df<-as.data.frame(table(exprs.df$PDL1))
colnames(PDL1.df)<-c("exprs","freq")
PDL1.df$exprs<-as.numeric(as.character(PDL1.df$exprs))
p.pdl1<-ggplot(PDL1.df,aes(x=exprs,y=freq))+geom_line()+geom_vline(xintercept =0.69,linetype='dashed')+
  theme_stata()+
  labs(x="Normalized mean intensity",y="cell numbers",title="PD-L1")+
  scale_x_continuous(breaks = c(0,0.69,1,2,3,4))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.3))
ggsave(paste0(outdire,"/pdl1.densityplot.png"),plot = p.pdl1,width = 3,height = 4)


CD31.df<-as.data.frame(table(exprs.df$CD31))
colnames(CD31.df)<-c("exprs","freq")
CD31.df$exprs<-as.numeric(as.character(CD31.df$exprs))
p.cd31<-ggplot(CD31.df,aes(x=exprs,y=freq))+geom_line()+geom_vline(xintercept =1,linetype='dashed')+
  theme_stata()+
  labs(x="Normalized mean intensity",y="cell numbers",title="CD31")+
  scale_x_continuous(breaks = c(0,0.69,1,2,3,4))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.6))
ggsave(paste0(outdire,"/cd31.densityplot.png"),plot = p.cd31,width = 3.5,height = 4)



FOXP3.df<-as.data.frame(table(exprs.df$FOXP3))
colnames(FOXP3.df)<-c("exprs","freq")
FOXP3.df$exprs<-as.numeric(as.character(FOXP3.df$exprs))
p.foxp3<-ggplot(FOXP3.df,aes(x=exprs,y=freq))+geom_line()+geom_vline(xintercept =0.69,linetype='dashed')+
  theme_stata()+
  labs(x="Normalized mean intensity",y="cell numbers",title="Foxp3")+
  scale_x_continuous(breaks = c(0,0.69,1,2,3,4))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.3))
# p.foxp3
ggsave(paste0(outdire,"/foxp3.densityplot.png"),plot = p.foxp3,width = 3,height = 4)


CD4.df<-as.data.frame(table(exprs.df$CD4))
colnames(CD4.df)<-c("exprs","freq")
CD4.df$exprs<-as.numeric(as.character(CD4.df$exprs))
p.cd4<-ggplot(CD4.df,aes(x=exprs,y=freq))+geom_line()+geom_vline(xintercept =2,linetype='dashed')+
  theme_stata()+
  labs(x="Normalized mean intensity",y="cell numbers",title="CD4")+
  scale_x_continuous(breaks = c(0,1,2,3,4))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.6))

ggsave(paste0(outdire,"/cd4.densityplot.png"),plot = p.cd4,width = 3,height = 4)


Ecad.df<-as.data.frame(table(exprs.df$Ecad))
colnames(Ecad.df)<-c("exprs","freq")
Ecad.df$exprs<-as.numeric(as.character(Ecad.df$exprs))
p.ecad<-ggplot(Ecad.df,aes(x=exprs,y=freq))+geom_line()+geom_vline(xintercept =0.69,linetype='dashed')+
  theme_stata()+
  labs(x="Normalized mean intensity",y="cell numbers",title = "Ecad")+
  scale_x_continuous(breaks = c(0,0.69,1,2,3,4))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.3))
p.ecad
ggsave(paste0(outdire,"/ecad.densityplot.png"),plot = p.ecad,width = 3,height = 4)

CD8a.df<-as.data.frame(table(exprs.df$CD8a))
colnames(CD8a.df)<-c("exprs","freq")
CD8a.df$exprs<-as.numeric(as.character(CD8a.df$exprs))
p.cd8<-ggplot(CD8a.df,aes(x=exprs,y=freq))+geom_line()+geom_vline(xintercept =1,linetype='dashed')+
  theme_stata()+
  labs(x="Normalized mean intensity",y="cell numbers",title="CD8a")+
  scale_x_continuous(breaks = c(0,0.69,1,2,3,4))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.3))
ggsave(paste0(outdire,"/CD8A.densityplot.png"),plot = p.cd8,width = 3,height = 4)

PD1.df<-as.data.frame(table(exprs.df$PD1))
colnames(PD1.df)<-c("exprs","freq")
PD1.df$exprs<-as.numeric(as.character(PD1.df$exprs))
p.PD1<-ggplot(PD1.df,aes(x=exprs,y=freq))+geom_line()+geom_vline(xintercept =0.69,linetype='dashed')+
  theme_stata()+
  labs(x="Normalized mean intensity",y="cell numbers",title="PD1")+
  scale_x_continuous(breaks = c(0,0.69,1,2,3,4))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.3))
p.PD1
ggsave(paste0(outdire,"/PD1.densityplot.png"),plot = p.PD1,width = 3,height = 4)


GZMB.df<-as.data.frame(table(exprs.df$GZMB))
colnames(GZMB.df)<-c("exprs","freq")
GZMB.df$exprs<-as.numeric(as.character(GZMB.df$exprs))
p.GZMB<-ggplot(GZMB.df,aes(x=exprs,y=freq))+geom_line()+geom_vline(xintercept =1,linetype='dashed')+
  theme_stata()+ labs(x="Normalized mean intensity",y="cell numbers",title="GZMB")+
  scale_x_continuous(breaks = c(0,1,2,3,4))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.6))
ggsave(paste0(outdire,"/GZMB.densityplot.png"),plot = p.GZMB,width = 3,height = 4)

KI67.df<-as.data.frame(table(exprs.df$KI67))
colnames(KI67.df)<-c("exprs","freq")
KI67.df$exprs<-as.numeric(as.character(KI67.df$exprs))
p.ki67<-ggplot(KI67.df,aes(x=exprs,y=freq))+geom_line()+geom_vline(xintercept =0.69,linetype='dashed')+
  theme_stata()+
  labs(x="Normalized mean intensity",y="cell numbers",title = "Ki67")+
  scale_x_continuous(breaks = c(0,0.69,1,2,3,4))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.3))
p.ki67
ggsave(paste0(outdire,"/KI67.densityplot.png"),plot = p.ki67,width = 3,height = 4)

CD3.df<-as.data.frame(table(exprs.df$CD3))
colnames(CD3.df)<-c("exprs","freq")
CD3.df$exprs<-as.numeric(as.character(CD3.df$exprs))
p.cd3<-ggplot(CD3.df,aes(x=exprs,y=freq))+geom_line()+geom_vline(xintercept =0.69,linetype='dashed')+
  theme_stata()+
  labs(x="Normalized mean intensity",y="cell numbers",title="CD3")+
  scale_x_continuous(breaks = c(0,0.69,1,2,3,4))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.3))

ggsave(paste0(outdire,"/CD3.densityplot.png"),plot = p.cd3,width = 3,height = 4)



CD45RO.df<-as.data.frame(table(exprs.df$CD45RO))
colnames(CD45RO.df)<-c("exprs","freq")
CD45RO.df$exprs<-as.numeric(as.character(CD45RO.df$exprs))
p.cd45ro<-ggplot(CD45RO.df,aes(x=exprs,y=freq))+geom_line()+geom_vline(xintercept =0.69,linetype='dashed')+
  theme_stata()+labs(x="Normalized mean intensity",y="cell numbers",title="CD45RO")+scale_x_continuous(breaks = c(0,0.69,1,2,3,4))+
  theme(axis.text.y  = element_text(size = 6,hjust =0.5),axis.text.x  = element_text(size =6,hjust =0.3))

ggsave(paste0(outdire,"/CD45RO.densityplot.png"),plot = p.cd45ro,width = 3,height = 4)

exprs.sc<-merge.sc@assays$SCT@data
barcode.b<-colnames(exprs.sc[,which(exprs.sc["CD20",]>=0.69)])
barcode.t<-colnames(exprs.sc[,which(exprs.sc["CD3",]>=0.69)])
barcode.cd4<-colnames(exprs.sc[,which(exprs.sc["CD4",]>=2)])
barcode.cd8<-colnames(exprs.sc[,which(exprs.sc["CD8a",]>=1)])
barcode.mar<-colnames(exprs.sc[,which(exprs.sc["CD68",]>=1.6)])
barcode.cd163<-colnames(exprs.sc[,which(exprs.sc["CD163",]>=1)])
barcode.cd31<-colnames(exprs.sc[,which(exprs.sc["CD31",]>=1)])
barcode.vim<-colnames(exprs.sc[,which(exprs.sc["Vim",]>=2)])
barcode.Ecad<-colnames(exprs.sc[,which(exprs.sc["Ecad",]>=0.69)])
barcode.aSMA<-colnames(exprs.sc[,which(exprs.sc["aSMA",]>=1.6)])
barcode.PanCK<-colnames(exprs.sc[,which(exprs.sc["PanCK",]>=1)])


merge.sc.mar<-SetIdent(merge.sc, cells = barcode.mar, value = 'Macrophage Cell')
merge.sc.t<-SetIdent(merge.sc.mar, cells = barcode.t, value = 'T Cell')

merge.sc.B<-SetIdent(merge.sc.t, cells = barcode.b, value = 'B Cell')
alltype.stat<-table(Idents(merge.sc.B))
barcode.nonimmune<-WhichCells(merge.sc.B,idents = names(alltype.stat)[4:length(alltype.stat)] )

#aSMA  Fibroblast
barcode.asma.keep<-intersect(barcode.nonimmune,barcode.aSMA)
merge.sc.fibr<-SetIdent(merge.sc.B, cells = barcode.asma.keep, value = 'Fibroblast')

#Endotheilal
barcode.cd31.keep<-intersect(barcode.nonimmune,barcode.cd31)
merge.sc.Endo<-SetIdent(merge.sc.fibr, cells = barcode.cd31.keep, value = 'Endothelial Cell')


#Vim Mesenchymal cells
barcode.vim.keep<-intersect(barcode.nonimmune,barcode.vim)
merge.sc.messchymal<-SetIdent(merge.sc.Endo, cells = barcode.vim.keep, value = 'Mesenchymal Cell')

#Ecad panck tumor cell
barcode.Ecad.keep<-intersect(barcode.nonimmune,barcode.Ecad)
barcode.panck.keep<-intersect(barcode.nonimmune,barcode.PanCK)
tumor.barcode<-c(barcode.Ecad,barcode.PanCK)
merge.sc.tumor<-SetIdent(merge.sc.messchymal, cells = tumor.barcode, value = 'Tumor Cell')

#unknown cell
idents.num<-names(table(Idents(merge.sc.tumor))[8:length(table(Idents(merge.sc.tumor)))])
barcode.unknown<-WhichCells(merge.sc.tumor,idents =idents.num)
merge.sc.annotate<-SetIdent(merge.sc.tumor, cells = barcode.unknown, value = 'Unknown')

#-----identity coventational immune cell subtype---#
celltype.df<-as.data.frame(Idents(merge.sc.annotate))
colnames(celltype.df)<-"Main.celltype"
celltype.df$minor.celltype<-as.character(celltype.df$Main.celltype)


#CD4 T cell
barcode.T<-WhichCells(merge.sc.annotate,idents = "T Cell")
CD4.T.barcode<-intersect(barcode.T,barcode.cd4)
celltype.df[CD4.T.barcode,"minor.celltype"]<-"T Cell-CD4"

#CD8 T cell
barcode.T<-WhichCells(merge.sc.annotate,idents = "T Cell")
CD8.T.barcode<-intersect(barcode.T,barcode.cd8)
celltype.df[CD8.T.barcode,"minor.celltype"]<-"T Cell-CD8"

#CD163 macrophage cell
barcode.macro<-WhichCells(merge.sc.annotate,idents = "Macrophage Cell")
CD163.Macro.barcode<-intersect(barcode.macro,barcode.cd163)
celltype.df[CD163.Macro.barcode,"minor.celltype"]<-"Macrophage Cell-CD163"

merge.sc.addmeta<-AddMetaData(merge.sc.annotate,metadata =celltype.df )
levels(merge.sc.addmeta)<-rev( names(table(Idents(merge.sc.addmeta))) )
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#1E90FF","#DDDDDD")
DimPlot(merge.sc.addmeta,cols = allcolour)

saveRDS(merge.sc.addmeta,paste0(outdire,"/",prefix,".sc.cellAnnotate.Rds"))
