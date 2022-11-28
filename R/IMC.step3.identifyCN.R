#!/usr/bin/Rscript
pars <- commandArgs(trailingOnly = TRUE)
if(length(pars)!=6){
  print("Description: identify CellNeighborhood by SPIAT  R package")
  print("Rscript IMC.step3.identifyCN.R histocatfile rdatafile mixscorefile panelfile outdire prefix")
  quit()
}
histocatfile<-pars[1]
rdatafile<-pars[2]
mixscorefile<-pars[3]
panelfile<-pars[4]
outdire<-pars[5]
prefix<-pars[6]


library(ggplot2)
library(SingleCellExperiment)
library(SpatialExperiment)
library(magrittr)
library(dplyr)



# histocatfile<-"example_data/P0_histocat_outfile.csv"
# rdatafile<-"1Cellannotate/P0.sc.cellAnnotate.Rds"
# mixscorefile<-"2Mixscore/P0.mixscore.csv"
# panelfile<-"example_data/panel.csv"
# outdire<-"3CellNeighborhood/"
# prefix<-"P0"


sc.spat<-readRDS(rdatafile)
mixscore.file<-read.table(mixscorefile,header = T,sep =",")
mypanel<-read.table(panelfile,sep = ",",header = T)
#--load source code of SPIAT R package---#
sourceDir <- function(path, trace = TRUE, ...) {
  op <- options(); on.exit(options(op)) # to reset after each 
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
    options(op)
  }
}
sourceDir("R/SPIAT-main/R/")
splitlabel<-function(x){
  unlist(strsplit(x,split = "_"))[c(4)]
}

mypheno<-sc.spat@meta.data
type.sample<-mixscore.file$Type[mixscore.file$patientID%in%prefix]
histocat.dt<-read.table(histocatfile,sep = ",",header = T,stringsAsFactors = F)
myspatialData<-DataFrame(histocat.dt[,c("Orientation","X_position","Y_position")])
rownames(myspatialData)<-paste(prefix,c(1:dim(histocat.dt)[1]),sep = "_")
mycounts<-as.matrix(t(histocat.dt[,3:(2+17)]))
count.rowname<-sapply(rownames(mycounts),splitlabel)
names(count.rowname)<-NULL
rownames(mycounts)<-count.rowname
colnames(mycounts)<-paste(prefix,c(1:dim(histocat.dt)[1]),sep = "_")
myrowData<-DataFrame(mypanel[c(3:19),c(3,4,5)])
mycolData<-data.frame(prefix=rep(prefix,dim(histocat.dt)[1]),histocat.dt[,c("X_position","Y_position")],
                        row.names =paste(prefix,c(1:dim(histocat.dt)[1]),sep = "_"),Cell.Type="")
colnames(mycolData)<-c("prefix","Cell.X.Position","Cell.Y.Position","Cell.Type")
mycolData$Cell.Type<-mypheno$minor.celltype[match(rownames(mycolData),rownames(mypheno),0L)]
metadata<-DataFrame(mypheno)
mycounts<-as(as.matrix(mycounts), "dgCMatrix")
#---creat spatial object---#
spe.obj<- SpatialExperiment(
    assay = mycounts, 
    colData = mycolData, 
    spatialCoords  = myspatialData)
pheno.spe<-colData(spe.obj)
allcelltype<- sort(unique(pheno.spe$Cell.Type))
if(type.sample=="Clustered"){
    target.cell<-allcelltype[!grepl("Endothelial|Mesenchymal|Fibroblast|Tumor|Unknown",allcelltype)]
  }else if(type.sample=="Scattered"){
    target.cell<-allcelltype[!grepl("Unknown",allcelltype)]
  }
CN.dt <- identify_neighborhoods.new(spe.obj, 
                                       method = "hierarchical",
                                       min_neighborhood_size = 10,
                                       cell_types_of_interest = target.cell, 
                                       radius = 10, feature_colname = "Cell.Type")
clusters.dt<-CN.dt$sc
coldata.isample<-clusters.dt@colData
coldata.isample.cluster<-coldata.isample[grep("Cluster",coldata.isample$Neighborhood),]
coldata.isample.cluster.df<-as.data.frame(coldata.isample.cluster)
neighorhoods_vis <- composition_of_neighborhoods(clusters.dt, feature_colname = "Cell.Type")
neighorhoods_vis <- neighorhoods_vis[neighorhoods_vis$Total_number_of_cells >=5,]
neighorhoods_vis1<-neighorhoods_vis[neighorhoods_vis$Cell.Type%in%target.cell,]
 #select cluster
neighorhoods_vis.out<-neighorhoods_vis1[grep("Cluster",neighorhoods_vis1$Neighborhood),]
#just have one Neighborhood
if(length(unique(neighorhoods_vis.out$Neighborhood))==1){
    cn.enrich<-paste0(paste(neighorhoods_vis.out$Cell.Type[neighorhoods_vis.out$Number_of_cells>5],collapse = ";")," Enriched")
    cn.annotate<-data.frame(Enriched=cn.enrich,Cluster="CN_1")
    write.table(cn.annotate,paste0(outdire,"/",prefix,"_CN-Enriched-annotate.tsv"),
                sep = "\t",col.names = T,row.names = F,quote = F)
      
      next
  }else{
     cn.annotate<-composition_cn_enrich_annotate(composition = neighorhoods_vis1,feature_colname = "Cell.Type",number = 30)
     cn.annotate<-cn.annotate[cn.annotate$Enrich!="",]
   #--------------------------plot--------------------#
    cells_in_clusters <-  as.data.frame(coldata.isample[grepl("Cluster",coldata.isample$Neighborhood),])
    cells_not_in_clusters <- as.data.frame(coldata.isample[!grepl("Cluster",coldata.isample$Neighborhood),])
    colnames(cells_in_clusters)[5]<-"Cluster"
    colnames(cells_not_in_clusters)[5]<-"Cluster"
    cluster.newannotate.df<-data.frame(ClusterEnriched=unique(cn.annotate$Enrich),cluster.new.num= rep(c(1:length(unique(cn.annotate$Enrich))))  )
    cells_in_clusters$ClusterEnriched<-factor(cells_in_clusters$Cluster,levels =cn.annotate$Cluster,labels = cn.annotate$Enrich )
    cells_in_clusters<-cells_in_clusters[is.na(cells_in_clusters$ClusterEnriched)==F,]
    cells_in_clusters$num.cluster<-factor(cells_in_clusters$ClusterEnriched,levels =cluster.newannotate.df$ClusterEnriched,labels = cluster.newannotate.df$cluster.new.num )
    number_of_clusters <- length(unique(cells_in_clusters$ClusterEnriched))
      
   
    #---------------no annotate CN---------------#
    if( length(which(cells_in_clusters$ClusterEnriched==""))!=0){
      cells_in_clusters.noriched<-cells_in_clusters[cells_in_clusters$ClusterEnriched=="",]
      cells_not_in_clusters1<-rbind(cells_not_in_clusters,cells_in_clusters.norich[,c(1:5)])
      
    }else{
      cells_not_in_clusters1<-cells_not_in_clusters
    }
    
    cells_in_clusters1<-cells_in_clusters[cells_in_clusters$ClusterEnriched!="",]
    
    #label the Cluster centre by averaging x and y
    cells_in_clusters1$ClusterEnriched<-as.character(cells_in_clusters1$ClusterEnriched)
    cells_in_clusters1$ClusterEnriched<-factor(cells_in_clusters1$ClusterEnriched,levels =cluster.newannotate.df$ClusterEnriched )
    cells_in_clusters1$num.cluster<-factor(cells_in_clusters1$num.cluster,levels =cluster.newannotate.df$cluster.new.num )
    cells_in_clusters<-cells_in_clusters1
    label_location <- vector()
    for (Clusternumber in seq_len(number_of_clusters)) {
        cells_in_Cluster <- cells_in_clusters[cells_in_clusters$num.cluster == Clusternumber, ]
        minX <- min(cells_in_Cluster$Cell.X.Position)
        maxX <- max(cells_in_Cluster$Cell.X.Position)
        minY <- min(cells_in_Cluster$Cell.Y.Position)
        maxY <- max(cells_in_Cluster$Cell.Y.Position)
        averageX <- (minX + maxX)/2
        averageY <- (minY + maxY)/2
        label_location <- rbind(label_location,c(Clusternumber, averageX, averageY))
      }
      label_location <- as.data.frame(label_location)
      colnames(label_location) <- c("Cluster", "Xpos", "Ypos")
      # use colourblind-friendly colours
      
      allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                  "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
                  "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
                  "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
      
      cluster_colours <- allcolour[seq_len(number_of_clusters)]
      cells_in_clusters$ClusterEnriched<-paste0(cells_in_clusters$ClusterEnriched," Enriched")
      q <- ggplot(cells_in_clusters, aes(x=Cell.X.Position, y=Cell.Y.Position))
      q <- q + geom_point(aes(color = ClusterEnriched,size=ClusterEnriched), size = 1)+
        guides(color = guide_legend(title = "Enriched"))
      q <- q + scale_color_manual(values=cluster_colours)
      if(dim(cells_not_in_clusters1)[1]!=0){
        q <- q + geom_point(data = cells_not_in_clusters1,  colour = "black",size = 0.5)
      }
      q <- q + xlab("Cell.X.Position") + ylab("Cell.Y.Position") +
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              legend.position = "none",panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
              axis.title = element_text( size=12,face="bold"),
              axis.text=  element_text( size=8,face="bold",hjust = 0.9))+scale_y_reverse()
      q <- q + xlab("Cell.X.Position") + ylab("Cell.Y.Position") +
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
              axis.title = element_text( size=12,face="bold"),
              axis.text=  element_text( size=8,face="bold",hjust = 0.9))+scale_y_reverse()
    ggsave(paste0(outdire,"/",prefix,"_cell_neighborhood-label-legend.png"),plot = q,width = 8,height = 6)
    
    cluster.newannotate.df$cluster.new.num<-paste("CN",cluster.newannotate.df$cluster.new.num,sep = "_")
    cluster.newannotate.df$ClusterEnriched<-paste(cluster.newannotate.df$ClusterEnriched,"Enriched",sep = " ")
    colnames(cluster.newannotate.df)<-c("Enriched","CN")
    write.table(cluster.newannotate.df,paste0(outdire,"/",prefix,"_CN-Enriched-annotate.tsv"),
                  sep = "\t",col.names = T,row.names = F,quote = F)
    }
  
  
 




