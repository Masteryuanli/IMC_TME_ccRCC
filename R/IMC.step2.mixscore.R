#!/usr/bin/Rscript
pars <- commandArgs(trailingOnly = TRUE)
if(length(pars)!=4){
  print("Description: caculate mixscore using mixing_score_summary function by SPIAT R package")
  print("Rscript IMC.step2.mixscore.R histocatfile rdatafile outdire prefix")
  quit()
}
library(ggplot2)
library(S4Vectors)

histocatfile<-pars[1]
rdatafile<-pars[2]
outdire<-pars[3]
prefix<-pars[4]

# histocatfile<-"example_data/P0_histocat_outfile.csv"
# rdatafile<-"1Cellannotate/sc.cellAnnotate.Rds"
# outdire<-"2Mixscore"
# prefix<-"P0"

sc.spat<-readRDS(rdatafile)
mypheno<-sc.spat@meta.data
exprs.sc<-as.matrix(sc.spat@assays$SCT@data)

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
sourceDir("/backup/niyl/Project/ZhangDW/IMC_Test/scripts/SPIAT-main/R/")

histocat.dt<-read.table(histocatfile,sep = ",",header = T,stringsAsFactors = F)
myspatialData<-DataFrame(histocat.dt[,c("Orientation","X_position","Y_position")])
rownames(myspatialData)<-paste(prefix,c(1:dim(histocat.dt)[1]),sep = "_")
mixscore.df<-data.frame(Phenotype=mypheno$Main.celltype,Cell.X.Position=rep("",nrow(mypheno)),Cell.Y.Position=rep("",nrow(mypheno)),Cell.Type=rep("",nrow(mypheno)),row.names =rownames(mypheno) )
mixscore.df$Cell.Type<-as.character(mixscore.df$Phenotype)
mixscore.df$Cell.Type[mixscore.df$Cell.Type%in%c("Unknown","Tumor Cell","Fibroblast","Endothelial Cell","Mesenchymal Cell")]<-"nonImmune"
mixscore.df$Cell.Type[mixscore.df$Cell.Type%in%c("B Cell","T Cell","Macrophage Cell")]<-"Immune"
myspatialData1<-myspatialData[rownames(mypheno),]
mixscore.df$Cell.X.Position<-myspatialData1$X_position
mixscore.df$Cell.Y.Position<-myspatialData1$Y_position
mixscore.final<-mixing_score_summary(mixscore.df,radius = 30,reference_celltype = "Immune",target_celltype = "nonImmune")
mixscore.final$patientID<-prefix
exprs.sub<-exprs.sc[,rownames(mypheno)]
mixscore.df.exprs<-cbind(mixscore.df,t(exprs.sub))
p1<-ggplot(mixscore.df,aes(x=Cell.X.Position,y=Cell.Y.Position,colour=Cell.Type ))+geom_point(size=0.5,alpha=0.6)+theme_bw()+
    labs(colour="CellType",x="Position of Cell X ",y="Position of Cell Y")+
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
          legend.text = element_text(size = 10,face = "bold"),
          legend.title = element_text(size = 10,face = "bold"),
          axis.text=element_text(size = 10,face = "bold"),axis.title = element_text(size = 10,face = "bold"))+
    scale_color_manual(values = c("#FF0000","#DDDDDD"),labels=c("B;T;Macrophage","Other"))+scale_y_reverse()
p1
ggsave(paste0(outdire,"/",prefix,".immume-nonimmune.png"),plot = p1,width = 6,height = 4)
mixscore.final$Type[mixscore.final$Mixing_score<=1.4]<-"Clustered"
write.table(mixscore.final,paste0(outdire,"/",prefix,".mixscore.csv"),col.names = T,row.names = F,sep = ",",quote = F)




