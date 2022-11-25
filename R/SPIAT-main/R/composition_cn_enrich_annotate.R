#' composition_cn_enrich_annotate
#' @description annotate according to the main cell type of each cluster
#' @param composition data.frame object in the form of the output of composition_of_neighborhoods function
#' @param feature_colname feature_colname String. Column with cell types.
#' @param number number for percentage of cell types
#' @import reshape2
#' @return An data.frame 


composition_cn_enrich_annotate<-function(composition,feature_colname,number=30){
  
  cluster_size <- unique(data.frame(Neighborhood = composition$Neighborhood,
                                    Total_cells = composition$Total_number_of_cells))
  rownames(cluster_size) <- cluster_size$Neighborhood
  cluster_size$Neighborhood <- NULL
  
  composition1 <- composition[,c(feature_colname, "Neighborhood", "Percentage")]
  composition2 <- reshape2::dcast(composition1, paste(feature_colname, "~", "Neighborhood"), value.var="Percentage")
  
  rownames(composition2) <- composition2[,feature_colname]
  composition2[,feature_colname] <- NULL
  composition2[is.na(composition2)] <- -1
  composition2 <- as.matrix(composition2)
  free.loc<-grep("Free_cell",colnames(composition2))
  composition2<-composition2[,-free.loc]
  composition3<-t(composition2)
  if( is.null(nrow(composition3))==T ){
    next
  }
  cn.annotate<-as.character()
  for(inow in c(1:nrow(composition3)) ){
    icluster<-sort(composition3[inow,],decreasing = T)
    icluster.celltypEnrich<-sort(names(icluster)[which(icluster>=number)])
    celltype.order<-c("B Cell","T Cell","T Cell-CD4","T Cell-CD8","Macrophage Cell",
                      "Macrophage Cell-CD163","Tumor Cell","Endothelial Cell","Mesenchymal Cell","Fibroblast")
    icluster.celltypEnrichs<-intersect(celltype.order,icluster.celltypEnrich)
    mat<-matrix(c(rownames(composition3)[inow],paste(icluster.celltypEnrichs,collapse = ";") ),1)
    cn.annotate<-rbind(cn.annotate,mat)
  }
  colnames(cn.annotate)<-c("Cluster","Enrich")
  cn.annotate<-as.data.frame(cn.annotate)
  return(cn.annotate)
}
