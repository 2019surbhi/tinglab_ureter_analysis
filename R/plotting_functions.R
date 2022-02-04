library(Seurat)

library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)

library(openxlsx)
library(dplyr)

##### Collapsed Heatmap Functions #####

### Get differential marker table function ###

## Extracts Differential markers from differential gene tables created by our Seurat pipeline. If your tables were generated differently, you may need to modify the code accordingly ##

# fname_list - vector of differential gene filenames (must include full path)
# n - number of genes to extract (per cluster) for printing on Heatmap
# clus_order - Specify cluster order on heatmap so that the genelist can be also sorted accordingly [By default set to 'default' denoting no sorting needed]

get_diff_marker_table<-function(fname_list,n,clus_order='default')
  
{
  
  files_list<-lapply(fname_list,read.table,header=TRUE)
  
  files_list<-lapply(files_list,head,n)
  
  name.split<-strsplit(fname_list,split='_')
  j<-grep('cluster',name.split[[1]])
  clus<-data.table::transpose(name.split)[[j]]
  names(files_list)<-clus
  
  genes<-lapply(X=1:length(files_list), function(x) { return(files_list[[x]]$gene)})
  
  names(genes)<-clus
  
  m<-vector()
  clus<-vector()
  
  for(i in 1:length(genes))
  {
    m<-append(m,unlist(genes[[i]]))
    clus<-append(clus,rep(names(genes[i]),n))
  }
  
  diff_tab<-data.frame("diff_markers"=m,"cluster"=clus)
  
  #Add numeric cluster to compare with marker ref
  clus<-gsub('cluster','',diff_tab$cluster) %>% as.numeric()
  diff_tab$cluster.num<-clus
  
  # Sort (optional)
  
  if(clus_order!='default')
  {o<-match(clus_order,diff_tab$cluster.num)
  o<-lapply(1:length(o), function(x){return(seq(o[x],(o[x]+(n-1))))}) %>% unlist()
  diff_tab<-diff_tab[o,]
  }
  return(diff_tab)
  
}


### Collapsed cluster matrix function ###

## This function extracts gene x cell matrix from input Seurat obj and returns the mean gene expression to collapse this matrx to gene x cluster matrix ##

# obj - seurat obj
# assay - Specify 'integrated','RNA' or other assay to extract gene exp from [Default=integrated]
# slot - data slot to use (options are 'counts', 'data' or 'scaled.data'. Default='data')

get_cluster_mat<-function(obj,assay='integrated',slot='data')
{
  
  DefaultAssay(obj)<-assay
  if(slot=='scaled.data')
  {
    if(dim(obj@assays[[assay]]@scale.data)==0)
    {
      obj<-ScaleData(obj)
    }
  }
  
  clus<-levels(obj)
  collapsed_mat<-vector()
  
  for(i in 1: length(clus))
  {
    sub<-NULL
    sub<-subset(obj,idents = clus[i])
    sub_mat<-GetAssayData(object = sub,slot = slot,assay=assay) %>% as.matrix()
    rm(sub)
    
    clus_sum<-rowSums(sub_mat)
    clus_avg<-clus_sum/ncol(sub_mat)
    rm(sub_mat)
    
    collapsed_mat<-cbind(collapsed_mat,clus_avg)
    rm(clus_sum)
    rm(clus_avg)
  }
  
  
  colnames(collapsed_mat)<-clus
  
  return(collapsed_mat)
  
}



### Print Heatmap function ### 

# obj - seurat obj
# assay - Specify 'integrated','RNA' or other assay to extract gene exp from [Default=integrated]
# slot - data slot to use (options are 'counts', 'data' or 'scaled.data'. Default='data')
# markers - differential markers sorted by cluster order
# umap_col - specify a vector of colors for cluster annotation
# out - output directory [Default='./']
# file_prefix - file prefix to add to output file name
# clus_order- specify vector of integers or names to denote cluster order (Set to 'default' which will take current cluster order in Seurat obj)
# cr - set to TRUE  to cluster heatmap rows, FALSE otherwise
# cc -set to TRUE  to cluster heatmapcolumns, FALSE otherwise
# w= - Heatmap width in inches [Default=11]
# h - Heatmap height in inches [Default=11]

print_mean_exp_per_cluster_heatmap<-function(obj,assay='integrated',slot='data',markers,umap_col,out='./',file_prefix,clus_order='default',cr=FALSE,cc=FALSE,w=11,h=11)
{
  
  # Get per cluster mean expression data
  clus.mat<-get_cluster_mat(obj,assay,slot)
  
  if(clus_order=='default')
  {
    clus_order<-levels(obj)
  }
  
  # Subset genes
  rows<-match(markers,rownames(clus.mat))
  sub.mat<- clus.mat[rows,]
  
  
  #Order columns
  o<-match(clus_order,colnames(sub.mat))
  sub.mat<-sub.mat[,o]
  
  # Heatmap annotations
  cell_type<-colnames(sub.mat)
  
  #Add name to umap_col
  names(umap_col)<-cell_type
  
  top<-HeatmapAnnotation(cell_types=cell_type,
                         col=list(cell_types=umap_col),
                         show_legend = TRUE)
  
  # Heatmap colors
  col<-c("#2C7BB6","#FFFFBF","#D7191C")
  hmcols<-colorRampPalette(col)(256)
  
  # Heatmap name
  name<-paste0('mean_expression (',slot,' slot)')
  fname<-paste0(out,file_prefix,name,'.png')
  
  # Print heatmap
  
  png(fname,width=w,height=h,units='in',res=300)
  draw(Heatmap(sub.mat, col= hmcols,
               show_heatmap_legend = TRUE,
               name = name ,
               cluster_rows = cr,
               cluster_columns = cc,
               top_annotation = top,
  ))
  
  dev.off()
  
}


