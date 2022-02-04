## Running Plotting functions

**Load plotting functions**

`source('./plotting_functions.R')`

**1. Printing collapsed heatmap - using our stromal subset as example**

**Specify cluster order** (We grouped/ordered these clusters in our pblication based on gene expression similarities/dissimilarities)

`clus_order_stromal<-c(1,6,2,5,4,3,0)`

**Specify colors**

`cols_s<-'#009FFF,#FF0000,#FF00BF,#FFBF00,#469990,#3CB44B,#BFEF45'`

`cols_s<-unlist(strsplit(cols_s,split=','))`

**Get differential gene table**

`diff_dir<-'/path/to/diff/dir/'`

`fname_list<-list.files(diff_dir,full.names = TRUE)`

`fname_list_stromal<-fname_list[grep('.tsv',fname_list)]`


**Get sorted gene list**

`diff_tab_s<-get_diff_marker_table(fname_list = fname_list_stromal,n = 5,clus_order=clus_order_stromal)`

**Get Heatmap**

`out_dir<-'./'`

`file_prefix_s<-obj_s@project.name`

`print_mean_exp_per_cluster_heatmap(obj_s, assay='integrated', slot = 'scale.data', markers=diff_tab_s$diff_markers, umap_col = cols_s, out=out_dir, file_prefix_s, clus_order=clus_order_stromal, w=11,h=11, cr=FALSE,cc=FALSE)`




