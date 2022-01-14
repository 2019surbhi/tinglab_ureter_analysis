# Tinglab ureter analysis


## Overview

This repository harbors analysis code and files relevent to scRNA sequencing analysis relevent to our recent work "**Ureter single-cell and spatial mapping reveal cell types, architecture, and signaling networks**". Preprint avaialable [here](https://www.biorxiv.org/content/10.1101/2021.12.22.473889v1)

This code repository utilizes a conventional and popular published scRNA analysis tool- [Seuratv3.2.1](https://www.sciencedirect.com/science/article/pii/S0092867419305598?via%3Dihub), and other available Bioinfirmatics tools, to explore the cellular heterogenity of the human ureters (n=10) using 10X Genomics 3' single cell RNA sequencing platform.

## Pipeline Overview

### Seurat-based pipeline 

The pipeline reads aligned data resulting from the cellranger's standard pipelines - [mkfastq](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#cr-mkfastq) and [count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#cr-count). It runs the following set of analysis using 

<p> (1) Load data and create Seurat object
<p> (2) Quality Control (Filter low quality cells based on abberent mitochondrial gene%, gene counts and UMI counts)
<p> (3) Data Pre-processing (Normalization and Find Variable genes)
<p> (4) Batch correction and data integration (for sample size>1)
<p> (5) Add metadata (age,sex,etc.) - optional (provide a metadata file in xlsx or csv format with appropriate column names)
<p> (6) Pre-clustering processing - Run PCA 
<p> (7) Clustering optimization - Generate Clustree and Silhouette plots (optional)
<p> (7) Clustering - Number of PCs and resolution can be defined else runs with default value of 50 PCs and clusters at resolution 0.5
<p> (8) Post clustering analysis
       <p> (a) Cell proportion table
       <p> (b) Differential gene expression
       <p> (c) Feature plots (for top 9-20 differential markers per cluster)
       <p> (d) Heatmap (for top15 cluster markers)
       <p> (e) Dendrogram

## Content

* [R](https://github.com/2019surbhi/tinglab_ureter_analysis/tree/main/R): `R` package code
  * [tinglab_scRNA_pipeline_functions.R](https://github.com/2019surbhi/tinglab_ureter_analysis/blob/main/R/tinglab_scRNA_pipeline_functions.R) - user defined functions of customized version of standard functions needed to run the pipeline (Please make sure to download this in the currecnt directory or add correct path to this script in the main script to source these functions). Each function provides a brief description of the required and optional variables with their default values also noted.
  *  [tinglab_scRNA_pipeline.R](https://github.com/2019surbhi/tinglab_ureter_analysis/blob/main/R/tinglab_scRNA_pipeline.R) - main script to run the pipeline from start to finish.
* [Documents](https://github.com/2019surbhi/tinglab_ureter_analysis/tree/main/Documentation): Contains instructions on running this pipeline     


