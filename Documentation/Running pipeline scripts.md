# Ting lab's Seurat based scRNA pipeline


Here are examples of how to run the script

**1) Using raw data from cellranger output**

  **(a) Inputs from Single flowcell**

`/path/to/script/tinglab_scRNA_pipeline.R -i /path/to/cellranger/ouput/flowcell1/ -f file_prefix -o path/to/output/dir/ -s sample1,sample2,sample3,sample4 -t mt,gene_low,gene_high,lib_low,lib_high -c 20 -v -q -u clustree -g all -d 1:30 -e 0.2`

*Notes:*
* Leave out -s argument to load all samples
* Leave out -q and -u to not generate plot qc and clustree (saves time)
* -t defines cell fitering thresholds for high mictochondrial % (mt), low gene count (gene_low), high gene count (gene_high), low library size (lib_low) and high library size (lib_high). Note that even if you don't want to specify one of the cutoffs, the script expects 5 values so e.g. in case of no high mitochondrial cutoff, use 100 for mt cutoff.
* -c is number of cores and change it according to dataset
* -m can be used to specify memory for globals (I use Seurat's suggested parallel processing in this pipeline so need to specify size global object passed to functions; usually set it to 50 or 100GB but may be higher depending on the dataset)
* -g is used to decide whether sample anchors, or all genes are integrated post batch correction. Integrating all genes is problematic for larger datasets.
* By default clustering is done using 1st 50 PCs at resolution 0.5; a different resolution can be specified used -e and different number of PCs using -d like illustrated in the example.

 **(b) Inputs from Multiple flowcells**

`/path/to/script/tinglab_scRNA_pipeline.R -i /path/to/cellranger/ouput/flowcell1/:path/to/cellranger/output/flowcells2/ -f file_prefix -o path/to/output/dir/ -s sample1,sample2,sample3:sample11,sample12,sample13 -t mt,gene_low,gene_high,lib_low,lib_high -c 20 -v -q -u clustree -g all`

*Notes:*
* again leave out -s argument if all samples are to be loaded from both (or more) flowcells

**2) Using Seurat object (Subset analysis)**

`/path/to/script/tinglab_scRNA_pipeline.R -b /path/to/seurat/object/seruat_object.rds -r file_prefix -o path/to/output/dir/ -l 0,1,2,5,6 -c 20 -v -u clustree -g all`

*Notes:*
* -l specifies clusters to be inlcuded in the subset analysis
* Use -t if additional filtering needs to be applied (cannot be lower than the threshold values used for the parent object)
