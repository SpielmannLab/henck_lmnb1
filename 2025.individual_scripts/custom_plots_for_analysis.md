# Scripts to make custom plots.

```bash
srun -p shortterm -c 1 -m 200GB --pty bash
module load singularity
singularity shell $WORK/singularity/varunkas-seurat5-plus-1.4.img
workdir=/data/humangen_adld/varun/custom_plots/
mkdir -p ${workdir} && cd ${workdir}
```

1. Now get the Seurat object here as a link, so that no mistakes are made

```bash
#
ln -s /data/humangen_adld/varun/2025.03.01_integratefastMNN_nocovars/clustered/E115-E185-P30-integrate.fastmnn-2500-30_clustering.rds ./sc_obj.rds
```

2. In R, start plotting
   2.1. Create Metadata plots of various QC stuff

```R
sc_obj <- read.rds("sc_obj.rds")
sc_obj$unspliced_ratio <- sc_obj$nCount_unspliced/sc_obj$nCount_spliced

features_to_plot <- c("pct_mt", "Malat1", "nCount_spliced", "nCount_unspliced", "unspliced_ratio", "pct_mt", "pct_rb", "nCount_RNA")

pdf("Metadata_split_by_samples.pdf", height = 10, width = 13)
VlnPlot(sc_obj, features = features_to_plot, group.by = "orig.ident", stack = TRUE, flip = TRUE)
dev.off()


pdf("Metadata_split_by_genotype.pdf", height = 10, width = 6)
VlnPlot(sc_obj, features = features_to_plot, group.by = "sampletype", stack = TRUE, flip = TRUE)
dev.off()

pdf("Metadata_split_by_cluster.pdf", height = 12, width = 15)
VlnPlot(sc_obj, features = features_to_plot, group.by = "CellType", stack = TRUE, flip = TRUE)
dev.off()
```
