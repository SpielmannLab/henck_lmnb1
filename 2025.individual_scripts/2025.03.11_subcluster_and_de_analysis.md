# Scripts to subcluster as needed

```bash
srun -p shortterm -c 1 --mem 200GB --pty bash
module load singularity
singularity shell $WORK/singularity/varunkas-seurat5_plus-1.4.img
export savedir=/data/humangen_adld/varun/2025.03.11_subcluster_and_de_genes/
mkdir -p ${savedir}
```

1.  Now get the Seurat object here as a link, so that no mistakes are made

    ```bash
    #
    cp /data/humangen_adld/varun/2025.03.01_integratefastMNN_nocovars/clustered/E115-E185-P30-integrate.fastmnn-2500-30_clustering.rds $SCRATCH/ && cd $SCRATCH
    ```

2.  Now start subclustering.

    ```R
    suppressPackageStartupMessages(library(Seurat))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(ggplot2))

    # Gather parameters
    file_sc_obj <- list.files(pattern = ".rds")

    # Create a subclustering table. Or read from CSV file
    subcluster_table <- data.frame(RNA_snn_res.0.03 = c("0", "1", "2", "5", "7", "8", "9", "10"),
        resolution = c(0.05, 0.01, 0.005, 0.003, 0.005, 0.03, 0.01, 0.01))
    graph_name = "RNA_snn"
    savedir <- Sys.getenv("savedir")

    # Read the seurat object
    sc_obj <- readRDS(file_sc_obj)

    # Function operates on the global sc_object. Be very careful of overwriting
    sequential_subcluster <- function(cluster, resolution, graph_name) {
        message("Subclustering cluster: ", cluster)
        sc <- FindSubCluster(sc_obj, cluster = cluster,
            graph.name = graph_name, resolution = resolution)
        Idents(sc) <- "sub.cluster"
        sc_obj <<- sc
    }

    Idents(sc_obj) <- names(subcluster_table)[1] # Set the idents to the name of the provided subclustering table

    # Subcluster for all requested values
    purrr::walk2(subcluster_table[,1], subcluster_table[,2], .f = sequential_subcluster, graph_name = graph_name)

    # Save the PDF
     filename <- paste0(savedir, gsub(file_sc_obj, pattern = ".rds", replacement = "_subcluster.pdf"))
    pdf(filename, width = 12, height = 9)
    DimPlot(sc_obj, group.by = "sub.cluster", label = TRUE) +
        theme_void() +
        theme(aspect.ratio = 1)
    dev.off()

    ```

    2.2 Save this seurat object

    ```R

     # This will take some time
     filename <- paste0(gsub(file_sc_obj, pattern = ".rds", replacement = "_subcluster.rds"))
     saveRDS(sc_obj, file = filename)

    ```

    2.3 Move the Seurat object to storage.
    **Now send the R instance to the background using ctrl+z**

    ```bash

    cp *_subcluster.rds ${savedir}

    ```

    2.4 Create a Diet version for plotting etc to save time and memory

    ```R
    sc_obj_diet <- DietSeurat(sc_obj,
        assays = "RNA",
        dimreducs = "umap"
    )

    filename <- paste0(gsub(file_sc_obj, pattern = ".rds", replacement = "_subcluster_diet.rds"))
    saveRDS(sc_obj_diet, file = filename)
    rm(sc_obj_diet)
    ```

    2.3 Move the Diet Seurat object to storage.
    **Now send the R instance to the background using ctrl+z**

    ```bash
    cp *_diet.rds ${savedir}
    ```

3.  Do differential expression on subclusters to find markers and plot markers

    3.1 Do DE

    ```R
    suppressPackageStartupMessages(library(Seurat))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(ggplot2))

    file_sc_obj <- list.files(pattern = "_diet.rds")

    test.use <- "wilcox"
    min.cells.group <- 5
    min.pct <- 0.05
    logfc.threshold <- 0.2
    meta_col <- "sub.cluster"

    # Read the seurat object
    sc_obj <- readRDS(file_sc_obj)

    Idents(sc_obj) <- meta_col

    markers <- FindAllMarkers(sc_obj,
        logfc.threshold = logfc.threshold,
        min.pct = min.pct,
        test.use = test.use,
        min.cells.group = min.cells.group
    )

    # Save the de genes
    savedir <- Sys.getenv("savedir")
    filename <- paste0(savedir, gsub(file_sc_obj, pattern = ".rds", replacement = "_subcluster_markers.tsv"))
    write.table(markers, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)

    ```

    3.2 Plot the DE genes

    ```R

    # Chose top 10 markers for the heatmap
    top_10_markers <- markers %>%
      group_by(cluster) %>%
      dplyr::filter(avg_log2FC > 1) %>%
      slice_head(n = 10) %>%
      ungroup()


    nc <- length(levels(sc_obj))

    # Marker expression visualization (heatmap) [AUC]
    filename <- paste0(savedir, gsub(file_sc_obj, pattern = ".rds", replacement = "_subcluster_heatmap.pdf"))
    pdf(filename, width = 200 * nc / 72, height = 200 * nc / 72)
    sc_obj <- ScaleData(sc_obj, features = top_10_markers$gene)
    DoHeatmap(
      sc_obj,
      features = top_10_markers$gene,
      cells = sample(Cells(sc_obj), size = 6000, replace = FALSE), # Plot only 6000 cells
      group.bar = TRUE,
      group.colors = NULL,
      disp.min = -2.5,
      slot = "scale.data",
      label = TRUE,
      raster = FALSE,
      draw.lines = TRUE,
      combine = TRUE
    ) +
      scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b"))
    dev.off()

    marker <- split(top_10_markers$gene, factor(top_10_markers$cluster, levels = unique(top_10_markers$cluster)))
    print(marker)

    # Plot the list of genes in UMAP space
    filename <- paste0(savedir, gsub(file_sc_obj, pattern = ".rds", replacement = "_subcluster_featureplot_markergenes.pdf"))
    pdf(filename, width = 7 * 4, height = 7 * 3) # 4 genes are plotted per row. 3 per column
    for (i in seq_along(marker)) {
      print(FeaturePlot(
        sc_obj,
        features = marker[[i]], order = TRUE,
        raster = TRUE
      ))
    }
    dev.off()

    ```
