#! /usr/bin/env Rscript
"Write Seurat Single cell Object as h5ad file for scVelo

Usage: rds_h5ad.R --scobject=<file> 

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --scobject=<file>    *.rds file LoomObject.

"-> doc

suppressMessages(library(docopt))
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)

# --- Dependencies

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))
suppressMessages(library(ggwordcloud))
suppressMessages(library(SeuratDisk))
suppressMessages(library(SeuratWrappers))


# ---------------
# --- Parameters
# ---------------
sc_file <- arguments$scobject
samplename <- gsub(pattern=".rds",replacement="",sc_file)
sc <- readRDS(sc_file)
DefaultAssay(sc) <- "RNA"
print(sc)
SaveH5Seurat(sc, filename=paste0(samplename,".h5Seurat"))
Convert(paste0(samplename,".h5Seurat"), dest="h5ad")
message(paste("The h5ad object",
	      samplename,".h5ad", "has been created."))
