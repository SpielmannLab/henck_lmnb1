# conda activate scVelocity
# The script is used to add spliced and unspliced assay counts to the Seurat object, because Jana made a mistake and used the RDS object without these details for downstream analysis
library(Seurat)
library(dplyr)

file_sc_obj <- "/data/humangen_mouse/Lmnb1_Jana/analysis/integrate/E115_E185_25k_30_30/0.1/relable_column/subcluster/cluster/E115_E185_25k_30_30_1-integrate-seurat_annotated.rds"
sc_obj <- readRDS(file_sc_obj)
individual_files <- c("/data/humangen_mouse/Lmnb1_Jana/analysis/E185/seurat/3_seurat_10x_n_velocyto/Del-E185-1-Del-E185_10x_velocyto.rds", 
"/data/humangen_mouse/Lmnb1_Jana/analysis/E185/seurat/3_seurat_10x_n_velocyto/Del-E185-2-Del-E185_10x_velocyto.rds", 
"/data/humangen_mouse/Lmnb1_Jana/analysis/E185/seurat/3_seurat_10x_n_velocyto/Dup-E185-1-Dup-E185_10x_velocyto.rds", 
"/data/humangen_mouse/Lmnb1_Jana/analysis/E185/seurat/3_seurat_10x_n_velocyto/Dup-E185-2-Dup-E185_10x_velocyto.rds", 
"/data/humangen_mouse/Lmnb1_Jana/analysis/E185/seurat/3_seurat_10x_n_velocyto/WT-E185-1-WT-E185_10x_velocyto.rds", 
"/data/humangen_mouse/Lmnb1_Jana/analysis/E185/seurat/3_seurat_10x_n_velocyto/WT-E185-2-WT-E185_10x_velocyto.rds", 
"/data/humangen_mouse/Lmnb1_Jana/analysis/E115/seurat/3_seurat_10x_n_velocyto/Del-E115-1-Del-E115_10x_velocyto.rds", 
"/data/humangen_mouse/Lmnb1_Jana/analysis/E115/seurat/3_seurat_10x_n_velocyto/Del-E115-2-Del-E115_10x_velocyto.rds", 
"/data/humangen_mouse/Lmnb1_Jana/analysis/E115/seurat/3_seurat_10x_n_velocyto/Dup-E115-1-Dup-E115_10x_velocyto.rds", 
"/data/humangen_mouse/Lmnb1_Jana/analysis/E115/seurat/3_seurat_10x_n_velocyto/Dup-E115-2-Dup-E115_10x_velocyto.rds", 
"/data/humangen_mouse/Lmnb1_Jana/analysis/E115/seurat/3_seurat_10x_n_velocyto/WT-E115-1-WT-E115_10x_velocyto.rds", 
"/data/humangen_mouse/Lmnb1_Jana/analysis/E115/seurat/3_seurat_10x_n_velocyto/WT-E115-2-WT-E115_10x_velocyto.rds")

# Load each orig.ident at a time and get spliced and unspliced counts from them
idents_in_data <- unique(sc_obj$orig.ident)
#Genes in RNA assay
DefaultAssay(sc_obj) <- "RNA"

# Gather spliced and unspliced counts matrix from the individual files
genes <- rownames(sc_obj)
# Initialize
spliced <- NULL
unspliced <- NULL
for (ident in idents_in_data) {
    message("Getting spliced and unspliced counts for ident: ", ident)
    individual_file <- grep(pattern = ident, x = individual_files, value = TRUE)
    barcodes <- sc_obj@meta.data %>%
        filter(orig.ident == ident) %>%
        rownames()
    file_id <- gsub(barcodes, pattern = "^.*_", replacement = "") %>% unique()
    cells <- gsub(barcodes, pattern = "_\\d*$", replacement = "")

    individual_sc_obj <- readRDS(individual_file)
    if(!all(cells %in% Cells(individual_sc_obj))) {stop(paste0("Not all cells found in : ", ident))}

    individual_unspliced <- GetAssayData(individual_sc_obj, slot = "counts", assay = "unspliced") %>%
        .[genes, cells]
    individual_spliced <- GetAssayData(individual_sc_obj, slot = "counts", assay = "spliced") %>%
        .[genes, cells]
    colnames(individual_spliced) <- paste0(colnames(individual_spliced), "_", file_id)
    colnames(individual_unspliced) <- paste0(colnames(individual_unspliced), "_", file_id)
    
    spliced <- cbind(spliced, individual_spliced)
    unspliced <- cbind(unspliced, individual_unspliced)
}

spliced_assay <- CreateAssayObject(counts = spliced)
unspliced_assay <- CreateAssayObject(counts = unspliced)

sc_obj[["spliced"]] <- spliced_assay
sc_obj[["unspliced"]] <- unspliced_assay

saveRDS(sc_obj, file = "/data/humangen_mouse/test_area/varun/E115_E185_25k_30_30_1-integrate-seurat_annotated_radialglia.rds")