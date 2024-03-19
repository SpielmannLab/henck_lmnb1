# Lioba Ridgeplot analysis notes

## Context

Lioba had performed Lochness analysis for the E115 and E185 timepoints separately at the main_cluster level and at the sub_cluster levels. However, the nextflow script currently plots each of them separately, which prohibits good overview, because the scales are so different across these plots. The goal is to plot them all together to see if we can see some patterns.

## File download
The lochness scores with umap co-ordinates were dowloaded to this directory using the following two bash lines

    rsync -rv --filter "- **.pdf" --filter "- **.txt" --filter "- **randomized**" --filter "- /*/*vs*/**" --prune-empty-dirs omics:/data/humangen_singlecell/Lioba/E115_E185_integrate-seurat_2500_20_0.07/subcluster/cell_abundance_analysis_E115 ./
    rsync -rv --filter "- **.pdf" --filter "- **.txt" --filter "- **randomized**" --filter "- /*/*vs*/**" --prune-empty-dirs omics:/data/humangen_singlecell/Lioba/E115_E185_integrate-seurat_2500_20_0.07/subcluster/cell_abundance_analysis_E185 ./ 
    # Note E115 has not cluster_6
    # Since the following files were out of place, they were deleted
    rm -r cell_abundance_analysis_E115/cluster_5/WT_vs*

## Plotting
Check out the combined_ridge_plot.R script
