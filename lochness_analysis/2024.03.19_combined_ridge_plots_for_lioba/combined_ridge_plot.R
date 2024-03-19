#!/usr/bin/env Rscript

# --- Load all the libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(yaml))

# usage ./combined_ridge_plot.R --args params_yaml_file timepoint out_filename_prefix
args <- commandArgs(trailingOnly = TRUE)
print(args)
params_yaml_file <- args[2]
timepoint <- args[3]
out_filename_prefix <- args[4]

# Read file list from yaml
input_list <- read_yaml(file = "params_file_list.yaml") %>%
  .[[timepoint]]

# Collect all the data into a list
data_list <- lapply(input_list, FUN = function(input) {
  if (!file.exists(input$file)) {
    stop("The file", input$file, "does not exist")
  }
  lochness <- read.delim(input$file) %>%
    mutate(cluster = input$cluster) %>%
    mutate(lochness_comparison = input$comparison) %>%
    filter(condition != "WT") %>%
    select(-contains(c("annotation"))) %>%
    # Rename the subclustering columns with 'subcluster'
    rename_with(.cols = starts_with(c("integrated", "RNA")), .fn = function(colname) {
      gsub(colname, pattern = "(RNA|integrated)_.*", replacement = "subcluster")
    })
})

# Combine them to a long dataframe
data <- do.call(rbind, data_list) %>%
  group_by(cluster, subcluster) %>%
  mutate(size = n()) %>%
  ungroup() %>%
  mutate(cluster_subcluster_size = as.factor(paste(cluster, subcluster, size, sep = ":")))

order_by_size <- levels(data$cluster_subcluster_size) %>%
  gsub(pattern = "\\d:\\d:", replacement = "") %>%
  as.numeric() %>%
  order()
data$cluster_subcluster_size <- factor(data$cluster_subcluster_size, levels = levels(data$cluster_subcluster_size)[order_by_size])

# Make combined ridgeplot
plot <- ggplot(data, aes(x = lochNESS, y = cluster_subcluster_size, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 2) +
  scale_fill_gradient2(
    low = muted("red"),
    mid = "white", high = muted("blue"), midpoint = 0, guide = "colorbar",
    name = "lochNESS score"
  ) +
  theme_classic() +
  theme(
    aspect.ratio = (0.3 + 0.1 * length(unique(data$cluster_subcluster_size))),
    panel.grid.major.y = element_line()
  ) +
  facet_grid(cols = vars(condition))

ggsave(plot,
  filename = paste0(out_filename_prefix, "_lochNESS_plot_ridge_plot.pdf"), width = 8,
  height = 12
)

quit("no")

# Make own density estimates
test_density <- density(data$lochNESS, weights = rep(1, length(data$lochNESS))) %>%
  .[c("x", "y")] %>%
  as.data.frame() %>%
  mutate(group = "full")
# Make own density estimates
test_density_small <- density(sample(data$lochNESS, size = 500), weights = rep(
  1,
  500
)) %>%
  .[c("x", "y")] %>%
  as.data.frame() %>%
  mutate(group = "sampled")
test_density <- rbind(test_density, test_density_small)


ggplot(test_density, aes(x = x, y = group, height = y)) +
  geom_density_ridges(stat = "identity")
