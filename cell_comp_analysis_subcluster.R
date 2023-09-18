# conda activate henck_lmnb1
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(openxlsx)
library(corrplot)
library(forcats)
library(countdata) #package that does betabinomial test on count-data

cellcount_excelfile <-  "202309_Lmnb1_cellcomposition_annotated.xlsx"

# Read the excel sheet with number of cells in main clusters
sheetnames <- getSheetNames(file = cellcount_excelfile)
main_cluster_comp <- read.xlsx(xlsxFile = cellcount_excelfile, sheet = 1, rowNames = TRUE) %>%
    t() %>%
    data.frame() %>%
    replace(is.na(.), 0)

# For each main cluster in the dataframe, replace with subcluster compositions, if it exists in the excel sheet
combined_cluster_comp <- main_cluster_comp
for (main_cluster in (rownames(combined_cluster_comp) %>% .[. %in% sheetnames])) {
    sub_cluster_comp <- read.xlsx(xlsxFile = cellcount_excelfile, sheet = main_cluster, rowNames = TRUE) %>%
        t() %>%
        data.frame() %>%
        mutate(across(everything(), as.integer)) %>%
        replace(is.na(.), 0)

    # Make sure the number of cells in the main cluster is the same as sum of no cells in all subclusters
    if (!identical(colSums(sub_cluster_comp, na.rm=TRUE) %>% as.numeric, combined_cluster_comp[main_cluster,] %>% as.numeric)){
        message("The numbers do not match for the main_Cluster: ", main_cluter)
    }

    # Replace the main cluster cell number with subcluster
    combined_cluster_comp <- combined_cluster_comp %>%
        filter(rownames(combined_cluster_comp) != main_cluster) %>%
        rbind(sub_cluster_comp)
}

# Check the total cell count after merging main_cluster_comp and sub_cluster_comp remains unchanged. Throw an error if not
if(!identical(colSums(main_cluster_comp), colSums(combined_cluster_comp))) stop("There is an error in merging subcluster level values onto the main cluster values")

###### ------------ Do statistical test using VGAM package like CX --------
library(VGAM)
# Requires conversion to dataframe
comp_df <- combined_cluster_comp %>%
    tibble::rownames_to_column(var = "celltype") %>%
    pivot_longer(cols = where(is.numeric), names_to = "sample", values_to = "celltype_n") %>%
    mutate(genotype = gsub(x = sample, pattern = "\\..+\\..+", replacement = "")) %>%
    arrange(desc(genotype)) %>% #Have wild type as the first
    mutate(age = gsub(x = sample, pattern = "\\.\\d*$", replacement = "")) %>%
    mutate(age = gsub(x = age, pattern = "^\\w*\\.", replacement = ""))

size_factor_df <- comp_df %>% 
    group_by(sample) %>% 
    summarize(size_factor = sum(celltype_n))

comp_df <- comp_df %>%
    left_join(size_factor_df, by = "sample") %>%
    mutate(celltype_n_norm = round(celltype_n/size_factor * 1E4)) 

total_norm_df <- comp_df %>%
    group_by(sample) %>% 
    summarize(total_n_norm = sum(celltype_n_norm))

comp_df <- comp_df %>%
    left_join(total_norm_df, by = "sample")

celltype_list <- unique(comp_df$celltype)
genotype_list <- unique(comp_df$genotype)
age_list <- unique(comp_df$age)

# Convert genotype to a factor to prevent alphabetical re-ordering
comp_df$genotype <- factor(comp_df$genotype, levels = genotype_list)

# Do statistical test for each embryonic age and each celltype separately
stat_test_vgam <- list()
for (age_i in age_list) {
    # Initialize a dataframe to store statistical test results
    res <- NULL
    res_error <- NULL

    # For every cell type, do the statistics separately.
    for (i in 1:length(celltype_list)) {
        celltype_i <- celltype_list[i]
        print(paste0(i, " / ", celltype_i))

        comp_df_sub <- comp_df %>%
            filter(celltype == celltype_i) %>%
            filter(age == age_i)

        count_df <- cbind(comp_df_sub$celltype_n_norm, comp_df_sub$total_n_norm - comp_df_sub$celltype_n_norm)

        # Run the regrression
        fit <- tryCatch(
            {vglm(count_df ~ genotype, data = comp_df_sub, family = betabinomial, trace = FALSE)},
            error = function(e) {
                return(NA)
            }
        )
        # Convert the fit results to a usable format
        if (!is.na(fit)) {
            tmp <- data.frame(coef(summary(fit)))
            tmp <- tmp[,c(1,4)]
            names(tmp) <- c("estimate", "pval")
            rownames(tmp) <- gsub('genotype', '', rownames(tmp))
            tmp <- tmp[rownames(tmp) %in% genotype_list,]
            tmp$mutant <- rownames(tmp)
            tmp$celltype <- celltype_i
            rownames(tmp) <- NULL
            res <- rbind(res, tmp)
        } else {
            res_error <- c(res_error, celltype_i)
        }
    }
    stat_test_vgam[[age_i]] <- res %>%
        mutate(fdr = p.adjust(pval, method = "fdr")) %>%
        select(mutant, celltype, fdr) %>%
        pivot_wider(names_from = celltype, values_from = fdr) %>%
        tibble::column_to_rownames(var = "mutant") %>%
        as.matrix()
}

##### ------------ Do statistical test using "countdata" package----------------
# use countdata package to calculate beta-binomial p-values
stat_test_Del.E115 <- combined_cluster_comp %>%
    select(starts_with(c("Del.E115", "WT.E115"))) %>%
    countdata::bb.test(., colSums(.), c("Del", "Del", "WT", "WT")) %>%
    .[["p.value"]] %>%
    p.adjust(method = "BH")
stat_test_Del.E185 <- combined_cluster_comp %>%
    select(starts_with(c("Del.E185", "WT.E185"))) %>%
    countdata::bb.test(., colSums(.), c("Del", "Del", "WT", "WT")) %>%
    .[["p.value"]] %>%
    p.adjust(method = "BH")
stat_test_Dup.E115 <- combined_cluster_comp %>%
    select(starts_with(c("Dup.E115", "WT.E115"))) %>%
    countdata::bb.test(., colSums(.), c("Dup", "Dup", "WT", "WT")) %>%
    .[["p.value"]] %>%
    p.adjust(method = "BH")
stat_test_Dup.E185 <- combined_cluster_comp %>%
    select(starts_with(c("Dup.E185", "WT.E185"))) %>%
    countdata::bb.test(., colSums(.), c("Dup", "Dup", "WT", "WT")) %>%
    .[["p.value"]] %>%
    p.adjust(method = "BH")

#  combine the stat test results based on the time point and add colnames and rownames
stat_test_countdata <- list()
stat_test_countdata[["E115"]] <- cbind(stat_test_Del.E115, stat_test_Dup.E115) %>%
    t()
stat_test_countdata[["E185"]] <- cbind(stat_test_Del.E185, stat_test_Dup.E185) %>%
    t()
colnames(stat_test_countdata[["E115"]]) <- rownames(combined_cluster_comp)
colnames(stat_test_countdata[["E185"]]) <- rownames(combined_cluster_comp)
rownames(stat_test_countdata[["E115"]]) <- gsub(rownames(stat_test_countdata[["E115"]]), pattern = "stat_test_", replacement = "")
rownames(stat_test_countdata[["E185"]]) <- gsub(rownames(stat_test_countdata[["E185"]]), pattern = "stat_test_", replacement = "")
rownames(stat_test_countdata[["E115"]]) <- gsub(rownames(stat_test_countdata[["E115"]]), pattern = ".E115", replacement = "")
rownames(stat_test_countdata[["E185"]]) <- gsub(rownames(stat_test_countdata[["E185"]]), pattern = ".E185", replacement = "")

# ------ Do Plotting ------
# Normalize cell comp and average between repeats
saved_rownames <- rownames(combined_cluster_comp)
combined_cluster_comp <- combined_cluster_comp %>%
    mutate(across(everything(), ~ round(.x/sum(.x)*10000))) %>%
    rowwise() %>% 
    mutate(Del.E115 = mean(c_across(cols=starts_with("Del.E115")))) %>%
    mutate(Del.E185 = mean(c_across(cols=starts_with("Del.E185")))) %>%
    mutate(Dup.E115 = mean(c_across(cols=starts_with("Dup.E115")))) %>%
    mutate(Dup.E185 = mean(c_across(cols=starts_with("Dup.E185")))) %>%
    mutate(WT.E115 = mean(c_across(cols=starts_with("WT.E115")))) %>%
    mutate(WT.E185 = mean(c_across(cols=starts_with("WT.E185")))) %>%
    data.frame() %>%
    select(matches("^\\w*\\.E\\d*$")) %>%
    mutate(across(everything(), round))
# Fix rownames
rownames(combined_cluster_comp) <- saved_rownames

# Remove cell types in which the mean cell number < 10
combined_cluster_comp <- combined_cluster_comp %>%
    mutate(meanE185 = rowMeans(dplyr::select(., ends_with("E185")))) %>%
    mutate(meanE115 = rowMeans(dplyr::select(., ends_with("E115")))) %>%
    filter(meanE185 > 10 & meanE115 > 10) %>%
    select(-meanE185, -meanE115)

# Split the dataframe to E115 and E185
combined_cluster_comp_E185 <- combined_cluster_comp %>%
    select(ends_with("E185"))
combined_cluster_comp_E115 <- combined_cluster_comp %>%
    select(ends_with("E115"))

# Get only the WT composition for reordering and bar plot
combined_cluster_comp_E185_wt <- combined_cluster_comp_E185 %>%
    select(starts_with("WT"))
combined_cluster_comp_E115_wt <- combined_cluster_comp_E115 %>%
    select(starts_with("WT"))

# Calculate log2FC 
combined_cluster_comp_E185 <- combined_cluster_comp_E185 %>%
    mutate(Del.E185 = log2(Del.E185 / WT.E185)) %>%
    mutate(Dup.E185 = log2(Dup.E185 / WT.E185)) %>%
    select(-starts_with("WT"))
combined_cluster_comp_E115 <- combined_cluster_comp_E115 %>%
    mutate(Dup.E115 = log2(Dup.E115 / WT.E115)) %>%
    mutate(Del.E115 = log2(Del.E115 / WT.E115)) %>%
    select(-starts_with("WT"))
# pheatmap(combined_cluster_comp_log2fc, display_numbers = TRUE, number_format = "%.2f", cluster_rows = FALSE, cluster_cols = FALSE, cellwidth = 50, cellheight = 50, filename = "Log2FC_to_WT_subcluster.pdf")
# pheatmap(combined_cluster_comp, display_numbers = TRUE, number_format = "%d", cluster_rows = FALSE, cluster_cols = FALSE, filename = "Size_numbers_averaged_subcluster.pdf")

### Plot using corrplot
# First convert values betweeen -1 and 1
combined_cluster_comp_E185 <- as.matrix(combined_cluster_comp_E185) %>%
    t() %>%
    replace((. > 2), 2) %>%
    replace((. < -2), -2)
rownames(combined_cluster_comp_E185) <- gsub(rownames(combined_cluster_comp_E185),
    pattern = ".\\w*$", 
    replacement = "")
combined_cluster_comp_E115 <- as.matrix(combined_cluster_comp_E115) %>%
    t() %>%
    replace((. > 2), 2) %>%
    replace((. < -2), -2)
rownames(combined_cluster_comp_E115) <- gsub(rownames(combined_cluster_comp_E115),
    pattern = ".\\w*$", 
    replacement = "")

# Re-order the log2fc columns based on the decreasing order of cell composition if required
reorder_clusters = TRUE 
if(reorder_clusters) {
    combined_cluster_comp_E115 <- combined_cluster_comp_E115[, rownames(arrange(combined_cluster_comp_E115_wt, across(everything(), desc)))]
    combined_cluster_comp_E185 <- combined_cluster_comp_E185[, rownames(arrange(combined_cluster_comp_E185_wt, across(everything(), desc)))]
}

# Get the statistical test matrix also in the same column order and row_order
stat_test_vgam[["E115"]] <- stat_test_vgam[["E115"]][rownames(combined_cluster_comp_E115), ]
stat_test_vgam[["E185"]] <- stat_test_vgam[["E185"]][rownames(combined_cluster_comp_E185), ]
stat_test_countdata[["E115"]] <- stat_test_countdata[["E115"]][rownames(combined_cluster_comp_E115), ]
stat_test_countdata[["E185"]] <- stat_test_countdata[["E185"]][rownames(combined_cluster_comp_E185), ]
stat_test_vgam[["E115"]] <- stat_test_vgam[["E115"]][, colnames(combined_cluster_comp_E115)]
stat_test_vgam[["E185"]] <- stat_test_vgam[["E185"]][, colnames(combined_cluster_comp_E185)]
stat_test_countdata[["E115"]] <- stat_test_countdata[["E115"]][, colnames(combined_cluster_comp_E115)]
stat_test_countdata[["E185"]] <- stat_test_countdata[["E185"]][, colnames(combined_cluster_comp_E185)]

if (reorder_clusters) {
    pdf("Log2FC_to_WT_pie_heatmap_subcluster_ordered.pdf", width = 10, height = 10)
} else {
    pdf("Log2FC_to_WT_pie_heatmap_subcluster.pdf", width = 10, height = 10)
}

# Plot the log2fc cell composition for E115
# plot the bar graph of cell numbers in WT
combined_cluster_comp_E115_wt %>% 
    ggplot(aes(x = fct_reorder(row.names(.), WT.E115, .desc = TRUE), y = WT.E115)) + 
    geom_col() +
    scale_y_log10() + 
    theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1),
        aspect.ratio = 0.3) +
    xlab("") + 
    ylab("Number of cells") +
    ggtitle("WT cell composition at E115")
corrplot(corr = combined_cluster_comp_E115,
    p.mat = stat_test_countdata[["E115"]],
    sig.level = -0.05, # set at -0.05 to print all pvalues
    insig = "p-value",
    number.digits = 2,
    method = "pie",
    title = "E115_with_my_stattest_package_countdata",
    is.corr = FALSE,
    mar = c(0,0,2,0),
    col.lim = c(-2, 2))
corrplot(corr = combined_cluster_comp_E115,
    p.mat = stat_test_vgam[["E115"]],
    sig.level = -0.05, # set at -0.05 to print all pvalues
    insig = "p-value",
    number.digits = 2,
    method = "pie",
    title = "E115_with_CX_stattest_package_VGAM",
    is.corr = FALSE,
    mar = c(0,0,2,0),
    col.lim = c(-2, 2))
combined_cluster_comp_E185_wt %>% 
    ggplot(aes(x = fct_reorder(row.names(.), WT.E185, .desc = TRUE), y = WT.E185)) + 
    geom_col() +
    scale_y_log10() + 
    theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1),
        aspect.ratio = 0.3) +
    xlab("") + 
    ylab("Number of cells") +
    ggtitle("WT cell composition at E185")
corrplot(corr = combined_cluster_comp_E185,
    p.mat = stat_test_countdata[["E185"]],
    sig.level = -0.05,# set at -0.05 to print all pvalues
    insig = "p-value",
    number.digits = 2,
    method = "pie",
    title = "E185_with_my_stattest_package_countdata",
    is.corr = FALSE,
    mar = c(0,0,2,0),
    col.lim = c(-2, 2))
corrplot(corr = combined_cluster_comp_E185,
    p.mat = stat_test_vgam[["E185"]],
    sig.level = -0.05,# set at -0.05 to print all pvalues
    insig = "p-value",
    number.digits = 2,
    method = "pie",
    title = "E185_with_CX_stattest_package_VGAM",
    is.corr = FALSE,
    mar = c(0,0,2,0),
    col.lim = c(-2, 2))
dev.off()