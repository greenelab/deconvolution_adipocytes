##########################################################################################
### 06_visualize_instaprism_output.R
### 
### This script creates several figures to visualize the deconvolution results generated
### in the previous script, and compare the results when run with and without adipocyte
### single nucleus RNA sequencing data in the reference data. It creates 100% stacked bar
### charts to visualize the total cell proportions per bulk dataset, 100% stacked bar
### charts showing cell proportions per sample in each dataset, and bar charts showing
### the absolute change of cell type proportions in total per dataset.
##########################################################################################

# Clear R's memory
rm(list = ls())
gc()

# Load environment
renv::load()

# Load required libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(here)
library(tidyr)

# Set base directories
# Be sure to open this script along with a corresponding project
# that is located in the same directory as the scripts and the 
# input_data folder and enclosed files.
output_data <- file.path(here(), "output_data")
instaprism_folder <- file.path(output_data, "instaprism")
instaprism_figs <- file.path(instaprism_folder, "instaprism_output_visualization")

# Ensure instaprism_figs directory exists
dir.create(instaprism_figs, recursive = TRUE, showWarnings = FALSE)

# Set seed for reproducibility
set.seed(5)

# Setting consistent colors corresponding to the cell types,
# with the first color value representing adipocytes
my_colors <- c(
  "#1f77b4",  # Blue
  "#ff7f0e",  # Orange
  "#2ca02c",  # Green
  "#d62728",  # Red
  "#9467bd",  # Purple
  "#8c564b",  # Brown
  "#e377c2",  # Pink
  "#7f7f7f",  # Gray
  "#bcbd22",  # Olive Green/Lime
  "#17becf",  # Teal/Cyan
  "#a6cee3",  # Light Blue
  "#ffbb78",  # Light Orange
  "#98df8a",  # Light Green
  "#ff9896",  # Light Red
  "#c5b0d5",  # Light Purple
  "#c49c94",  # Tan
  "#f7b6d2"   # Light Pink
)

##########################################################
# 1) Read Instaprism output files
##########################################################

# Read files (both with and without adipocytes)
# These have cell types as columns and samples as rows
dataset_list <- c("SchildkrautB", "SchildkrautW", "TCGA_bulk", "TCGA_microarray", "Tothill", "Yoshihara")

for (ds in dataset_list) {
  assign(paste0("instaprism_output_",ds,"_with_adipocytes"),
         fread(file.path(instaprism_folder,
                         paste0("instaprism_outputs/instaprism_output_",ds,"_with_adipocytes.csv")),
               header = TRUE, data.table = FALSE))
}

for (ds in dataset_list) {
  assign(paste0("instaprism_output_",ds,"_no_adipocytes"),
         fread(file.path(instaprism_folder,
                         paste0("instaprism_outputs/instaprism_output_",ds,"_no_adipocytes.csv")),
               header = TRUE, data.table = FALSE))
}

# Set the values of column 1 (sample IDs) as the rownames
for (ds in dataset_list) {
  for (value in c("with","no")){
    assign(paste0("instaprism_output_",ds,"_",value,"_adipocytes"),
           data.frame(get(paste0("instaprism_output_",ds,"_",value,"_adipocytes"))[,-1],
                      row.names=get(paste0("instaprism_output_",ds,"_",value,"_adipocytes"))[,1])) 
  }
}

##########################################################
# 2) Create 100% stacked bar charts to visualize the
#    total cell proportions per bulk dataset
##########################################################

# Find the proportion of each cell type for both deconvolution results of each bulk dataset
# We will end up with 12 "proportions" objects containing cell types as rownames and total 
# cell proportion values as the first and only column.
for (ds in dataset_list) {
  assign(paste0("sums_",ds,"_with_adipocytes"),
         data.frame(colSums(get(paste0("instaprism_output_",ds,"_with_adipocytes")))))
  assign(paste0("proportions_",ds,"_with_adipocytes"),
         get(paste0("sums_",ds,"_with_adipocytes")) / sum((get(paste0("sums_",ds,"_with_adipocytes"))[,1])))
  assign(paste0("sums_",ds,"_no_adipocytes"),
         data.frame(colSums(get(paste0("instaprism_output_",ds,"_no_adipocytes")))))
  assign(paste0("proportions_",ds,"_no_adipocytes"),
         get(paste0("sums_",ds,"_no_adipocytes")) / sum((get(paste0("sums_",ds,"_no_adipocytes"))[,1])))
}

# Make the colnames (cell types) a column and add a column specifying the dataset
# Now, the proportions objects will contain 3 columns entitled "Proporion",
# "CellType", and "Dataset". They will also still contain the cell types as rownames.
for (ds in dataset_list) {
  for (value in c("with","no")){
    proportions <- get(paste0("proportions_",ds,"_",value,"_adipocytes"))
    colnames(proportions) <- c("Proportion")
    proportions$CellType <- rownames(proportions)
    proportions$Dataset <- ds
    assign(paste0("proportions_",ds,"_",value,"_adipocytes"), proportions)
  }
}

# Bind the dataset proportions together
proportions_with_adipocytes_data <- rbind(proportions_SchildkrautB_with_adipocytes, proportions_SchildkrautW_with_adipocytes,
                                          proportions_TCGA_bulk_with_adipocytes, proportions_TCGA_microarray_with_adipocytes,
                                          proportions_Tothill_with_adipocytes, proportions_Yoshihara_with_adipocytes)
proportions_no_adipocytes_data <- rbind(proportions_SchildkrautB_no_adipocytes, proportions_SchildkrautW_no_adipocytes,
                                        proportions_TCGA_bulk_no_adipocytes, proportions_TCGA_microarray_no_adipocytes,
                                        proportions_Tothill_no_adipocytes, proportions_Yoshihara_no_adipocytes)

# Create the stacked bar chart

proportions_with_adipocytes <- ggplot(proportions_with_adipocytes_data, aes(x = Dataset, y = Proportion, fill = CellType)) +
                                  geom_bar(stat = "identity", position = "stack") +
                                  scale_fill_manual(values = my_colors) +
                                  labs(title = "InstaPrism deconvolution with adipocyte snRNAseq reference data",
                                       x = "Dataset",
                                       y = "Proportion",
                                       fill = "Cell Type") +
                                  theme_bw() +
                                  theme(axis.text.x = element_text(angle = 0, hjust = 1))

proportions_no_adipocytes <- ggplot(proportions_no_adipocytes_data, aes(x = Dataset, y = Proportion, fill = CellType)) +
                                  geom_bar(stat = "identity", position = "stack") +
                                  scale_fill_manual(values = my_colors[-1]) +
                                  labs(title = "InstaPrism deconvolution without adipocyte snRNAseq reference data",
                                       x = "Dataset",
                                       y = "Proportion",
                                       fill = "Cell Type") +
                                  theme_bw() +
                                  theme(axis.text.x = element_text(angle = 0, hjust = 1))

# Write final files
dir.create(file.path(instaprism_figs,"stacked_bar_total"), recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(instaprism_figs,"stacked_bar_total/cell_proportions_with_adipocytes.pdf"),
       proportions_with_adipocytes, width = 11, height = 8.5, units = "in", device = "pdf")
ggsave(file.path(instaprism_figs,"stacked_bar_total/cell_proportions_no_adipocytes.pdf"),
       proportions_no_adipocytes, width = 11, height = 8.5, units = "in", device = "pdf")

##########################################################
# 3) Create 100% stacked bar charts to visualize the
#    cell proportions per sample in each bulk dataset
##########################################################

# Currently the Instaprism output objects have cell types as columns and samples as rows.
# I would like to create a table for each of the 12 objects that has columns "SampleID",
# "CellType", and "Proportion."
make_long_data <- function(ds, value){
  data <- get(paste0("instaprism_output_",ds,"_",value,"_adipocytes"))
  # Get cell types
  cell_types <- colnames(data)
  # Add a column with sample ID
  data$SampleID <- paste0(ds, 1:nrow(data))
  # Use tidyr to pivot longer
  out_data <- pivot_longer(data, cols = all_of(cell_types))
  # Rename the columns
  colnames(out_data) <- c("SampleID","CellType","Proportion")
  return(out_data)
}

for (ds in dataset_list) {
  for (value in c("with","no")){
    assign(paste0("long_instaprism_output_",ds,"_",value,"_adipocytes"),
           make_long_data(ds, value))
  }
}

# Create lists to store the stacked bar charts
per_sample_stacked_no_adipos <- list()
per_sample_stacked_with_adipos <- list()

# Create the charts for Instaprism outputs without adipocytes
for (ds in dataset_list){
  df <- get(paste0("long_instaprism_output_",ds,"_no_adipocytes"))
  
  p <- ggplot(df, aes(x = SampleID, y = Proportion, fill = CellType)) +
    geom_col(position = "fill", color = NA) +
    scale_fill_manual(values = my_colors[-1]) +
    labs(title = paste0("InstaPrism deconvolution without adipocyte snRNAseq reference data: ",ds)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.spacing.x = unit(0, "lines")) +
    labs(x = "Samples", y = "Proportion", fill = "Cell Type")
  
  per_sample_stacked_no_adipos[[ds]] <- p
}

# Create the charts for Instaprism outputs with adipocytes 
for (ds in dataset_list){
  df <- get(paste0("long_instaprism_output_",ds,"_with_adipocytes"))
  
  p <- ggplot(df, aes(x = SampleID, y = Proportion, fill = CellType)) +
    geom_col(position = "fill", color = NA) +
    scale_fill_manual(values = my_colors) +
    labs(title = paste0("InstaPrism deconvolution with adipocyte snRNAseq reference data: ",ds)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.spacing.x = unit(0, "lines")) +
    labs(x = "Samples", y = "Proportion", fill = "Cell Type")
  
  per_sample_stacked_with_adipos[[ds]] <- p
}

# Write final files
dir.create(file.path(instaprism_figs,"stacked_bar_per_sample"), recursive = TRUE, showWarnings = FALSE)
for (ds in dataset_list){
  ggsave(file.path(instaprism_figs,paste0("stacked_bar_per_sample/per_sample_stacked_bar_",ds,"_no_adipocytes.pdf")),
         per_sample_stacked_no_adipos[[ds]],
         width = 11, height = 8.5, , dpi = 1000, units = "in", device = "pdf")
  ggsave(file.path(instaprism_figs,paste0("stacked_bar_per_sample/per_sample_stacked_bar_",ds,"_with_adipocytes.pdf")),
         per_sample_stacked_with_adipos[[ds]],
         width = 11, height = 8.5, , dpi = 1000, units = "in", device = "pdf")
}

##########################################################
# 4) Create bar charts showing the absolute change of cell
#    type proportion in total per dataset before and after
#    adding adipocytes to deconvolution reference data
##########################################################

# Create a dataframe that stores the absolute change of cell proportion per dataset
# before and after adipocyte snRNA seq addition to InstaPrism reference data.
# Of note this is the absolute change (ex. 40% proportion -> 30% proportion is
# expressed as a -10% change, not a -25% change.)
absolute_change <- (proportions_with_adipocytes_data[(proportions_with_adipocytes_data$CellType != "Adipocytes"),"Proportion"] -
                    proportions_no_adipocytes_data[,"Proportion"])
absolute_change_df <- data.frame(proportions_no_adipocytes_data[,c("CellType", "Dataset")])
absolute_change_df$AbsChangeProportion <- absolute_change

# Create list to store the bar charts
abs_change_plots <- list()

# Create bar charts to visualize absolute change of cell types per dataset
for (ds in dataset_list) {
  df_subset <- absolute_change_df %>% filter(Dataset == ds)
  p <- ggplot(df_subset, aes(x = CellType, y = AbsChangeProportion, fill = CellType)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colors[-1]) +
    labs(title = paste0("InstaPrism deconvolution after addition of adipocyte reference: ",ds), x = "Cell Type", y = "Absolute Change in Cell Type Proportion") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(limits = c(-0.25, 0.05))
  
  abs_change_plots[[ds]] <- p
}

# Write final files
dir.create(file.path(instaprism_figs,"absolute_change_cell_types_with_without_adipos"), recursive = TRUE, showWarnings = FALSE)
for (ds in dataset_list){
  ggsave(file.path(instaprism_figs,paste0("absolute_change_cell_types_with_without_adipos/absolute_change_cell_types_with_without_adipos_",ds,".pdf")),
         abs_change_plots[[ds]],
         width = 11, height = 7.5, , dpi = 1000, units = "in", device = "pdf")
}
