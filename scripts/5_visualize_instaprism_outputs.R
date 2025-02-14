##################
## The following script visualizes the raw InstaPrism deconvolution outputs (cell type proportions).
## It also creates and plots a PCA analysis of the InstaPrism deconvolution outputs.
##################

# Load required libraries
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# Define working directory / file directory
file_dir <- file.path("/projects/gakatsu@xsede.org/hgsoc_adipocyte_deconvolution")

# Set seed for reproducibility
set.seed(13)

##################
## 1) Read and arrange the files
##################

# Read files (both with and without adipocytes)
dataset_list <- c("SchildkrautB", "SchildkrautW", "TCGA", "Mayo", "Tothill", "Yoshihara")

for (ds in dataset_list) {
  assign(paste0("instaprism_output_",ds,"_with_adipocytes"),
         fread(file.path(file_dir,paste0("instaprism_output/instaprism_output_",ds,"_with_adipocytes.csv")), header = TRUE, data.table = FALSE))
}

for (ds in dataset_list) {
  assign(paste0("instaprism_output_",ds,"_no_adipocytes"),
         fread(file.path(file_dir,paste0("instaprism_output/instaprism_output_",ds,"_no_adipocytes.csv")), header = TRUE, data.table = FALSE))
}

# Set the values of column 1 as the rownames
for (ds in dataset_list) {
  for (value in c("with","no")){
    assign(paste0("instaprism_output_",ds,"_",value,"_adipocytes"),
           data.frame(get(paste0("instaprism_output_",ds,"_",value,"_adipocytes"))[,-1],
                      row.names=get(paste0("instaprism_output_",ds,"_",value,"_adipocytes"))[,1])) 
  }
}

##################
## 2) Create 100% stacked bar charts to visualize InstaPrism raw output (cell type proportions)
##################

# Find the proportion of each cell type in every sample
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
                                     proportions_TCGA_with_adipocytes, proportions_Mayo_with_adipocytes,
                                     proportions_Tothill_with_adipocytes, proportions_Yoshihara_with_adipocytes)
proportions_no_adipocytes_data <- rbind(proportions_SchildkrautB_no_adipocytes, proportions_SchildkrautW_no_adipocytes,
                                   proportions_TCGA_no_adipocytes, proportions_Mayo_no_adipocytes,
                                   proportions_Tothill_no_adipocytes, proportions_Yoshihara_no_adipocytes)

# Create the stacked bar chart
my_palette <- brewer.pal(9, "RdBu") # 9 is a common max for sequential
my_colors <- colorRampPalette(my_palette)(17)  # Interpolate to get 17 colors

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
dir.create(file.path(file_dir,"instaprism_output/figures"), recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(file_dir,"instaprism_output/figures/cell_proportions_with_adipocytes.pdf"),
       proportions_with_adipocytes, width = 11, height = 8.5, units = "in", device = "pdf")
ggsave(file.path(file_dir,"instaprism_output/figures/cell_proportions_no_adipocytes.pdf"),
       proportions_no_adipocytes, width = 11, height = 8.5, units = "in", device = "pdf")

##################
## 3) Principal component analysis
##################

# Principal component analysis of each dataset
for (ds in dataset_list) {
  for (value in c("with","no")){
    # First perform the PCA 
    assign(paste0("all_pca_",ds,"_",value,"_adipocytes"),
           prcomp(get(paste0("instaprism_output_",ds,"_",value,"_adipocytes")), scale. = TRUE))
    # Then select the top 15 prinicpal components
    assign(paste0("pca_",ds,"_",value,"_adipocytes"),
           get(paste0("all_pca_",ds,"_",value,"_adipocytes"))$x[, 1:15])
  }
}

##################
## 4) K-means clustering for k=3 and k=4
##################

# Perform k-means clustering for k=3 and k=4
for (ds in dataset_list) {
  for (value in c("with","no")){
    # K-means clustering with k=3
    assign(paste0("kmeans_3_",ds,"_",value,"_adipocytes"),
           kmeans(get(paste0("pca_",ds,"_",value,"_adipocytes")), centers = 3))
    # K-means clustering with k=4
    assign(paste0("kmeans_4_",ds,"_",value,"_adipocytes"),
           kmeans(get(paste0("pca_",ds,"_",value,"_adipocytes")), centers = 4))
  }
}

##################
## 5) Get centroids and then center datasets
##################

# Get the centroids
for (ds in dataset_list) {
  for (value in c("with","no")){
    assign(paste0("kmeans_3_centroid_",ds,"_",value,"_adipocytes"),
           as.data.frame(get(paste0("kmeans_3_",ds,"_",value,"_adipocytes"))$centers))
    assign(paste0("kmeans_4_centroid_",ds,"_",value,"_adipocytes"),
           as.data.frame(get(paste0("kmeans_4_",ds,"_",value,"_adipocytes"))$centers))
  }
}

# Calculate the centroid for each dataset (in aggregate) and subtract that from all values
# to set the center of the dataset to (0,0,0,0...). This will center the datasets.
for (ds in dataset_list) {
  for (value in c("with","no")){
    for (k in c(3,4)){
      # First find the coordinates of the dataset's centroid 
      assign(paste0("kmeans_",k,"_coords_",ds,"_",value,"_adipocytes"),
             colMeans(get(paste0("kmeans_",k,"_centroid_",ds,"_",value,"_adipocytes"))))
      # Then subtract the this centroid coordinate from all values in the dataframe
      assign(paste0("kmeans_",k,"_centered_",ds,"_",value,"_adipocytes"),
             get(paste0("kmeans_",k,"_centroid_",ds,"_",value,"_adipocytes")) - get(paste0("kmeans_",k,"_coords_",ds,"_",value,"_adipocytes")))
    }
  }
}

# Add Cluster and Dataset columns to the centroid dataframes
for (ds in dataset_list) {
  for (value in c("with","no")){
    # k=3
    kmeans_3_centroid <- get(paste0("kmeans_3_centered_",ds,"_",value,"_adipocytes"))
    kmeans_3_centroid$Cluster <- 1:3
    kmeans_3_centroid$Dataset <- ds
    assign(paste0("kmeans_3_centered_",ds,"_",value,"_adipocytes"), kmeans_3_centroid)
    
    # k=4
    kmeans_4_centroid <- get(paste0("kmeans_4_centered_",ds,"_",value,"_adipocytes"))
    kmeans_4_centroid$Cluster <- 1:4
    kmeans_4_centroid$Dataset <- ds
    assign(paste0("kmeans_4_centered_",ds,"_",value,"_adipocytes"), kmeans_4_centroid)
  }
}

# Bind the centroid results together
# k=3
kmeans_3_centroid_with_adipocytes <- rbind(kmeans_3_centered_SchildkrautB_with_adipocytes, kmeans_3_centered_SchildkrautW_with_adipocytes,
                                           kmeans_3_centered_TCGA_with_adipocytes, kmeans_3_centered_Mayo_with_adipocytes,
                                           kmeans_3_centered_Tothill_with_adipocytes, kmeans_3_centered_Yoshihara_with_adipocytes)
kmeans_3_centroid_no_adipocytes <- rbind(kmeans_3_centered_SchildkrautB_no_adipocytes, kmeans_3_centered_SchildkrautW_no_adipocytes,
                                         kmeans_3_centered_TCGA_no_adipocytes, kmeans_3_centered_Mayo_no_adipocytes,
                                         kmeans_3_centered_Tothill_no_adipocytes, kmeans_3_centered_Yoshihara_no_adipocytes)

# k=4
kmeans_4_centroid_with_adipocytes <- rbind(kmeans_4_centered_SchildkrautB_with_adipocytes, kmeans_4_centered_SchildkrautW_with_adipocytes,
                                           kmeans_4_centered_TCGA_with_adipocytes, kmeans_4_centered_Mayo_with_adipocytes,
                                           kmeans_4_centered_Tothill_with_adipocytes, kmeans_4_centered_Yoshihara_with_adipocytes)
kmeans_4_centroid_no_adipocytes <- rbind(kmeans_4_centered_SchildkrautB_no_adipocytes, kmeans_4_centered_SchildkrautW_no_adipocytes,
                                         kmeans_4_centered_TCGA_no_adipocytes, kmeans_4_centered_Mayo_no_adipocytes,
                                         kmeans_4_centered_Tothill_no_adipocytes, kmeans_4_centered_Yoshihara_no_adipocytes)

##################
## 6) Create PCA plots
##################

# Define cluster colors
cluster_colors <- c("1" = "indianred2", "2" = "chartreuse3", "3" = "deepskyblue2", "4" = "slateblue1")

# Create centroid plots for the different clusters and datasets (k=3) with adipocytes
gg12_3_with <- ggplot(kmeans_3_centroid_with_adipocytes, aes(x = PC2, y = PC1, color = factor(Cluster), shape = Dataset)) +
  geom_point() +
  theme_bw() + 
  theme(legend.position = "none") +
  xlab("PC2 After Centroid Centering") + 
  ylab("PC1 After Centroid Centering") +
  scale_color_manual(values = cluster_colors[1:3])

gg13_3_with <- ggplot(kmeans_3_centroid_with_adipocytes, aes(x = PC3, y = PC1, color = factor(Cluster), shape = Dataset)) +
  geom_point() +
  theme_bw() +
  xlab("PC3 After Centroid Centering") + 
  ylab("PC1 After Centroid Centering") +
  scale_color_manual(values = cluster_colors[1:3])

gg23_3_with <- ggplot(kmeans_3_centroid_with_adipocytes, aes(x = PC2, y = PC3, color = factor(Cluster), shape = Dataset)) +
  geom_point() +
  theme_bw() + 
  theme(legend.position = "none") +
  xlab("PC2 After Centroid Centering") + 
  ylab("PC3 After Centroid Centering") +
  scale_color_manual(values = cluster_colors[1:3])

gg_3_with = ggarrange(gg12_3_with, gg13_3_with, gg23_3_with, ncol = 3, nrow=1, common.legend = TRUE, legend="bottom")

# Create centroid plots for the different clusters and datasets (k=4) with adipocytes
gg12_4_with <- ggplot(kmeans_4_centroid_with_adipocytes, aes(x = PC2, y = PC1, color = factor(Cluster), shape = Dataset)) +
  geom_point() +
  theme_bw() + 
  theme(legend.position = "none") +
  xlab("PC2 After Centroid Centering") + 
  ylab("PC1 After Centroid Centering") +
  scale_color_manual(values = cluster_colors)

gg13_4_with <- ggplot(kmeans_4_centroid_with_adipocytes, aes(x = PC3, y = PC1, color = factor(Cluster), shape = Dataset)) +
  geom_point() +
  theme_bw() +
  xlab("PC3 After Centroid Centering") + 
  ylab("PC1 After Centroid Centering") +
  scale_color_manual(values = cluster_colors)

gg23_4_with <- ggplot(kmeans_4_centroid_with_adipocytes, aes(x = PC2, y = PC3, color = factor(Cluster), shape = Dataset)) +
  geom_point() +
  theme_bw() + 
  theme(legend.position = "none") +
  xlab("PC2 After Centroid Centering") + 
  ylab("PC3 After Centroid Centering") +
  scale_color_manual(values = cluster_colors)

gg_4_with = ggarrange(gg12_4_with, gg13_4_with, gg23_4_with, ncol = 3, nrow=1, common.legend = TRUE, legend="bottom")

# Creating final 6 plots for kmeans = 3,4 for deconvolution with adipocytes
gg_with = ggarrange(gg_3_with, gg_4_with, ncol = 1, nrow = 2)
gg_with <- annotate_figure(gg_with, top = text_grob("Instaprism deconvolution with adipocyte snRNAseq data", 
                                      color = "slateblue4", face = "bold", size = 14))

# Create centroid plots for the different clusters and datasets (k=3) without adipocytes
gg12_3_no <- ggplot(kmeans_3_centroid_no_adipocytes, aes(x = PC2, y = PC1, color = factor(Cluster), shape = Dataset)) +
  geom_point() +
  theme_bw() + 
  theme(legend.position = "none") +
  xlab("PC2 After Centroid Centering") + 
  ylab("PC1 After Centroid Centering") +
  scale_color_manual(values = cluster_colors[1:3])

gg13_3_no <- ggplot(kmeans_3_centroid_no_adipocytes, aes(x = PC3, y = PC1, color = factor(Cluster), shape = Dataset)) +
  geom_point() +
  theme_bw() +
  xlab("PC3 After Centroid Centering") + 
  ylab("PC1 After Centroid Centering") +
  scale_color_manual(values = cluster_colors[1:3])

gg23_3_no <- ggplot(kmeans_3_centroid_no_adipocytes, aes(x = PC2, y = PC3, color = factor(Cluster), shape = Dataset)) +
  geom_point() +
  theme_bw() + 
  theme(legend.position = "none") +
  xlab("PC2 After Centroid Centering") + 
  ylab("PC3 After Centroid Centering") +
  scale_color_manual(values = cluster_colors[1:3])

gg_3_no = ggarrange(gg12_3_no, gg13_3_no, gg23_3_no, ncol = 3, nrow=1, common.legend = TRUE, legend="bottom")

# Create centroid plots for the different clusters and datasets (k=4) without adipocytes
gg12_4_no <- ggplot(kmeans_4_centroid_no_adipocytes, aes(x = PC2, y = PC1, color = factor(Cluster), shape = Dataset)) +
  geom_point() +
  theme_bw() + 
  theme(legend.position = "none") +
  xlab("PC2 After Centroid Centering") + 
  ylab("PC1 After Centroid Centering") +
  scale_color_manual(values = cluster_colors)

gg13_4_no <- ggplot(kmeans_4_centroid_no_adipocytes, aes(x = PC3, y = PC1, color = factor(Cluster), shape = Dataset)) +
  geom_point() +
  theme_bw() +
  xlab("PC3 After Centroid Centering") + 
  ylab("PC1 After Centroid Centering") +
  scale_color_manual(values = cluster_colors)

gg23_4_no <- ggplot(kmeans_4_centroid_no_adipocytes, aes(x = PC2, y = PC3, color = factor(Cluster), shape = Dataset)) +
  geom_point() +
  theme_bw() + 
  theme(legend.position = "none") +
  xlab("PC2 After Centroid Centering") + 
  ylab("PC3 After Centroid Centering") +
  scale_color_manual(values = cluster_colors)

gg_4_no = ggarrange(gg12_4_no, gg13_4_no, gg23_4_no, ncol = 3, nrow=1, common.legend = TRUE, legend="bottom")

# Creating final 6 plots for kmeans = 3,4 for deconvolution without adipocytes
gg_no = ggarrange(gg_3_no, gg_4_no, ncol = 1, nrow = 2)
gg_no <- annotate_figure(gg_no, top = text_grob("Instaprism deconvolution without adipocyte snRNAseq data", 
                                         color = "slateblue4", face = "bold", size = 14))

# Write final files
dir.create(file.path(file_dir,"instaprism_output/figures"), recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(file_dir,"instaprism_output/figures/centroids_k3_k4_with_adipocytes.pdf"),
       gg_with, width = 11, height = 8.5, units = "in", device = "pdf")
ggsave(file.path(file_dir,"instaprism_output/figures/centroids_k3_k4_no_adipocytes.pdf"),
       gg_no, width = 11, height = 8.5, units = "in", device = "pdf")







