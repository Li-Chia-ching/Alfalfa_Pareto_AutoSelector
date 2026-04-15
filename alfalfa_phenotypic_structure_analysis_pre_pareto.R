# Alfalfa Phenotypic Population Structure & Diversity Analysis (Pre-Pareto Module)
# Objective: Characterize phenotypic substructure and distribution patterns of plant architecture traits
# Key Features:
#   - Trait correlation analysis for redundancy detection
#   - Principal Component Analysis (PCA) for dimensionality reduction
#   - K-means clustering for subpopulation stratification
#   - Automated CSV export and high-resolution (600 DPI) figure generation
# Updates:
#   - Added Family_ID parsing from Plant_ID
#   - Implemented KNN-based missing value imputation
#   - Integrated optimal cluster number (K) evaluation

# ------------------- 1. Environment Setup & Package Loading -------------------
rm(list = ls())
options(stringsAsFactors = FALSE)

required_packages <- c("dplyr", "tidyr", "ggplot2", "readr", "factoextra", "corrplot", "cluster", "stringr", "VIM")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(factoextra)
library(corrplot)
library(cluster)
library(stringr)
library(VIM)

# Create output directory with timestamp
timestamp <- format(Sys.time(), "%Y%m%d")
out_dir <- paste0("Pheno_Structure_Analysis_", timestamp)
dir.create(out_dir, showWarnings = FALSE)

# ------------------- 2. Data Loading, Imputation & Scaling -------------------
cat("Loading and merging phenotypic datasets...\n")

# Read input CSV files (ensure they are in the working directory or provide full paths)
df_mar  <- read_csv("Height_March.csv", show_col_types = FALSE)     # Plant height at regrowth stage (March)
df_apr  <- read_csv("Height_April.csv", show_col_types = FALSE)     # Plant height at early flowering stage (April)
df_multi <- read_csv("Multifoliate_April.csv", show_col_types = FALSE) # Multifoliate trait score

# Standardize column names (assumes two columns: Plant_ID and trait value)
colnames(df_mar)   <- c("Plant_ID", "Height_Mar")
colnames(df_apr)   <- c("Plant_ID", "Height_Apr")
colnames(df_multi) <- c("Plant_ID", "Multi_Score")

# Perform full join to retain all individuals across datasets
df_pheno <- df_mar %>%
  full_join(df_apr, by = "Plant_ID") %>%
  full_join(df_multi, by = "Plant_ID")

cat("Total number of unique individuals after merging:", nrow(df_pheno), "\n")

# Extract Family_ID from Plant_ID (assumes format like "L1-A", prefix before "-" is family)
df_pheno <- df_pheno %>%
  mutate(Family_ID = str_extract(Plant_ID, "^[^-]+")) %>%
  relocate(Family_ID, .after = Plant_ID)

# KNN-based imputation for missing trait values
# Note: If a trait is entirely missing for an individual, prediction is based on distances from other traits
df_imputed <- kNN(
  df_pheno,
  variable = c("Height_Mar", "Height_Apr", "Multi_Score"),
  k = 5,           # Number of nearest neighbors (tunable parameter)
  imp_var = FALSE  # Do not generate additional indicator columns
)

# Z-score standardization of traits
df_scaled <- df_imputed %>%
  mutate(across(c(Height_Mar, Height_Apr, Multi_Score), scale, .names = "{.col}_z"))

# Optional: export merged raw dataset
write_csv(df_pheno, file.path(out_dir, "Data_00_Merged_Raw.csv"))

# ------------------- 3. Phenotypic Correlation Analysis -------------------
cat("Computing phenotypic correlation matrix...\n")

cor_matrix <- cor(df_scaled %>% select(ends_with("_z")))

colnames(cor_matrix) <- c("March Height", "April Height", "Multi-Score")
rownames(cor_matrix) <- colnames(cor_matrix)

pdf(file.path(out_dir, "Plot_01_Correlation_Matrix.pdf"), width = 6, height = 6)
corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45,
  title = "Phenotypic Trait Correlations",
  mar = c(0,0,2,0)
)
dev.off()

# ------------------- 4. Principal Component Analysis (PCA) -------------------
cat("Performing PCA for plant architecture traits...\n")

pca_res <- prcomp(
  df_scaled %>% select(ends_with("_z")),
  center = FALSE,
  scale. = FALSE
)

# PCA biplot visualization
p_pca <- fviz_pca_biplot(
  pca_res,
  geom.ind = "point",
  pointshape = 21,
  pointsize = 1.5,
  fill.ind = "gray70",
  alpha.ind = 0.6,
  col.var = "#E41A1C",
  arrowsize = 0.8,
  labelsize = 4,
  repel = TRUE,
  title = "PCA Biplot of Plant Architecture Traits"
) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(out_dir, "Plot_02_PCA_Biplot.pdf"), p_pca, width = 8, height = 6, dpi = 600)

# Append PCA scores to dataset
df_scaled$PC1 <- pca_res$x[, 1]
df_scaled$PC2 <- pca_res$x[, 2]

# ------------------- 5. Optimal K Selection & K-means Clustering -------------------
cat("Evaluating optimal cluster number and performing K-means clustering...\n")

# Elbow method for determining optimal K
p_elbow <- fviz_nbclust(
  df_scaled %>% select(ends_with("_z")),
  kmeans,
  method = "wss"
) +
  geom_vline(xintercept = 3, linetype = 2, color = "red") +
  labs(
    title = "Optimal Number of Clusters (Elbow Method)",
    subtitle = "Identification of the elbow point for subpopulation stratification"
  ) +
  theme_bw()

ggsave(file.path(out_dir, "Plot_00_Optimal_K_Elbow.pdf"), p_elbow, width = 6, height = 5, dpi = 600)

# Perform K-means clustering (adjust K based on elbow plot)
optimal_k <- 3
set.seed(2026)

km_res <- kmeans(
  df_scaled %>% select(ends_with("_z")),
  centers = optimal_k,
  nstart = 25
)

df_scaled$Subpopulation <- as.factor(paste0("Cluster_", km_res$cluster))

# Visualization of clusters in PCA space
p_cluster <- ggplot(df_scaled, aes(x = PC1, y = PC2, fill = Subpopulation, color = Subpopulation)) +
  geom_point(shape = 21, alpha = 0.7, size = 2, stroke = 0.5) +
  stat_ellipse(geom = "polygon", alpha = 0.15, linetype = "dashed") +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Phenotypic Subpopulations (K-means)",
    subtitle = "Diversity of plant architecture prior to selection",
    x = paste0("Dim1 (", round(pca_res$sdev[1]^2/sum(pca_res$sdev^2)*100, 1), "%)"),
    y = paste0("Dim2 (", round(pca_res$sdev[2]^2/sum(pca_res$sdev^2)*100, 1), "%)")
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "Plot_03_Subpopulations.pdf"), p_cluster, width = 8, height = 6, dpi = 600)

# ------------------- 6. Data Export -------------------
cat("Exporting analysis outputs...\n")

# Full dataset including metadata, traits, PCA scores, and cluster assignments
write_csv(df_scaled, file.path(out_dir, "Data_01_Pheno_Structure_Full.csv"))

# PCA eigenvalues and explained variance
eig_val <- get_eigenvalue(pca_res)
write.csv(eig_val, file.path(out_dir, "Data_02_PCA_Variance.csv"))

# Cluster-level summary statistics
cluster_profile <- df_scaled %>%
  group_by(Subpopulation) %>%
  summarise(
    Count = n(),
    Mean_Height_Mar = mean(Height_Mar),
    Mean_Height_Apr = mean(Height_Apr),
    Mean_Multi_Score = mean(Multi_Score)
  )

write_csv(cluster_profile, file.path(out_dir, "Data_03_Cluster_Profiles.csv"))

cat("\n========================================================\n")
cat("Pre-Pareto analysis completed. Outputs saved to:", out_dir, "\n")
cat("========================================================\n")
