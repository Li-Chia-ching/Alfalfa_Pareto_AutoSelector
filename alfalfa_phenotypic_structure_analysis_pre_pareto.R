# Alfalfa Phenotypic Population Structure & Diversity Analysis (Pre-Pareto Module)
# Objective: Identify phenotypic subgroups and distribution patterns of plant architecture traits
# Features:
#   - Correlation matrix for detecting trait redundancy
#   - Principal Component Analysis (PCA) for dimensionality reduction
#   - K-means clustering for subpopulation classification
#   - Automated CSV export and 600 DPI publication-quality visualization
# Updates:
#   - Fixed matrix-column export bug (scale() output issue)
#   - Suppressed expected NA warnings during numeric coercion
#   - Optimized KNN input format (data.frame instead of matrix)

# ------------------- 1. Environment & Package Loading -------------------
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

# ------------------- 2. Robust Matrix-to-Long Data Loading -------------------
# Convert matrix-style CSV (field layout) into tidy long format
matrix_to_long <- function(file_path, value_name) {
  if (!file.exists(file_path)) stop(paste("File not found:", file_path))
  
  # Read raw CSV as character to avoid early coercion issues
  df_raw <- read_csv(file_path, col_names = FALSE, col_types = cols(.default = col_character()))
  
  # Extract row and column identifiers
  row_names <- df_raw[[1]][-1]                # First column (excluding header)
  col_names <- as.character(df_raw[1, -1])    # First row (excluding first cell)
  data_mat <- df_raw[-1, -1]                  # Core data matrix
  
  # Clean and convert to numeric
  data_num <- data_mat %>%
    mutate(across(everything(), ~ ifelse(trimws(.) == "" | . == "NA" | . == "NULL", NA, .))) %>%
    # Suppress expected warnings from numeric coercion
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
    as.matrix()
  
  rownames(data_num) <- row_names
  colnames(data_num) <- col_names
  
  # Convert to long format and construct Plant_ID
  df_long <- data_num %>%
    as.data.frame() %>%
    mutate(Row_ID = rownames(.)) %>%
    pivot_longer(cols = -Row_ID, names_to = "Col_ID", values_to = "Value") %>%
    mutate(Plant_ID = paste0(Row_ID, "-", Col_ID)) %>%
    select(Plant_ID, Row_ID, Col_ID, Value) %>%
    rename(!!value_name := Value)
  
  return(df_long)
}

cat("Loading data using matrix_to_long()...\n")

df_mar   <- matrix_to_long("Height_March.csv", "Height_Mar")
df_apr   <- matrix_to_long("Height_April.csv", "Height_Apr")
df_multi <- matrix_to_long("Multifoliate_April.csv", "Multi_Score")

# Merge all traits
df_raw <- full_join(df_mar, df_apr, by = c("Plant_ID", "Row_ID", "Col_ID")) %>%
  full_join(df_multi, by = c("Plant_ID", "Row_ID", "Col_ID"))

cat("Total number of plants after merging:", nrow(df_raw), "\n")

# Assign Family_ID (row-based grouping)
df_raw <- df_raw %>%
  mutate(Family_ID = Row_ID) %>%
  relocate(Family_ID, .after = Plant_ID)

# ------------------- 3. Missing Data Handling -------------------
# Compute missing rate per trait
missing_prop <- df_raw %>%
  summarise(
    Mar_NA = mean(is.na(Height_Mar)),
    Apr_NA = mean(is.na(Height_Apr)),
    Multi_NA = mean(is.na(Multi_Score))
  )

cat("\nMissing proportions:\n")
print(missing_prop)

# Conditional strategy based on missing rate
if (max(missing_prop) < 0.05) {
  cat("Missing rate < 5%: removing rows with missing values.\n")
  df_clean <- df_raw %>% drop_na(Height_Mar, Height_Apr, Multi_Score)
} else {
  cat("Missing rate >= 5%: applying KNN imputation (k=5)...\n")
  set.seed(2025)
  
  # KNN prefers data.frame input
  df_to_imp <- as.data.frame(df_raw)
  
  invisible(capture.output(
    imp_result <- kNN(df_to_imp,
                     k = 5,
                     variable = c("Height_Mar", "Height_Apr", "Multi_Score"),
                     imp_var = FALSE)
  ))
  
  df_clean <- imp_result
}

# Ensure no missing values remain
stopifnot(all(!is.na(df_clean[, c("Height_Mar", "Height_Apr", "Multi_Score")])))

# ------------------- 4. Z-score Standardization -------------------
# Force numeric conversion to avoid matrix-column issues from scale()
df_scaled <- df_clean %>%
  mutate(across(c(Height_Mar, Height_Apr, Multi_Score),
                ~ as.numeric(scale(.)),
                .names = "{.col}_z"))

# ------------------- 5. Phenotypic Correlation Analysis -------------------
cat("Generating phenotypic correlation matrix...\n")

cor_matrix <- cor(df_scaled %>% select(ends_with("_z")))
colnames(cor_matrix) <- c("March Height", "April Height", "Multi-Score")
rownames(cor_matrix) <- colnames(cor_matrix)

pdf(file.path(out_dir, "Plot_01_Correlation_Matrix.pdf"), width = 6, height = 6)
corrplot(cor_matrix,
         method = "color",
         type = "upper",
         addCoef.col = "black",
         tl.col = "black",
         tl.srt = 45,
         title = "Phenotypic Trait Correlations",
         mar = c(0,0,2,0))
dev.off()

# ------------------- 6. Principal Component Analysis (PCA) -------------------
cat("Performing PCA for plant architecture traits...\n")

pca_res <- prcomp(df_scaled %>% select(ends_with("_z")),
                  center = FALSE,
                  scale. = FALSE)

p_pca <- fviz_pca_biplot(pca_res,
                         geom.ind = "point",
                         pointshape = 21,
                         pointsize = 1.5,
                         fill.ind = "gray70",
                         alpha.ind = 0.6,
                         col.var = "#E41A1C",
                         arrowsize = 0.8,
                         labelsize = 4,
                         repel = TRUE,
                         title = "PCA Biplot of Plant Architecture Traits") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(out_dir, "Plot_02_PCA_Biplot.pdf"), p_pca, width = 8, height = 6, dpi = 600)

# Append PCA scores
df_scaled$PC1 <- pca_res$x[, 1]
df_scaled$PC2 <- pca_res$x[, 2]

# ------------------- 7. Optimal K Evaluation & K-means Clustering -------------------
cat("Evaluating optimal cluster number and performing K-means clustering...\n")

p_elbow <- fviz_nbclust(df_scaled %>% select(ends_with("_z")), kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2, color = "red") +
  labs(title = "Optimal Number of Clusters (Elbow Method)",
       subtitle = "Identification of the elbow point for subpopulation structure") +
  theme_bw()

ggsave(file.path(out_dir, "Plot_00_Optimal_K_Elbow.pdf"), p_elbow, width = 6, height = 5, dpi = 600)

optimal_k <- 3
set.seed(2026)

km_res <- kmeans(df_scaled %>% select(ends_with("_z")),
                 centers = optimal_k,
                 nstart = 25)

df_scaled$Subpopulation <- as.factor(paste0("Cluster_", km_res$cluster))

# Note: Occasional MASS::cov.trob warnings may occur due to outliers when fitting ellipses.
# These do not affect clustering results.
p_cluster <- ggplot(df_scaled, aes(x = PC1, y = PC2, fill = Subpopulation, color = Subpopulation)) +
  geom_point(shape = 21, alpha = 0.7, size = 2, stroke = 0.5) +
  stat_ellipse(geom = "polygon", alpha = 0.15, linetype = "dashed") +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Phenotypic Subpopulations (K-means)",
       subtitle = "Diversity of plant architecture prior to selection",
       x = paste0("Dim1 (", round(pca_res$sdev[1]^2/sum(pca_res$sdev^2)*100, 1), "%)"),
       y = paste0("Dim2 (", round(pca_res$sdev[2]^2/sum(pca_res$sdev^2)*100, 1), "%)")) +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "Plot_03_Subpopulations.pdf"), p_cluster, width = 8, height = 6, dpi = 600)

# ------------------- 8. Data Export -------------------
cat("Exporting results...\n")

write_csv(df_scaled, file.path(out_dir, "Data_01_Pheno_Structure_Full.csv"))

eig_val <- get_eigenvalue(pca_res)
write.csv(eig_val, file.path(out_dir, "Data_02_PCA_Variance.csv"))

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
cat("Pre-Pareto analysis complete. Results saved to:", out_dir, "\n")
cat("========================================================\n")
