# ==============================================================================
# Alfalfa Agronomic Population Structure & Diversity Analysis (Pre-Pareto Module)
# Objective: Identify agronomic subgroups and distribution patterns of plant architecture
# Updates:
#   - FIXED: Implemented STRICT 1050 Master Grid to eliminate NA.xxx phantom plants
#   - UPGRADED: Spearman correlation for robust handling of ordinal (Multi_Score) vs continuous (Height) traits
#   - UPGRADED: Terminology shifted to "Agronomic Traits" and "Plant Architecture" for publication readiness
#   - Principal Component Analysis (PCA) for dimensionality reduction
#   - K-means clustering for subpopulation classification
#   - Automated CSV export and 600 DPI publication-quality visualization
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ------------------- 1. Environment & Package Loading -------------------
required_packages <- c("dplyr", "tidyr", "ggplot2", "readr", "factoextra", 
                       "corrplot", "cluster", "stringr", "VIM", "viridis")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "https://cran.rstudio.com/")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(factoextra)
  library(corrplot)
  library(cluster)
  library(stringr)
  library(VIM)
  library(viridis)
})

# Create output directory with timestamp (Updated naming to Agronomic)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- paste0("Agronomic_Structure_Analysis_", timestamp)
dir.create(out_dir, showWarnings = FALSE)

# ------------------- 2. Define Strict Master Grid (1050 Plants) -------------------
# This is the ultimate defense against Excel phantom rows/columns
master_grid <- expand.grid(
  Row_ID = paste0("L", 1:50),
  Col_ID = LETTERS[1:21],
  stringsAsFactors = FALSE
) %>%
  mutate(Plant_ID = paste0(Row_ID, "-", Col_ID)) %>%
  select(Plant_ID, Row_ID, Col_ID)

# ------------------- 3. Robust Matrix-to-Long Data Loading -------------------
matrix_to_long <- function(file_path, value_name) {
  if (!file.exists(file_path)) stop(paste("File not found:", file_path))
  
  # Read raw CSV as character to avoid early coercion issues
  df_raw <- read_csv(file_path, col_names = FALSE, col_types = cols(.default = col_character()))
  
  row_names <- trimws(df_raw[[1]][-1])
  col_names <- trimws(as.character(df_raw[1, -1]))
  data_mat <- df_raw[-1, -1]
  
  colnames(data_mat) <- col_names
  data_mat$Row_ID <- row_names
  
  df_long <- data_mat %>%
    pivot_longer(cols = -Row_ID, names_to = "Col_ID", values_to = "Value") %>%
    mutate(
      Row_ID = trimws(Row_ID),
      Col_ID = trimws(Col_ID),
      Value = suppressWarnings(as.numeric(ifelse(Value %in% c("", "NA", "NULL"), NA, Value)))
    ) %>%
    # STRICT FILTER: Only allow valid identifiers
    filter(Row_ID %in% paste0("L", 1:50), Col_ID %in% LETTERS[1:21]) %>%
    mutate(Plant_ID = paste0(Row_ID, "-", Col_ID)) %>%
    select(Plant_ID, Row_ID, Col_ID, Value) %>%
    rename(!!value_name := Value)
  
  return(df_long)
}

cat("Loading data...\n")
df_mar   <- matrix_to_long("Height_March.csv", "Height_Mar")
df_apr   <- matrix_to_long("Height_April.csv", "Height_Apr")
df_multi <- matrix_to_long("Multifoliate_April.csv", "Multi_Score")

# Merge via Master Grid to guarantee exactly 1050 rows
df_raw <- master_grid %>%
  left_join(df_mar, by = c("Plant_ID", "Row_ID", "Col_ID")) %>%
  left_join(df_apr, by = c("Plant_ID", "Row_ID", "Col_ID")) %>%
  left_join(df_multi, by = c("Plant_ID", "Row_ID", "Col_ID")) %>%
  mutate(Family_ID = Row_ID) %>%
  relocate(Family_ID, .after = Plant_ID)

cat("Total valid plants fixed to experimental design:", nrow(df_raw), "(Target: 1050)\n")

# ------------------- 4. Missing Data Handling -------------------
missing_prop <- df_raw %>%
  summarise(
    Mar_NA = mean(is.na(Height_Mar)),
    Apr_NA = mean(is.na(Height_Apr)),
    Multi_NA = mean(is.na(Multi_Score))
  )

cat("\nMissing proportions:\n"); print(missing_prop)

if (max(missing_prop) < 0.05) {
  cat("Missing rate < 5%: removing rows with missing values.\n")
  df_clean <- df_raw %>% drop_na(Height_Mar, Height_Apr, Multi_Score)
} else {
  cat("Missing rate >= 5%: applying KNN imputation (k=5)...\n")
  set.seed(2025)
  df_to_imp <- as.data.frame(df_raw)
  invisible(capture.output(
    imp_result <- kNN(df_to_imp, k = 5, variable = c("Height_Mar", "Height_Apr", "Multi_Score"), imp_var = FALSE)
  ))
  df_clean <- imp_result
}
stopifnot(all(!is.na(df_clean[, c("Height_Mar", "Height_Apr", "Multi_Score")])))

# ------------------- 5. Z-score Standardization -------------------
df_scaled <- df_clean %>%
  mutate(across(c(Height_Mar, Height_Apr, Multi_Score),
                ~ as.numeric(scale(.)),
                .names = "{.col}_z"))

# ------------------- 6. Agronomic Correlation Analysis (Spearman) -------------------
cat("Generating agronomic correlation matrix using Spearman method...\n")
# UPGRADED: Changed to Spearman to handle ordinal score vs continuous heights properly
cor_matrix <- cor(df_scaled %>% select(ends_with("_z")), method = "spearman")
colnames(cor_matrix) <- c("March Height", "April Height", "Multi-Score")
rownames(cor_matrix) <- colnames(cor_matrix)

pdf(file.path(out_dir, "Plot_01_Correlation_Matrix.pdf"), width = 6, height = 6)
corrplot(cor_matrix, method = "color", type = "upper",
         addCoef.col = "black", tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("#440154FF", "white", "#FDE725FF"))(200), # Viridis-inspired
         title = "Agronomic Trait Correlations", mar = c(0,0,2,0)) # UPGRADED Title
dev.off()

# ------------------- 7. Principal Component Analysis (PCA) -------------------
cat("Performing PCA for plant architecture traits...\n")
pca_res <- prcomp(df_scaled %>% select(ends_with("_z")), center = FALSE, scale. = FALSE)

p_pca <- fviz_pca_biplot(pca_res,
                         geom.ind = "point", pointshape = 21, pointsize = 1.5,
                         fill.ind = "gray80", alpha.ind = 0.5,
                         col.var = "#E41A1C", arrowsize = 0.8, labelsize = 4,
                         repel = TRUE, title = "PCA Biplot of Plant Architecture Traits") +
  theme_bw() + theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(out_dir, "Plot_02_PCA_Biplot.pdf"), p_pca, width = 8, height = 6, dpi = 600)

df_scaled$PC1 <- pca_res$x[, 1]
df_scaled$PC2 <- pca_res$x[, 2]

# ------------------- 8. Optimal K Evaluation & K-means Clustering -------------------
cat("Evaluating optimal cluster number and performing K-means clustering...\n")
p_elbow <- fviz_nbclust(df_scaled %>% select(ends_with("_z")), kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2, color = "red") +
  labs(title = "Optimal Number of Clusters (Elbow Method)",
       subtitle = "Identification of the elbow point for subpopulation structure") +
  theme_bw()

ggsave(file.path(out_dir, "Plot_00_Optimal_K_Elbow.pdf"), p_elbow, width = 6, height = 5, dpi = 600)

optimal_k <- 3
set.seed(2026)
km_res <- kmeans(df_scaled %>% select(ends_with("_z")), centers = optimal_k, nstart = 25)
df_scaled$Subpopulation <- as.factor(paste0("Cluster_", km_res$cluster))

# Plot clusters with robust styling
p_cluster <- ggplot(df_scaled, aes(x = PC1, y = PC2, fill = Subpopulation, color = Subpopulation)) +
  geom_point(shape = 21, alpha = 0.7, size = 2, stroke = 0.5) +
  stat_ellipse(geom = "polygon", alpha = 0.15, linetype = "dashed") +
  scale_fill_viridis_d(option = "cividis") +
  scale_color_viridis_d(option = "cividis") +
  # UPGRADED Title
  labs(title = "Agronomic Subpopulations (K-means)",
       subtitle = "Diversity of plant architecture prior to selection",
       x = paste0("Dim1 (", round(pca_res$sdev[1]^2/sum(pca_res$sdev^2)*100, 1), "%)"),
       y = paste0("Dim2 (", round(pca_res$sdev[2]^2/sum(pca_res$sdev^2)*100, 1), "%)")) +
  theme_bw() + theme(legend.position = "right", plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "Plot_03_Subpopulations.pdf"), p_cluster, width = 8, height = 6, dpi = 600)

# ------------------- 9. Data Export -------------------
cat("Exporting results...\n")
# UPGRADED Filename
write_csv(df_scaled, file.path(out_dir, "Data_01_Agronomic_Structure_Full.csv"))

eig_val <- get_eigenvalue(pca_res)
write.csv(eig_val, file.path(out_dir, "Data_02_PCA_Variance.csv"))

cluster_profile <- df_scaled %>%
  group_by(Subpopulation) %>%
  summarise(
    Count = n(),
    Mean_Height_Mar = mean(Height_Mar),
    Mean_Height_Apr = mean(Height_Apr),
    Mean_Multi_Score = mean(Multi_Score),
    .groups = "drop"
  )
write_csv(cluster_profile, file.path(out_dir, "Data_03_Cluster_Profiles.csv"))

cat("\n========================================================\n")
cat("Pre-Pareto analysis complete. Results saved to:", out_dir, "\n")
cat("========================================================\n")
