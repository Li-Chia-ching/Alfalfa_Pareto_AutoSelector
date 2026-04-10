# ==============================================================================
# Alfalfa Pareto Selector V3.5 (March–April 2025)
# - Robust CSV reader for matrix-format data (rows L1–L50, columns A–U)
# - Traits: March height, April height, multifoliate score (1–5)
# - KNN imputation, Z-score standardization, Pareto ranking + crowding distance
# - NEW: Parallel coordinates visualization with category labels displayed on the plot
# - NEW: Anti-clipping plot margins and 600 DPI publication-ready output
# ==============================================================================

required_packages <- c("dplyr", "tidyr", "ggplot2", "stringr", "readr",
                       "devEMF", "viridis", "emoa", "VIM", "ggrepel")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "https://cran.rstudio.com/")

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(devEMF)
library(viridis)
library(emoa)
library(VIM)
library(ggrepel)

# Convert matrix-format data to long format (robust parser)
matrix_to_long <- function(file_path, value_name) {
  if (!file.exists(file_path)) stop(paste("File not found:", file_path))
  
  df_raw <- read_csv(file_path, col_names = FALSE, col_types = cols(.default = col_character()))
  row_names <- df_raw[[1]][-1]
  col_names <- as.character(df_raw[1, -1])
  data_mat <- df_raw[-1, -1]
  
  data_num <- data_mat %>%
    mutate(across(everything(), ~ ifelse(. == "" | . == "NA" | . == "NULL", NA, .))) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()
  
  rownames(data_num) <- row_names
  colnames(data_num) <- col_names
  
  df_long <- data_num %>%
    as.data.frame() %>%
    mutate(Row_ID = rownames(.)) %>%
    pivot_longer(cols = -Row_ID, names_to = "Col_ID", values_to = "Value") %>%
    mutate(Plant_ID = paste0(Row_ID, "-", Col_ID)) %>%
    select(Plant_ID, Row_ID, Col_ID, Value) %>%
    rename(!!value_name := Value)
  
  return(df_long)
}

# Load input data
cat("Loading data files...\n")
df_mar <- matrix_to_long("Height_March.csv", "Height_Mar")
df_apr <- matrix_to_long("Height_April.csv", "Height_Apr")
df_multi <- matrix_to_long("Multifoliate_April.csv", "Multi_Score")

df_raw <- full_join(df_mar, df_apr, by = c("Plant_ID", "Row_ID", "Col_ID")) %>%
  full_join(df_multi, by = c("Plant_ID", "Row_ID", "Col_ID"))

cat("Total plants:", nrow(df_raw), "\n")

# Handle missing data
missing_prop <- df_raw %>%
  summarise(
    Mar_NA = mean(is.na(Height_Mar)),
    Apr_NA = mean(is.na(Height_Apr)),
    Multi_NA = mean(is.na(Multi_Score))
  )
cat("\nMissing proportions:\n"); print(missing_prop)

if (max(missing_prop) < 0.05) {
  cat("Removing rows with missing values.\n")
  df_clean <- df_raw %>% na.omit()
} else {
  cat("Applying KNN imputation (k = 5)...\n")
  set.seed(2025)
  mat_imp <- df_raw %>% select(Height_Mar, Height_Apr, Multi_Score) %>% as.matrix()
  invisible(capture.output(
    imp_result <- kNN(mat_imp, k = 5, variable = colnames(mat_imp), imp_var = FALSE)
  ))
  df_clean <- df_raw
  df_clean$Height_Mar <- imp_result[, "Height_Mar"]
  df_clean$Height_Apr <- imp_result[, "Height_Apr"]
  df_clean$Multi_Score <- imp_result[, "Multi_Score"]
}
stopifnot(all(!is.na(df_clean[, c("Height_Mar", "Height_Apr", "Multi_Score")])))

# Standardization (Z-score)
standardize <- function(x) (x - mean(x)) / sd(x)
df_std <- df_clean %>%
  mutate(
    Height_Mar_std = standardize(Height_Mar),
    Height_Apr_std = standardize(Height_Apr),
    Multi_Score_std = standardize(Multi_Score)
  )

# Perform Pareto ranking
traits_matrix <- as.matrix(df_std[, c("Height_Mar_std", "Height_Apr_std", "Multi_Score_std")])
df_std$Pareto_Rank <- nds_rank(t(-traits_matrix))

# Compute crowding distance and composite score
calc_crowding <- function(df_rank) {
  if (nrow(df_rank) < 3) {
    df_rank$Crowd_Dist <- Inf
    df_rank$Composite_Score <- 0
    return(df_rank)
  }
  df_rank <- arrange(df_rank, Height_Apr_std)
  rng <- c(diff(range(df_rank$Height_Mar_std)), diff(range(df_rank$Height_Apr_std)), diff(range(df_rank$Multi_Score_std)))
  rng[rng == 0] <- 1
  n <- nrow(df_rank)
  crowd <- rep(0, n)
  crowd[1] <- crowd[n] <- Inf
  for (i in 2:(n-1)) {
    crowd[i] <- (df_rank$Height_Mar_std[i+1] - df_rank$Height_Mar_std[i-1])/rng[1] +
      (df_rank$Height_Apr_std[i+1] - df_rank$Height_Apr_std[i-1])/rng[2] +
      (df_rank$Multi_Score_std[i+1] - df_rank$Multi_Score_std[i-1])/rng[3]
  }
  df_rank$Crowd_Dist <- crowd
  df_rank$Composite_Score <- (df_rank$Height_Mar_std + df_rank$Height_Apr_std + df_rank$Multi_Score_std) * (1 + crowd)
  return(df_rank)
}

df_scored <- df_std %>%
  group_by(Pareto_Rank) %>%
  group_modify(~ calc_crowding(.)) %>%
  ungroup() %>%
  arrange(Pareto_Rank, desc(Composite_Score)) %>%
  mutate(Initial_Rank = row_number())

initial_selected <- df_scored %>% slice_head(n = 50)

# Force inclusion of user-specified elite plants (fixed and robust mechanism)
spec_list <- unique(c("A10","B30","B9","B14","C10","C13","E16","H4","I19","I21","I25","J26",
                      "M20","Q41","C31","C39","D41","D43","D44","E34","I41","K25","B31","Q42",
                      "S33","T37","P41","D11","D12","F6","H5","H8","C11","C9","D9","K3","N1"))

# Convert IDs to "Row-Col" format
spec_ids <- unique(paste0(substr(spec_list, 1, 1), "-", as.numeric(substr(spec_list, 2, nchar(spec_list)))))
existing_spec <- intersect(spec_ids, df_scored$Plant_ID)

cat("\nSpecified elite plants found in dataset:", length(existing_spec), "\n")

# Step 1: Extract all specified plants present in the dataset (reserve slots)
forced_df <- df_scored %>% filter(Plant_ID %in% existing_spec)

# Step 2: Determine remaining slots needed to reach 50 selections
remaining_needed <- 50 - nrow(forced_df)

if (remaining_needed > 0) {
  # Step 3: Fill remaining slots using Pareto rank and composite score
  filler_df <- df_scored %>%
    filter(!Plant_ID %in% existing_spec) %>%
    arrange(Pareto_Rank, desc(Composite_Score)) %>%
    slice_head(n = remaining_needed)
  
  # Step 4: Merge and globally sort results
  final_selected <- bind_rows(forced_df, filler_df) %>%
    arrange(Pareto_Rank, desc(Composite_Score))
  
} else {
  # Fallback: if specified plants ≥ 50, retain the top 50 among them
  final_selected <- forced_df %>%
    arrange(Pareto_Rank, desc(Composite_Score)) %>%
    slice_head(n = 50)
}

# Diagnostic output to verify selection composition
cat("Final selection composition:\n")
cat(" - Forced Elites:", sum(final_selected$Plant_ID %in% existing_spec), "\n")
cat(" - Algorithm Selected:", 50 - sum(final_selected$Plant_ID %in% existing_spec), "\n")

# Annotate final selection status
df_final <- df_scored %>%
  mutate(Selected_Final = ifelse(Plant_ID %in% final_selected$Plant_ID, "Selected", "Not Selected"))

plot_data <- df_final %>%
  mutate(Label = ifelse(Plant_ID %in% existing_spec & Selected_Final == "Selected", Plant_ID, ""))

# ==================== Visualization Module (Enhanced Robustness) ====================

# 0. Defensive check: ensure valid plotting data
if (nrow(plot_data) == 0) stop("No data available for plotting.")
required_cols <- c("Height_Apr", "Multi_Score", "Height_Mar", "Composite_Score", 
                   "Pareto_Rank", "Selected_Final", "Plant_ID")
stopifnot(all(required_cols %in% colnames(plot_data)))

# 1. Ensure Pareto rank is treated as an ordered factor (preserve legend order)
plot_data <- plot_data %>%
  mutate(Pareto_Rank_factor = factor(Pareto_Rank, levels = sort(unique(Pareto_Rank))))

# 2. Define Pareto front line (Rank 1, sorted by April height)
pareto_front_line <- plot_data %>%
  filter(Pareto_Rank == 1) %>%
  arrange(Height_Apr)
if (nrow(pareto_front_line) == 0) {
  warning("No Pareto front (Rank 1) found. The front line will not be drawn.")
}

# 3. Publication-style theme (legend on the right, vertically aligned)
academic_theme <- theme_bw() +
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "#E41A1C", face = "italic", size = 11),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )

# Plot 3: Parallel coordinates (with category labels displayed on the plot)
# Prepare label positions (mean scaled value at March Height for each group)
label_data <- pcp_data %>%
  filter(Trait == "Height_Mar") %>%
  group_by(Plot_Group) %>%
  summarise(Value_Scaled = mean(Value_Scaled, na.rm = TRUE), .groups = "drop")

p3_pcp <- ggplot(pcp_data, aes(x = Trait, y = Value_Scaled, group = Plant_ID, color = Plot_Group)) +
  geom_line(aes(alpha = Plot_Group, linewidth = Plot_Group)) +
  # Key modification: disable aesthetic inheritance; use only label_data columns
  geom_text_repel(
    data = label_data,
    aes(x = "Height_Mar", y = Value_Scaled, label = Plot_Group),
    inherit.aes = FALSE,          # Prevent inheriting Plant_ID from global aesthetics
    nudge_x = -0.2,
    direction = "y",
    hjust = 1,
    size = 3.5,
    segment.color = "gray50",
    box.padding = 0.5,
    point.padding = 0.3,
    show.legend = FALSE
  )
