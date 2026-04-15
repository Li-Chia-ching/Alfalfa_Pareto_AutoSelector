# ==============================================================================
# Alfalfa Pareto Selector V3.6 (March–April 2025)
# - Robust CSV reader with STRICT 1050 Master Grid (L1-L50 x A-U)
# - Traits: March height, April height, multifoliate score (1–5)
# - KNN imputation, Z-score standardization, Pareto ranking + crowding distance
# - NEW: Parallel coordinates visualization with category labels displayed on the plot
# - NEW: Anti-clipping plot margins and 600 DPI publication-ready output
# ==============================================================================

required_packages <- c("dplyr", "tidyr", "ggplot2", "stringr", "readr",
                       "devEMF", "viridis", "emoa", "VIM", "ggrepel")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "https://cran.rstudio.com/")

suppressPackageStartupMessages({
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
})

# ==================== 1. Define Strict Master Grid (1050 Plants) ====================
# Create a strict 50 x 21 grid to reject any CSV parsing garbage
master_grid <- expand.grid(
  Row_ID = paste0("L", 1:50),
  Col_ID = LETTERS[1:21],
  stringsAsFactors = FALSE
) %>%
  mutate(Plant_ID = paste0(Row_ID, "-", Col_ID)) %>%
  select(Plant_ID, Row_ID, Col_ID)

# ==================== 2. Robust CSV Parser ====================
# Convert matrix-format data to long format and filter out garbage
matrix_to_long <- function(file_path, value_name) {
  if (!file.exists(file_path)) stop(paste("File not found:", file_path))
  
  # Read as characters to avoid parsing failures on weird Excel artifacts
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
      # Convert strictly to numeric, coerce spaces/NULLs to NA
      Value = as.numeric(ifelse(Value %in% c("", "NA", "NULL"), NA, Value))
    ) %>%
    # STRICT FILTER: Keep only recognized rows and columns
    filter(Row_ID %in% paste0("L", 1:50), Col_ID %in% LETTERS[1:21]) %>%
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

# ==================== 3. Force Merge onto Master Grid ====================
# This guarantees exactly 1050 rows. Missing data from CSV becomes clean NAs.
df_raw <- master_grid %>%
  left_join(df_mar, by = c("Plant_ID", "Row_ID", "Col_ID")) %>%
  left_join(df_apr, by = c("Plant_ID", "Row_ID", "Col_ID")) %>%
  left_join(df_multi, by = c("Plant_ID", "Row_ID", "Col_ID"))

cat("Total valid plants fixed to experimental design:", nrow(df_raw), "(Target: 1050)\n")
if(nrow(df_raw) != 1050) warning("Row count is not 1050! Check Master Grid logic.")

# ==================== 4. Handle Missing Data ====================
missing_prop <- df_raw %>%
  summarise(
    Mar_NA = mean(is.na(Height_Mar)),
    Apr_NA = mean(is.na(Height_Apr)),
    Multi_NA = mean(is.na(Multi_Score))
  )
cat("\nMissing proportions:\n"); print(missing_prop)

if (max(missing_prop) < 0.05) {
  cat("Removing rows with missing values (Missing < 5%).\n")
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

# ==================== 5. Standardization & Pareto ====================
standardize <- function(x) (x - mean(x)) / sd(x)
df_std <- df_clean %>%
  mutate(
    Height_Mar_std = standardize(Height_Mar),
    Height_Apr_std = standardize(Height_Apr),
    Multi_Score_std = standardize(Multi_Score)
  )

traits_matrix <- as.matrix(df_std[, c("Height_Mar_std", "Height_Apr_std", "Multi_Score_std")])
df_std$Pareto_Rank <- nds_rank(t(-traits_matrix))

calc_crowding <- function(df_rank) {
  if (nrow(df_rank) < 3) {
    df_rank$Crowd_Dist <- 10  # [Fix] Use 10 instead of Inf to avoid downstream -Inf/NaN
    df_rank$Composite_Score <- (df_rank$Height_Mar_std + df_rank$Height_Apr_std + df_rank$Multi_Score_std) * (1 + 10)
    return(df_rank)
  }
  df_rank <- arrange(df_rank, Height_Apr_std)
  rng <- c(diff(range(df_rank$Height_Mar_std)), diff(range(df_rank$Height_Apr_std)), diff(range(df_rank$Multi_Score_std)))
  rng[rng == 0] <- 1
  n <- nrow(df_rank)
  crowd <- rep(0, n)
  
  for (i in 2:(n-1)) {
    crowd[i] <- (df_rank$Height_Mar_std[i+1] - df_rank$Height_Mar_std[i-1])/rng[1] +
      (df_rank$Height_Apr_std[i+1] - df_rank$Height_Apr_std[i-1])/rng[2] +
      (df_rank$Multi_Score_std[i+1] - df_rank$Multi_Score_std[i-1])/rng[3]
  }
  
  # [Fix] Use 10 instead of Inf for boundary points
  crowd[1] <- crowd[n] <- 10 
  
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

# ==================== 6. Elite Plants Injection (FIXED MATCHING) ====================
spec_list <- unique(c("A10","B30","B9","B14","C10","C13","E16","H4","I19","I21","I25","J26",
                      "M20","Q41","C31","C39","D41","D43","D44","E34","I41","K25","B31","Q42",
                      "S33","T37","P41","D11","D12","F6","H5","H8","C11","C9","D9","K3","N1"))

# Convert "A10" -> Col "A", Row "L10" -> Plant_ID "L10-A"
spec_cols <- substr(spec_list, 1, 1)
spec_rows <- paste0("L", as.numeric(substr(spec_list, 2, nchar(spec_list))))
spec_ids <- unique(paste0(spec_rows, "-", spec_cols))

existing_spec <- intersect(spec_ids, df_scored$Plant_ID)

cat("\nSpecified elite plants found in dataset:", length(existing_spec), "\n")

forced_df <- df_scored %>% filter(Plant_ID %in% existing_spec)
remaining_needed <- 50 - nrow(forced_df)

if (remaining_needed > 0) {
  filler_df <- df_scored %>%
    filter(!Plant_ID %in% existing_spec) %>%
    arrange(Pareto_Rank, desc(Composite_Score)) %>%
    slice_head(n = remaining_needed)
  
  final_selected <- bind_rows(forced_df, filler_df) %>%
    arrange(Pareto_Rank, desc(Composite_Score))
} else {
  final_selected <- forced_df %>%
    arrange(Pareto_Rank, desc(Composite_Score)) %>%
    slice_head(n = 50)
}

cat("Final selection composition:\n")
cat(" - Forced Elites:", sum(final_selected$Plant_ID %in% existing_spec), "\n")
cat(" - Algorithm Selected:", 50 - sum(final_selected$Plant_ID %in% existing_spec), "\n")

df_final <- df_scored %>%
  mutate(Selected_Final = ifelse(Plant_ID %in% final_selected$Plant_ID, "Selected", "Not Selected"))

# ==================== 7. Prepare Data for Plotting ====================
plot_data <- df_final %>%
  mutate(Label = ifelse(Plant_ID %in% existing_spec & Selected_Final == "Selected", Plant_ID, ""))

# Defensive check
if (nrow(plot_data) == 0) stop("No data available for plotting.")
required_cols <- c("Height_Apr", "Multi_Score", "Height_Mar", "Composite_Score", 
                   "Pareto_Rank", "Selected_Final", "Plant_ID")
stopifnot(all(required_cols %in% colnames(plot_data)))

plot_data <- plot_data %>%
  mutate(Pareto_Rank_factor = factor(Pareto_Rank, levels = sort(unique(Pareto_Rank))))

# Pareto front line (Rank 1, sorted by April height)
pareto_front_line <- plot_data %>%
  filter(Pareto_Rank == 1) %>%
  arrange(Height_Apr)

# ==================== 8. Publication-style Theme ====================
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

# ==================== 9. Plot 1: Pareto Front Scatter ====================
p1 <- ggplot(plot_data, aes(x = Height_Apr, y = Multi_Score)) +
  
  # Layer 1: Draw all points, colored by Pareto rank
  # [Fix] Compute safe size directly inside aes() to eliminate "object 'Plot_Size' not found" errors
  geom_point(aes(color = Pareto_Rank_factor, 
                 size = Composite_Score - min(Composite_Score, na.rm = TRUE) + 1), 
             alpha = 0.7, shape = 16) +
  
  # Layer 2: Highlight selected points with a red outline circle
  geom_point(data = filter(plot_data, Selected_Final == "Selected"),
             aes(size = Composite_Score - min(Composite_Score, na.rm = TRUE) + 1),
             shape = 21,        # shape 21 has separate fill and color aesthetics
             fill = NA,         # transparent interior, preserves Pareto color from layer 1
             color = "#FF0000", # pure red border for strong contrast and colorblind visibility
             stroke = 1.5) +    # [Fix] Use stroke to control outline thickness; larger value = thicker circle
  
  # Pareto front trend line layer
  geom_line(data = pareto_front_line, aes(x = Height_Apr, y = Multi_Score), 
            color = "black", linetype = "dashed", linewidth = 1, inherit.aes = FALSE) +
  
  # Color and size scales
  scale_color_viridis_d(option = "magma", direction = -1, name = "Pareto Rank") +
  scale_size_continuous(name = "Composite Score Size", range = c(1, 6)) + 
  labs(
    title = "Pareto Front: April Height vs. Multi-Score",
    subtitle = paste("Selected:", sum(plot_data$Selected_Final == "Selected"), "plants"),
    x = "April Height (original scale)", y = "Multi-Score (1-5)"
  ) +
  academic_theme

# ==================== 10. Plot 2: Growth Dynamics (March vs April Height) ====================
# [Modification] Replace manual colors with viridis cividis palette for colorblind-friendliness

p2 <- ggplot(plot_data, aes(x = Height_Mar, y = Height_Apr)) +
  geom_point(aes(color = Selected_Final, size = Multi_Score), alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  
  # [Modification] Use colorblind‑optimized cividis palette
  scale_color_viridis_d(option = "cividis", name = "Final Selection") +
  scale_size_continuous(name = "Multi-Score") +
  labs(
    title = "March vs. April Height Growth",
    subtitle = "Diagonal line indicates equal height",
    x = "March Height", y = "April Height"
  ) +
  academic_theme

# ==================== 11. Plot 3: Parallel Coordinates ====================
# Data preparation
pcp_data <- plot_data %>%
  select(Plant_ID, Pareto_Rank, Selected_Final, Height_Mar, Height_Apr, Multi_Score) %>%
  mutate(Plot_Group = case_when(
    Pareto_Rank == 1 ~ "1. Rank 1 Front",
    Selected_Final == "Selected" ~ "2. Other Selected",
    TRUE ~ "3. Background"
  )) %>%
  pivot_longer(cols = c(Height_Mar, Height_Apr, Multi_Score),
               names_to = "Trait", values_to = "Value") %>%
  group_by(Trait) %>%
  mutate(Value_Scaled = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(Trait = factor(Trait,
                        levels = c("Height_Mar", "Height_Apr", "Multi_Score"))) %>%
  # [Modification] Ensure Plot_Group is an ordered factor for proper viridis mapping
  mutate(Plot_Group = factor(Plot_Group,
                             levels = c("1. Rank 1 Front", "2. Other Selected", "3. Background")))

# Group mean lines
mean_lines <- pcp_data %>%
  group_by(Plot_Group, Trait) %>%
  summarise(Value_Scaled = mean(Value_Scaled, na.rm = TRUE), .groups = "drop")

# Background sampling (reduce clutter)
set.seed(2025)
bg_data <- pcp_data %>%
  filter(Plot_Group == "3. Background") %>%
  sample_frac(0.2)

fg_data <- pcp_data %>%
  filter(Plot_Group != "3. Background")

# Label positions (at March Height)
label_data <- pcp_data %>%
  filter(Trait == "Height_Mar") %>%
  group_by(Plot_Group) %>%
  summarise(Value_Scaled = mean(Value_Scaled, na.rm = TRUE), .groups = "drop")

# Build parallel coordinates plot
# [Modification] Replace manual colors with viridis plasma palette for colorblind-friendliness

p3_pcp <- ggplot() +
  # Background layer (thin, light)
  geom_line(data = bg_data,
            aes(x = Trait, y = Value_Scaled, group = Plant_ID,
                color = Plot_Group, alpha = Plot_Group, linewidth = Plot_Group)) +
  # Foreground layer (Rank 1 and Selected)
  geom_line(data = fg_data,
            aes(x = Trait, y = Value_Scaled, group = Plant_ID,
                color = Plot_Group, alpha = Plot_Group, linewidth = Plot_Group)) +
  # Mean trend lines (dashed)
  geom_line(data = mean_lines,
            aes(x = Trait, y = Value_Scaled, group = Plot_Group, color = Plot_Group),
            linewidth = 1.8, linetype = "dashed", inherit.aes = FALSE) +
  # Category labels
  geom_text_repel(data = label_data,
                  aes(x = "Height_Mar", y = Value_Scaled, label = Plot_Group),
                  inherit.aes = FALSE,
                  nudge_x = -0.2, direction = "y", hjust = 1,
                  size = 3.5, segment.color = "gray50",
                  box.padding = 0.5, point.padding = 0.3,
                  show.legend = FALSE) +
  
  # [Modification] Use viridis plasma palette
  scale_color_viridis_d(option = "plasma", name = "Plant Cohort", direction = -1) +
  
  scale_alpha_manual(values = c("1. Rank 1 Front" = 0.9,
                                "2. Other Selected" = 0.6,
                                "3. Background" = 0.15)) +
  scale_linewidth_manual(values = c("1. Rank 1 Front" = 1.5,
                                    "2. Other Selected" = 0.8,
                                    "3. Background" = 0.3)) +
  scale_x_discrete(labels = c("March Height", "April Height", "Multi-Score")) +
  labs(title = "Multi-Trait Trade-offs (Parallel Coordinates)",
       x = "Selection Traits", y = "Min–Max Scaled Value",
       color = "Plant Cohort") +
  guides(alpha = "none", linewidth = "none") +
  academic_theme +
  theme(legend.position = "bottom", legend.direction = "horizontal")

# ==================== 12. Save Outputs ====================
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- paste0("Alfalfa_Pareto_Select_", timestamp)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(df_raw, file.path(out_dir, "00_Raw_Merged_Data.csv"))
write_csv(df_final, file.path(out_dir, "01_Full_Processed_Data.csv"))
write_csv(final_selected, file.path(out_dir, "02_Final_Selected_50.csv"))
write_csv(plot_data, file.path(out_dir, "03_Plotting_Data.csv"))
write_csv(pcp_data, file.path(out_dir, "04_Parallel_Coords_Data.csv"))

ggsave(file.path(out_dir, "Plot_Pareto_Selection.pdf"), p1, 
       width = 11, height = 7.5, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Plot_Growth_Dynamics.pdf"), p2, 
       width = 10, height = 7.5, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Plot_Parallel_Coords.pdf"), p3_pcp, 
       width = 9.5, height = 6.5, dpi = 600, limitsize = FALSE)

cat("\nAnalysis complete! Results saved in:", out_dir, "\n")
