# ==============================================================================
# Alfalfa Pareto Selector V3.5.1 (March-April 2025)
# - Robust CSV reader for matrix format (rows L1~L50, cols A~U)
# - Traits: March height, April height, Multifoliate score (1-5)
# - KNN imputation, Z-score, Pareto + crowding distance
# - NEW: Parallel Coordinates visualization with category labels on the plot
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

# Robust matrix to long format
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

# Load data
cat("Loading data files...\n")
df_mar <- matrix_to_long("Height_March.csv", "Height_Mar")
df_apr <- matrix_to_long("Height_April.csv", "Height_Apr")
df_multi <- matrix_to_long("Multifoliate_April.csv", "Multi_Score")

df_raw <- full_join(df_mar, df_apr, by = c("Plant_ID", "Row_ID", "Col_ID")) %>%
  full_join(df_multi, by = c("Plant_ID", "Row_ID", "Col_ID"))

cat("Total plants:", nrow(df_raw), "\n")

# Missing data handling
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
  cat("Applying KNN imputation (k=5)...\n")
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

# Standardization
standardize <- function(x) (x - mean(x)) / sd(x)
df_std <- df_clean %>%
  mutate(
    Height_Mar_std = standardize(Height_Mar),
    Height_Apr_std = standardize(Height_Apr),
    Multi_Score_std = standardize(Multi_Score)
  )

# Pareto sorting
traits_matrix <- as.matrix(df_std[, c("Height_Mar_std", "Height_Apr_std", "Multi_Score_std")])
df_std$Pareto_Rank <- nds_rank(t(-traits_matrix))

# Crowding distance and composite score
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

# Force inclusion of specified plants (Fixed & Robust Mechanism)
spec_list <- unique(c("A10","B30","B9","B14","C10","C13","E16","H4","I19","I21","I25","J26",
                      "M20","Q41","C31","C39","D41","D43","D44","E34","I41","K25","B31","Q42",
                      "S33","T37","P41","D11","D12","F6","H5","H8","C11","C9","D9","K3","N1"))

# Convert format to "Row-Col"
spec_ids <- unique(paste0(substr(spec_list, 1, 1), "-", as.numeric(substr(spec_list, 2, nchar(spec_list)))))
existing_spec <- intersect(spec_ids, df_scored$Plant_ID)

cat("\nSpecified elite plants found in dataset:", length(existing_spec), "\n")

# Step 1: Extract all specified plants that exist in the dataset unconditionally (reserve slots for preferred)
forced_df <- df_scored %>% filter(Plant_ID %in% existing_spec)

# Step 2: Calculate how many slots remain to reach a total of 50
remaining_needed <- 50 - nrow(forced_df)

if (remaining_needed > 0) {
  # Step 3: Fill remaining slots by selecting the best candidates based on Pareto rank and composite score
  filler_df <- df_scored %>%
    filter(!Plant_ID %in% existing_spec) %>%
    arrange(Pareto_Rank, desc(Composite_Score)) %>%
    slice_head(n = remaining_needed)
  
  # Step 4: Merge the two groups and sort globally for cleaner output
  final_selected <- bind_rows(forced_df, filler_df) %>%
    arrange(Pareto_Rank, desc(Composite_Score))
  
} else {
  # Defensive fallback: if the number of specified plants ≥ 50, keep the top 50 among them
  final_selected <- forced_df %>%
    arrange(Pareto_Rank, desc(Composite_Score)) %>%
    slice_head(n = 50)
}

# (Diagnostic output to verify whether the mechanism is functioning)
cat("Final selection composition:\n")
cat(" - Forced Elites:", sum(final_selected$Plant_ID %in% existing_spec), "\n")
cat(" - Algorithm Selected:", 50 - sum(final_selected$Plant_ID %in% existing_spec), "\n")

# Mark final selection
df_final <- df_scored %>%
  mutate(Selected_Final = ifelse(Plant_ID %in% final_selected$Plant_ID, "Selected", "Not Selected"))

plot_data <- df_final %>%
  mutate(Label = ifelse(Plant_ID %in% existing_spec & Selected_Final == "Selected", Plant_ID, ""))

# ==================== Visualization Module (Enhanced Robustness) ====================

# 0. Defensive check: ensure plotting data is valid
if (nrow(plot_data) == 0) stop("No data available for plotting.")
required_cols <- c("Height_Apr", "Multi_Score", "Height_Mar", "Composite_Score", 
                   "Pareto_Rank", "Selected_Final", "Plant_ID")
stopifnot(all(required_cols %in% colnames(plot_data)))

# 1. Ensure Pareto_Rank is an ordered factor (to preserve legend order)
plot_data <- plot_data %>%
  mutate(Pareto_Rank_factor = factor(Pareto_Rank, levels = sort(unique(Pareto_Rank))))

# 2. Define Pareto front line (Rank 1 sorted by April height)
pareto_front_line <- plot_data %>%
  filter(Pareto_Rank == 1) %>%
  arrange(Height_Apr)
if (nrow(pareto_front_line) == 0) {
  warning("No Pareto front (Rank 1) found. Front line will not be drawn.")
}

# 3. Publication-style theme (legend on the right, vertically stacked)
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

# Plot 1: Pareto Selection
p1 <- ggplot(plot_data, aes(x = Height_Apr, y = Multi_Score)) +
  geom_point(aes(color = Pareto_Rank_factor, size = Composite_Score), alpha = 0.5) +
  { if(nrow(pareto_front_line) > 0) geom_path(data = pareto_front_line, 
                                              aes(x = Height_Apr, y = Multi_Score), color = "#E41A1C", linewidth = 1, alpha = 0.5) } +
  geom_point(data = filter(plot_data, Selected_Final == "Selected"),
             shape = 1, color = "red", size = 4.5, stroke = 1.2) +
  geom_text_repel(aes(label = Label), size = 3.5, box.padding = 0.8, 
                  point.padding = 0.5, max.overlaps = 20, min.segment.length = 0.1) +
  scale_color_viridis_d(name = "Pareto Rank", direction = -1) +
  scale_size_continuous(name = "Composite Score", range = c(1, 5)) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  coord_cartesian(clip = "off") +
  labs(x = "April Height (cm)", y = "Multifoliate Score (1-5)", 
       title = "Pareto Selection: Height vs Multifoliate",
       subtitle = "Red target circles indicate the final 50 selected plants",
       caption = "Triangles represent Pareto Rank 1 Front") +
  academic_theme +
  guides(color = guide_legend(order = 1), size = guide_legend(order = 2))

# Plot 2: Growth Dynamics
p2 <- ggplot(plot_data, aes(x = Height_Mar, y = Height_Apr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
  geom_point(aes(color = Pareto_Rank_factor, size = Composite_Score), alpha = 0.5) +
  geom_point(data = filter(plot_data, Selected_Final == "Selected"),
             shape = 1, color = "red", size = 4.5, stroke = 1.2) +
  scale_color_viridis_d(name = "Pareto Rank", direction = -1) +
  scale_size_continuous(name = "Composite Score", range = c(1, 5)) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  coord_cartesian(clip = "off") +
  labs(x = "March Height (cm)", y = "April Height (cm)", 
       title = "Growth Dynamics: March vs April Height",
       subtitle = "Red target circles indicate the final 50 selected plants") +
  academic_theme

# Plot 3: Parallel coordinates (with category labels on the plot)
# 准备平行坐标数据
pcp_data <- plot_data %>%
  select(Plant_ID, Pareto_Rank, Selected_Final, Height_Mar, Height_Apr, Multi_Score) %>%
  mutate(Plot_Group = case_when(
    Pareto_Rank == 1 ~ "1. Rank 1 Front",
    Selected_Final == "Selected" ~ "2. Other Selected",
    TRUE ~ "3. Background"
  )) %>%
  pivot_longer(cols = c(Height_Mar, Height_Apr, Multi_Score), names_to = "Trait", values_to = "Value") %>%
  group_by(Trait) %>%
  mutate(Value_Scaled = (Value - min(Value, na.rm = TRUE)) / 
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(Trait = factor(Trait, levels = c("Height_Mar", "Height_Apr", "Multi_Score")))

# 1. 计算各组均值线（基于完整数据）
mean_lines <- pcp_data %>%
  group_by(Plot_Group, Trait) %>%
  summarise(Value_Scaled = mean(Value_Scaled, na.rm = TRUE), .groups = "drop")

# 2. 背景采样：只保留 20% 的 Background 线条，减少视觉混乱
set.seed(2025)  # 保证可重复性
bg_data <- pcp_data %>%
  filter(Plot_Group == "3. Background") %>%
  sample_frac(0.2)
fg_data <- pcp_data %>%
  filter(Plot_Group != "3. Background")  # 前景全部保留

# 3. 标签数据（用于在图上显示分类名称）
label_data <- pcp_data %>%
  filter(Trait == "Height_Mar") %>%
  group_by(Plot_Group) %>%
  summarise(Value_Scaled = mean(Value_Scaled, na.rm = TRUE), .groups = "drop")

# 4. 构建图形（分层绘制）
p3_pcp <- ggplot() +
  # 背景层（采样后，透明度高、线宽细）
  geom_line(data = bg_data,
            aes(x = Trait, y = Value_Scaled, group = Plant_ID, color = Plot_Group,
                alpha = Plot_Group, linewidth = Plot_Group)) +
  # 前景层（Rank 1 和 Other Selected，全部保留）
  geom_line(data = fg_data,
            aes(x = Trait, y = Value_Scaled, group = Plant_ID, color = Plot_Group,
                alpha = Plot_Group, linewidth = Plot_Group)) +
  # 均值线层（粗虚线，展示各组中心趋势）
  geom_line(data = mean_lines,
            aes(x = Trait, y = Value_Scaled, group = Plot_Group, color = Plot_Group),
            linewidth = 1.8, linetype = "dashed", inherit.aes = FALSE) +
  # 分类标签（使用 ggrepel 避免重叠）
  geom_text_repel(data = label_data,
                  aes(x = "Height_Mar", y = Value_Scaled, label = Plot_Group),
                  inherit.aes = FALSE,
                  nudge_x = -0.2, direction = "y", hjust = 1,
                  size = 3.5, segment.color = "gray50",
                  box.padding = 0.5, point.padding = 0.3,
                  show.legend = FALSE) +
  # 手动颜色、透明度、线宽映射
  scale_color_manual(values = c("1. Rank 1 Front" = "#E41A1C",
                                "2. Other Selected" = "#377EB8",
                                "3. Background" = "gray80")) +
  scale_alpha_manual(values = c("1. Rank 1 Front" = 0.9,
                                "2. Other Selected" = 0.6,
                                "3. Background" = 0.15)) +
  scale_linewidth_manual(values = c("1. Rank 1 Front" = 1.5,   # 突出 Pareto front
                                    "2. Other Selected" = 0.8,
                                    "3. Background" = 0.3)) +
  scale_x_discrete(labels = c("March Height", "April Height", "Multi-Score")) +
  labs(title = "Multi-Trait Trade-offs (Parallel Coordinates)",
       x = "Selection Traits", y = "Min-Max Scaled Value",
       color = "Plant Cohort") +
  guides(alpha = "none", linewidth = "none") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.title = element_text(face = "bold", size = 15),
        axis.title = element_text(face = "bold"))

# Output
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- paste0("Alfalfa_Pareto_Select_", timestamp)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(df_raw, file.path(out_dir, "00_Raw_Merged_Data.csv"))
write_csv(df_final, file.path(out_dir, "01_Full_Processed_Data.csv"))
write_csv(initial_selected, file.path(out_dir, "02_Initial_Top50.csv"))
write_csv(final_selected, file.path(out_dir, "03_Final_Selected_50.csv"))
write_csv(plot_data, file.path(out_dir, "04_Plotting_Data.csv"))
write_csv(pcp_data, file.path(out_dir, "05_Parallel_Coords_Data.csv"))

# Save plots (ensure directory exists and width accommodates right-side legend)
ggsave(file.path(out_dir, "Plot_Pareto_Selection.pdf"), p1, 
       width = 11, height = 7.5, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Plot_Growth_Dynamics.pdf"), p2, 
       width = 10, height = 7.5, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Plot_Parallel_Coords.pdf"), p3_pcp, 
       width = 9.5, height = 6.5, dpi = 600, limitsize = FALSE)

cat("\nAnalysis complete! Results saved in:", out_dir, "\n")
