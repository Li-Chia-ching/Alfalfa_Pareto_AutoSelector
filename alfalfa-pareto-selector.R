# ==============================================================================
# Alfalfa Pareto Selector V3.4 (March-April 2025)
# - Robust CSV reader for matrix format (rows L1~L50, cols A~U)
# - Traits: March height, April height, Multifoliate score (1-5)
# - KNN imputation, Z-score, Pareto + crowding distance
# - NEW: Parallel Coordinates visualization for multi-objective trade-offs
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

# Force inclusion of specified plants
spec_list <- unique(c("A10","B30","B9","B14","C10","C13","E16","H4","I19","I21","I25","J26",
                      "M20","Q41","C31","C39","D41","D43","D44","E34","I41","K25","B31","Q42",
                      "S33","T37","P41","D11","D12","F6","H5","H8","C11","C9","D9","K3","N1"))
spec_ids <- unique(paste0(substr(spec_list,1,1), "-", as.numeric(substr(spec_list,2,100))))
existing_spec <- intersect(spec_ids, df_scored$Plant_ID)

if (length(existing_spec) > 0) {
  not_in_initial <- setdiff(existing_spec, initial_selected$Plant_ID)
  if (length(not_in_initial) > 0) {
    combined <- bind_rows(initial_selected, df_scored %>% filter(Plant_ID %in% not_in_initial)) %>%
      distinct(Plant_ID, .keep_all = TRUE) %>%
      arrange(Pareto_Rank, desc(Composite_Score))
    if (nrow(combined) >= 50) {
      final_selected <- combined[1:50, ]
    } else {
      extra <- df_scored %>% filter(!Plant_ID %in% combined$Plant_ID) %>%
        arrange(Pareto_Rank, desc(Composite_Score)) %>% slice_head(n = 50 - nrow(combined))
      final_selected <- bind_rows(combined, extra)
    }
  } else {
    final_selected <- initial_selected
  }
} else {
  final_selected <- initial_selected
}

# Mark final selection
df_final <- df_scored %>%
  mutate(Selected_Final = ifelse(Plant_ID %in% final_selected$Plant_ID, "Selected", "Not Selected"))

plot_data <- df_final %>%
  mutate(Label = ifelse(Plant_ID %in% existing_spec & Selected_Final == "Selected", Plant_ID, ""))

# ==================== Visualization Module (Enhanced Robustness) ====================

# 0. 防御性检查：确保绘图数据有效
if (nrow(plot_data) == 0) stop("No data available for plotting.")
required_cols <- c("Height_Apr", "Multi_Score", "Height_Mar", "Composite_Score", 
                   "Pareto_Rank", "Selected_Final", "Plant_ID")
stopifnot(all(required_cols %in% colnames(plot_data)))

# 1. 确保 Pareto_Rank 为有序因子（保证颜色图例顺序）
plot_data <- plot_data %>%
  mutate(Pareto_Rank_factor = factor(Pareto_Rank, levels = sort(unique(Pareto_Rank))))

# 2. 定义 Pareto 前沿线（Rank 1 按 April Height 排序）
pareto_front_line <- plot_data %>%
  filter(Pareto_Rank == 1) %>%
  arrange(Height_Apr)
if (nrow(pareto_front_line) == 0) {
  warning("No Pareto front (Rank 1) found. Front line will not be drawn.")
}

# 3. 学术主题（右侧图例，垂直排列）
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

# Plot 3: Parallel Coordinates (图例置于底部，因为横向空间更充裕)
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

p3_pcp <- ggplot(pcp_data, aes(x = Trait, y = Value_Scaled, group = Plant_ID, color = Plot_Group)) +
  geom_line(aes(alpha = Plot_Group, linewidth = Plot_Group)) +
  scale_color_manual(values = c("1. Rank 1 Front" = "#E41A1C", 
                                "2. Other Selected" = "#377EB8", 
                                "3. Background" = "gray80")) +
  scale_alpha_manual(values = c("1. Rank 1 Front" = 0.9, "2. Other Selected" = 0.6, "3. Background" = 0.15)) +
  scale_linewidth_manual(values = c("1. Rank 1 Front" = 1.2, "2. Other Selected" = 0.8, "3. Background" = 0.3)) +
  scale_x_discrete(labels = c("March Height", "April Height", "Multi-Score")) +
  labs(title = "Multi-Trait Trade-offs (Parallel Coordinates)",
       x = "Selection Traits", y = "Min-Max Scaled Value",
       color = "Plant Cohort") +
  guides(alpha = "none", linewidth = "none") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.title = element_text(face = "bold", size = 15),
        axis.title = element_text(face = "bold"))

# 输出
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- paste0("Alfalfa_Pareto_Select_", timestamp)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(df_raw, file.path(out_dir, "00_Raw_Merged_Data.csv"))
write_csv(df_final, file.path(out_dir, "01_Full_Processed_Data.csv"))
write_csv(initial_selected, file.path(out_dir, "02_Initial_Top50.csv"))
write_csv(final_selected, file.path(out_dir, "03_Final_Selected_50.csv"))
write_csv(plot_data, file.path(out_dir, "04_Plotting_Data.csv"))
write_csv(pcp_data, file.path(out_dir, "05_Parallel_Coords_Data.csv"))

# 保存图片（确保目录存在，且宽度足够容纳右侧图例）
ggsave(file.path(out_dir, "Plot_Pareto_Selection.pdf"), p1, 
       width = 11, height = 7.5, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Plot_Growth_Dynamics.pdf"), p2, 
       width = 10, height = 7.5, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Plot_Parallel_Coords.pdf"), p3_pcp, 
       width = 9.5, height = 6.5, dpi = 600, limitsize = FALSE)

cat("\nAnalysis complete! Results saved in:", out_dir, "\n")
