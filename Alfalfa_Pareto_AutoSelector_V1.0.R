# ==============================================================================
# Alfalfa Pareto AutoSelector V1.0: Unsupervised Pareto Selection
# Efficient non-dominated sorting (emoa package), vectorized crowding distance
# ==============================================================================

# 1. Environment Preparation -----------------------------------------------------------
required_packages <- c("dplyr", "tidyr", "ggplot2", "stringr", "readr",
                       "devEMF", "viridis", "emoa")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos="https://cran.rstudio.com/")
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(devEMF)
library(viridis)
library(emoa)

# 2. Data Loading and Cleaning ---------------------------------------------------------
clean_pheno_data <- function(df_raw, value_name) {
  colnames(df_raw)[1] <- "Row_ID"
  df_long <- df_raw %>%
    filter(str_detect(Row_ID, "^L\\d+$")) %>%
    mutate(across(-Row_ID, as.character)) %>%
    pivot_longer(
      cols = -Row_ID,
      names_to = "Col_ID",
      values_to = "Value",
      values_drop_na = FALSE
    ) %>%
    mutate(
      Value = suppressWarnings(as.numeric(Value)),
      Plant_ID = paste0(Row_ID, "-", Col_ID)
    ) %>%
    select(Plant_ID, Row_ID, Col_ID, Value) %>%
    rename(!!value_name := Value)
  return(df_long)
}

df_analysis <- expand.grid(Row_ID = paste0("L", 1:50), Col_ID = LETTERS[1:21]) %>%
  mutate(Plant_ID = paste0(Row_ID, "-", Col_ID)) %>%
  left_join(clean_pheno_data(recovering_period_2025, "Height_early"), by = c("Row_ID", "Col_ID", "Plant_ID")) %>%
  left_join(clean_pheno_data(initial_flowering_2025, "Height_late"), by = c("Row_ID", "Col_ID", "Plant_ID")) %>%
  left_join(clean_pheno_data(Alfalfa_Multi_202504, "Multi_Score"), by = c("Row_ID", "Col_ID", "Plant_ID"))

# Impute missing values globally using means (unsupervised)
mean_early <- mean(df_analysis$Height_early, na.rm = TRUE)
mean_late <- mean(df_analysis$Height_late, na.rm = TRUE)
mean_multi <- mean(df_analysis$Multi_Score, na.rm = TRUE)
mean_early <- ifelse(is.na(mean_early), 30, mean_early)
mean_late <- ifelse(is.na(mean_late), 50, mean_late)
mean_multi <- ifelse(is.na(mean_multi), 3, mean_multi)
df_analysis <- df_analysis %>%
  mutate(
    Height_early = ifelse(is.na(Height_early), mean_early, Height_early),
    Height_late = ifelse(is.na(Height_late), mean_late, Height_late),
    Multi_Score = ifelse(is.na(Multi_Score), mean_multi, Multi_Score)
  )

# 3. Efficient Pareto Non-Dominated Sorting --------------------------------------------
cat("Starting non-dominated sorting...\n")
start_time <- Sys.time()
df_valid <- df_analysis %>%
  filter(!is.na(Height_early) & !is.na(Height_late) & !is.na(Multi_Score))
if (nrow(df_valid) == 0) stop("No valid data points")
traits_matrix <- as.matrix(df_valid[, c("Height_early", "Height_late", "Multi_Score")])
traits_for_sorting <- -traits_matrix  # Negate to convert to minimization
df_valid$Rank <- nds_rank(t(traits_for_sorting))  # emoa requires transposed matrix
cat("Non-dominated sorting complete. Time taken:", difftime(Sys.time(), start_time, units = "secs"), "seconds. Total fronts:", length(unique(df_valid$Rank)), "\n")

# 4. Crowding Distance and Composite Scoring (Optimized) -------------------------------
cat("Calculating crowding distance and composite scores...\n")
calc_crowding_and_score_fast <- function(df_rank) {
  if (nrow(df_rank) < 3) {
    return(mutate(df_rank,
                  Crowd_Dist = Inf,
                  Composite_Score = 0,
                  Height_early_norm = 0,
                  Height_late_norm = 0,
                  Multi_norm = 0))
  }
  
  # Sort by Height_late
  df_rank <- arrange(df_rank, Height_late)
  
  # Calculate ranges (prevent division by zero)
  ranges <- sapply(c("Height_early", "Height_late", "Multi_Score"), function(col) {
    r <- max(df_rank[[col]]) - min(df_rank[[col]])
    ifelse(r == 0, 1, r)
  })
  names(ranges) <- c("e", "l", "m")
  
  n <- nrow(df_rank)
  crowd_dist <- rep(0, n)
  crowd_dist[c(1, n)] <- Inf
  
  # Vectorized calculation (i+1 - i-1)
  for (trait in c("Height_early", "Height_late", "Multi_Score")) {
    vals <- df_rank[[trait]]
    diffs <- (vals[3:n] - vals[1:(n-2)]) / ranges[substr(trait, 1, 1)]
    crowd_dist[2:(n-1)] <- crowd_dist[2:(n-1)] + diffs
  }
  
  # Normalization and composite score
  df_rank <- df_rank %>%
    mutate(
      Height_early_norm = (Height_early - min(Height_early)) / ranges["e"],
      Height_late_norm = (Height_late - min(Height_late)) / ranges["l"],
      Multi_norm = (Multi_Score - min(Multi_Score)) / ranges["m"],
      Crowd_Dist = crowd_dist,
      Composite_Score = (Height_early_norm + Height_late_norm + Multi_norm) * (1 + Crowd_Dist)
    )
  
  return(df_rank)
}

df_valid <- df_valid %>%
  group_by(Rank) %>%
  group_modify(~calc_crowding_and_score_fast(.)) %>%
  ungroup()

# 5. Select Top 50 (Unsupervised) ------------------------------------------------------
selected <- df_valid %>%
  arrange(Rank, desc(Composite_Score)) %>%
  slice_head(n = 50)

df_analysis <- df_analysis %>%
  mutate(
    Status = case_when(
      Plant_ID %in% selected$Plant_ID ~ "Selected",
      is.na(Height_early) | is.na(Height_late) | is.na(Multi_Score) ~ "Dead/Missing",
      TRUE ~ "Not Selected"
    ),
    Composite_Score = selected$Composite_Score[match(Plant_ID, selected$Plant_ID)]
  )

# 6. Visualization ---------------------------------------------------------------------
p_pareto <- ggplot(df_valid, aes(x = Height_late, y = Multi_Score)) +
  geom_point(aes(color = as.factor(Rank), size = Composite_Score), alpha = 0.85) +
  geom_path(data = df_valid %>% filter(Rank == 1) %>% arrange(Height_late),
            color = "#E41A1C", size = 1, linetype = "solid") +
  scale_color_viridis_d(name = "Pareto Rank", direction = 1) +  # Low ranks darker
  scale_size_continuous(name = "Composite Score", range = c(1, 6)) +
  annotate("segment", x = max(df_valid$Height_late, na.rm = TRUE)*0.85,
           xend = max(df_valid$Height_late, na.rm = TRUE),
           y = 1.5, yend = 2.0,
           arrow = arrow(type = "closed", length = unit(0.3, "cm")),
           color = "black") +
  annotate("text", x = max(df_valid$Height_late, na.rm = TRUE)*0.92, y = 1.4,
           label = "Pareto Trade-off Gradient",
           color = "black", hjust = 0.5, size = 3.5, fontface = "bold") +
  labs(
    title = "Unsupervised Pareto Multi-Trait Selection Analysis",
    subtitle = "Selected Top 50 based on non-dominated sorting and crowding distance",
    x = "Plant Height (Late, cm)",
    y = "Multifoliate Score (1-5)",
    caption = "Points sized by Composite Score."
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.title = element_text(face = "bold", size = 14)
  )

# 7. Save Results ----------------------------------------------------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- paste0("Alfalfa_AutoSelect_", timestamp)
dir.create(out_dir, recursive = TRUE)
write_csv(df_analysis, file.path(out_dir, "00_All_Plants_With_Scores.csv"))
write_csv(selected, file.path(out_dir, "01_Pareto_Selection_Top50.csv"))
ggsave(file.path(out_dir, "Plot_Pareto_Selection.pdf"), p_pareto,
       width = 10, height = 7, dpi = 300, device = "pdf")
emf(file = file.path(out_dir, "Plot_Pareto_Selection.emf"), width = 10, height = 7)
print(p_pareto)
dev.off()

# 8. Output Summary Report -------------------------------------------------------------
cat(paste0("\n", strrep("=", 70), "\n"))
cat("UNSUPERVISED PARETO MULTI-TRAIT SELECTION ANALYSIS COMPLETE\n")
cat(paste0(strrep("=", 70), "\n\n"))
cat("Algorithm Settings:\n")
cat(" • Target Traits: Early Height (Height_early), Late Height (Height_late), Multifoliate Score (Multi_Score)\n")
cat(" • Selection Goal: Top 50 Plants\n")
cat(" • Optimization: emoa Non-Dominated Sorting + Vectorized Crowding Distance\n\n")
cat("Results Statistics:\n")
cat(" • Valid Data Plants:", nrow(df_valid), "\n")
cat(" • Selected Top:", nrow(selected), "\n\n")
cat("Output Directory:", out_dir, "\n")
cat(paste0(strrep("=", 70), "\n"))
