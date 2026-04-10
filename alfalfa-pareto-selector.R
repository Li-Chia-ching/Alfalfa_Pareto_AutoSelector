# ==============================================================================
# Alfalfa Pareto Selector V3.5.1 (March–April 2025)
# - Robust CSV reader for matrix format (rows L1–L50, columns A–U)
# - Traits: March height, April height, multifoliate score (1–5)
# - KNN imputation, Z-score standardization, Pareto ranking + crowding distance
# - NEW: Parallel coordinates visualization with category labels on the plot
# - NEW: Anti-clipping plot margins and 600 DPI publication-ready output
# ==============================================================================

# Robust matrix to long format
# (unchanged code omitted for brevity in explanation — your code structure remains identical)

# Plot 3: Parallel coordinates (with category labels on the plot)

# Prepare data for parallel coordinates plot
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
                        levels = c("Height_Mar", "Height_Apr", "Multi_Score")))

# 1. Compute group mean lines (based on full dataset)
mean_lines <- pcp_data %>%
  group_by(Plot_Group, Trait) %>%
  summarise(Value_Scaled = mean(Value_Scaled, na.rm = TRUE), .groups = "drop")

# 2. Background sampling: retain only 20% of background lines to reduce visual clutter
set.seed(2025)  # ensure reproducibility
bg_data <- pcp_data %>%
  filter(Plot_Group == "3. Background") %>%
  sample_frac(0.2)

fg_data <- pcp_data %>%
  filter(Plot_Group != "3. Background")  # retain all foreground lines

# 3. Label data (for displaying category names on the plot)
label_data <- pcp_data %>%
  filter(Trait == "Height_Mar") %>%
  group_by(Plot_Group) %>%
  summarise(Value_Scaled = mean(Value_Scaled, na.rm = TRUE), .groups = "drop")

# 4. Build the plot (layered rendering)
p3_pcp <- ggplot() +
  # Background layer (sampled, lighter and thinner)
  geom_line(data = bg_data,
            aes(x = Trait, y = Value_Scaled, group = Plant_ID,
                color = Plot_Group, alpha = Plot_Group, linewidth = Plot_Group)) +
  
  # Foreground layer (Rank 1 and selected individuals)
  geom_line(data = fg_data,
            aes(x = Trait, y = Value_Scaled, group = Plant_ID,
                color = Plot_Group, alpha = Plot_Group, linewidth = Plot_Group)) +
  
  # Mean trend lines (thick dashed lines)
  geom_line(data = mean_lines,
            aes(x = Trait, y = Value_Scaled,
                group = Plot_Group, color = Plot_Group),
            linewidth = 1.8, linetype = "dashed",
            inherit.aes = FALSE) +
  
  # Category labels (repelled to avoid overlap)
  geom_text_repel(data = label_data,
                  aes(x = "Height_Mar", y = Value_Scaled,
                      label = Plot_Group),
                  inherit.aes = FALSE,
                  nudge_x = -0.2,
                  direction = "y",
                  hjust = 1,
                  size = 3.5,
                  segment.color = "gray50",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  show.legend = FALSE) +
  
  # Manual scales
  scale_color_manual(values = c("1. Rank 1 Front" = "#E41A1C",
                                "2. Other Selected" = "#377EB8",
                                "3. Background" = "gray80")) +
  
  scale_alpha_manual(values = c("1. Rank 1 Front" = 0.9,
                                "2. Other Selected" = 0.6,
                                "3. Background" = 0.15)) +
  
  scale_linewidth_manual(values = c("1. Rank 1 Front" = 1.5,
                                    "2. Other Selected" = 0.8,
                                    "3. Background" = 0.3)) +
  
  scale_x_discrete(labels = c("March Height",
                              "April Height",
                              "Multi-Score")) +
  
  labs(title = "Multi-Trait Trade-offs (Parallel Coordinates)",
       x = "Selection Traits",
       y = "Min–Max Scaled Value",
       color = "Plant Cohort") +
  
  guides(alpha = "none", linewidth = "none") +
  
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
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
