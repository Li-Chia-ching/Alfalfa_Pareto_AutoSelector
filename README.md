# Alfalfa Pareto AutoSelector

## Overview
This R script performs unsupervised multi-objective selection of alfalfa plants based on three phenotypic traits: early plant height (`Height_early`), late plant height (`Height_late`), and multifoliate score (`Multi_Score`). It uses Pareto non-dominated sorting (from the `emoa` package) to rank plants, computes crowding distance for diversity, and selects the top 50 plants via a composite score. The process is fully automated, with no reliance on prior manual selections.

The algorithm balances trade-offs between maximizing all traits, ensuring a diverse set of high-performing plants without human bias.

## Features
- **Efficient Sorting**: Uses `emoa::nds_rank` for fast Pareto optimization.
- **Vectorized Computations**: Crowding distance and normalization are optimized for speed.
- **Global Imputation**: Missing values are filled with trait means for all plants.
- **Visualization**: Generates a Pareto front plot showing ranks and scores.
- **Outputs**: Saves CSV files for all plants, selected top 50, and plots (PDF and EMF) to a timestamped directory.

## Requirements
- R version 4.0 or higher.
- Required packages: `dplyr`, `tidyr`, `ggplot2`, `stringr`, `readr`, `devEMF`, `viridis`, `emoa`. These are automatically installed if missing.
- Input data frames: `recovering_period_2025`, `initial_flowering_2025`, `Alfalfa_Multi_202504` (assumed to be loaded in the environment; these should contain raw phenotypic data in a grid format with rows as "L1" to "L50" and columns as "A" to "U").

## How to Run
1. Load your phenotypic data frames into the R environment (e.g., via `read_csv` or similar).
2. Source the script: `source("alfalfa_pareto_autoselector.R")`.
3. The script will run automatically, outputting results to a new directory like `Alfalfa_AutoSelect_YYYYMMDD_HHMMSS`.

## Input Data Format
- Each data frame should have the first column as row IDs (e.g., "L1" to "L50") and subsequent columns as column IDs (e.g., "A" to "U") with numeric trait values.
- Plants are identified as "L<row>-<col>" (e.g., "L1-A").

## Output Files
- `00_All_Plants_With_Scores.csv`: All plants with traits, status ("Selected", "Not Selected", or "Dead/Missing"), and composite scores.
- `01_Pareto_Selection_Top50.csv`: Details of the top 50 selected plants.
- `Plot_Pareto_Selection.pdf` / `.emf`: Visualization of the Pareto analysis.
- Console summary: Prints key stats and output directory.

## Limitations
- Assumes a fixed grid of 50 rows x 21 columns; modify `expand.grid` if needed.
- No error handling for invalid data; ensure inputs are numeric and properly formatted.
- For larger datasets, memory usage may increase due to matrix operations.

## Version History
- V1.0: Initial unsupervised version, derived from supervised V4.5 (January 2026).

For questions or contributions, contact the maintainer.
