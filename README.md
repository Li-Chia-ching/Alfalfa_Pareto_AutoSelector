# Alfalfa Pareto Selector V3.5.1.

## Overview

Alfalfa Pareto Selector V3.4 is an R-based pipeline for multi-trait selection in alfalfa (Medicago sativa), integrating data preprocessing, multi-objective optimization, and publication-quality visualization.

The framework is designed for matrix-format field data and supports robust selection of elite individuals using Pareto ranking and diversity-aware scoring.

---

## Features

- Robust matrix-format CSV parsing (L1–L50 × A–U layout)
- Automatic missing data handling:
  - Row removal (<5% missing)
  - KNN imputation (k = 5)
- Z-score standardization
- Multi-objective optimization:
  - Pareto ranking (non-dominated sorting)
  - Crowding distance
- Composite scoring system
- Forced inclusion of elite plants
- Publication-ready visualization:
  - Pareto scatter plot
  - Growth dynamics plot
  - Parallel coordinates plot
- Anti-clipping layout and 600 DPI export

---

## Input Data

Three CSV files are required:

| File | Description |
|------|------------|
| Height_March.csv | March plant height |
| Height_April.csv | April plant height |
| Multifoliate_April.csv | Multifoliate score (1–5) |

### Format

- First row: column IDs (A–U)
- First column: row IDs (L1–L50)
- Values: numeric

---

## Workflow

1. Data import and reshaping
2. Missing data handling
3. Z-score normalization
4. Pareto ranking
5. Crowding distance calculation
6. Elite selection (Top 50 + forced inclusion)
7. Visualization

---

## Output

A timestamped directory will be created:
Alfalfa_Pareto_Select_YYYYMMDD_HHMMSS/

### Data files

- Raw merged data
- Processed dataset
- Initial and final selections
- Plotting datasets

### Figures

- Pareto selection plot
- Growth dynamics plot
- Parallel coordinates plot

---

## Parallel Coordinates Plot (Key Feature)

The PCP visualization includes:

- Background sampling (20%) to reduce clutter
- Full foreground display for selected individuals
- Group mean trend lines (dashed)
- Direct category labeling on plot

This allows intuitive interpretation of:

- Trait trade-offs
- Selection bias
- Pareto front structure

---

## Installation

```r
install.packages(c("dplyr","tidyr","ggplot2","stringr","readr",
                   "devEMF","viridis","emoa","VIM","ggrepel"))
```
---

## Usage
source("Alfalfa_Pareto_Selector_V3.4.R")

---

## License

MIT License

---
