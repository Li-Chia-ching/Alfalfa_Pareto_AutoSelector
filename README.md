# Alfalfa Pareto Selector V3.5

## Overview

**Alfalfa Pareto Selector V3.5** is an integrated R pipeline for multi-trait selection in alfalfa (Medicago sativa), designed for field phenotyping datasets in matrix format. The workflow combines data cleaning, imputation, normalization, multi-objective optimization (Pareto ranking), and publication-quality visualization.

This version introduces a robust parallel coordinates plot with explicit category labeling and improved graphical output suitable for scientific publication.

---

## Key Features

* Robust CSV parser for matrix-format field data (e.g., L1–L50 × A–U)
* Automated missing data handling:

  * Row removal (<5% missing)
  * KNN imputation (k = 5) for larger gaps
* Z-score standardization across traits
* Multi-objective optimization:

  * Pareto ranking (non-dominated sorting)
  * Crowding distance estimation
* Composite scoring for prioritization
* Forced inclusion mechanism for elite individuals
* Publication-ready visualization:

  * Pareto scatter plot
  * Growth dynamics plot
  * Parallel coordinates plot with category labels
* 600 DPI high-resolution output

---

## Input Data Format

Input files must be CSV matrices with the following structure:

* First row: column identifiers (e.g., A–U)
* First column: row identifiers (e.g., L1–L50)
* Remaining cells: numeric trait values

### Required Files

| File Name                | Trait Description        |
| ------------------------ | ------------------------ |
| `Height_March.csv`       | Plant height in March    |
| `Height_April.csv`       | Plant height in April    |
| `Multifoliate_April.csv` | Multifoliate score (1–5) |

---

## Workflow

1. **Data Import & Reshaping**

   * Convert matrix-format data into long format
   * Generate unique `Plant_ID`

2. **Missing Data Handling**

   * Remove rows or apply KNN imputation

3. **Standardization**

   * Z-score normalization per trait

4. **Pareto Optimization**

   * Non-dominated sorting (Pareto rank)
   * Crowding distance calculation

5. **Selection Strategy**

   * Top 50 individuals selected
   * Forced inclusion of user-defined elite plants

6. **Visualization**

   * Multi-panel plotting for trait relationships and trade-offs

---

## Output Files

All outputs are saved in a timestamped directory:

```
Alfalfa_Pareto_Select_YYYYMMDD_HHMMSS/
```

### Data Outputs

| File                          | Description                        |
| ----------------------------- | ---------------------------------- |
| `00_Raw_Merged_Data.csv`      | Merged raw dataset                 |
| `01_Full_Processed_Data.csv`  | Processed dataset with scores      |
| `02_Initial_Top50.csv`        | Initial top 50 selection           |
| `03_Final_Selected_50.csv`    | Final selection after constraints  |
| `04_Plotting_Data.csv`        | Data used for visualization        |
| `05_Parallel_Coords_Data.csv` | Data for parallel coordinates plot |

### Figures

| File                        | Description            |
| --------------------------- | ---------------------- |
| `Plot_Pareto_Selection.pdf` | Pareto scatter plot    |
| `Plot_Growth_Dynamics.pdf`  | March vs April growth  |
| `Plot_Parallel_Coords.pdf`  | Multi-trait trade-offs |

---

## Installation

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "stringr", "readr",
                   "devEMF", "viridis", "emoa", "VIM", "ggrepel"))
```

---

## Usage

1. Place input CSV files in the working directory
2. Run the script:

```r
source("Alfalfa_Pareto_Selector_V3.4.R")
```

3. Outputs will be generated automatically

---

## Methodological Notes

* **Pareto ranking** identifies non-dominated individuals across traits
* **Crowding distance** ensures diversity within Pareto fronts
* **Composite score** integrates performance and diversity
* **Parallel coordinates plot** enables intuitive visualization of multi-trait trade-offs

---

## License

 **MIT License** 

---
