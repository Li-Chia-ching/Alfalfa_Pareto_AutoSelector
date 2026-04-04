# Alfalfa Pareto Selector (V3.4)

## Overview

This pipeline implements a multi-trait selection framework for alfalfa (Medicago sativa) using Pareto optimization. It integrates early growth traits and a morphological score to identify elite individuals under a multi-objective decision space.

Key features include:

* Robust ingestion of matrix-formatted CSV data
* Adaptive missing data handling (deletion or KNN imputation)
* Trait standardization (Z-score)
* Pareto ranking using non-dominated sorting
* Diversity preservation via crowding distance
* Forced inclusion of predefined elite genotypes
* Publication-quality visualization (600 DPI)

---

## Input Data Format

The script expects three CSV files in matrix format:

* `Height_March.csv`
* `Height_April.csv`
* `Multifoliate_April.csv`

### Structure

* First row: column identifiers (e.g., A–U)
* First column: row identifiers (e.g., L1–L50)
* Remaining cells: numeric values

Each matrix is internally reshaped into long format and merged by a unique `Plant_ID` defined as:

```
Plant_ID = Row_ID-Column_ID
```

---

## Traits Used

* March Height (cm)
* April Height (cm)
* Multifoliate Score (1–5)

All traits are treated as **maximization objectives**.

---

## Workflow

### 1. Data Loading and Reshaping

* Reads matrix-style CSV files
* Converts to long format
* Merges all traits into a unified dataset

### 2. Missing Data Handling

* If missing rate < 5% → remove incomplete rows
* Otherwise → apply KNN imputation (k = 5)

### 3. Standardization

* Z-score normalization applied to all traits

### 4. Pareto Optimization

* Non-dominated sorting (`emoa::nds_rank`)
* Rank 1 = Pareto front

### 5. Diversity Preservation

* Crowding distance calculated within each Pareto rank
* Composite score:

```
Composite Score = Sum(Z-scores) × (1 + Crowding Distance)
```

### 6. Initial Selection

* Top 50 individuals selected based on:

  * Pareto rank (ascending)
  * Composite score (descending)

### 7. Forced Inclusion Mechanism

* Predefined list of elite genotypes is always retained if present
* Remaining slots are filled using algorithmic ranking

---

## Outputs

A timestamped directory is generated:

```
Alfalfa_Pareto_Select_YYYYMMDD_HHMMSS/
```

### Data Files

* `00_Raw_Merged_Data.csv` — merged raw dataset
* `01_Full_Processed_Data.csv` — standardized and ranked data
* `02_Initial_Top50.csv` — initial selection (pure algorithm)
* `03_Final_Selected_50.csv` — final selection (after forced inclusion)
* `04_Plotting_Data.csv` — visualization dataset
* `05_Parallel_Coords_Data.csv` — parallel coordinates data

### Figures (600 DPI)

* `Plot_Pareto_Selection.pdf`
* `Plot_Growth_Dynamics.pdf`
* `Plot_Parallel_Coords.pdf`

---

## Visualization

### 1. Pareto Selection Plot

* X-axis: April height
* Y-axis: Multifoliate score
* Color: Pareto rank
* Size: Composite score
* Highlights final selected individuals

### 2. Growth Dynamics Plot

* March vs April height
* Includes 1:1 reference line

### 3. Parallel Coordinates Plot

* Displays multi-trait trade-offs
* Highlights Pareto front and selected individuals

---

## Dependencies

Required R packages:

```
dplyr
tidyr
ggplot2
stringr
readr
devEMF
viridis
emoa
VIM
ggrepel
```

The script automatically installs missing packages from CRAN.

---

## Reproducibility

* Random seed fixed (`set.seed(2025)`) for KNN imputation
* All intermediate datasets are exported

---

## Notes and Limitations

* Assumes all traits are positively associated with selection targets
* KNN imputation is sensitive to trait scaling
* Crowding distance is an approximation of NSGA-II
* Forced inclusion introduces domain-driven bias into the selection process

---

## Intended Use

This tool is designed for:

* Early-generation selection in breeding programs
* Multi-trait decision support
* Visualization of trade-offs among agronomic traits

It is best suited for structured experimental datasets with consistent layout.

---

## Version

**V3.4 (March–April 2025)**

Includes enhanced visualization, improved robustness, and a hybrid selection mechanism combining optimization and expert knowledge.
