# Alfalfa Pareto Selector V3.5.1.

## Overview

Alfalfa Pareto Selector V3.5.1. is an R-based pipeline for multi-trait selection in alfalfa (*Medicago sativa*), integrating data preprocessing, multi-objective optimization, and publication-quality visualization.

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

# Alfalfa Phenotypic Structure & Diversity Analysis (Pre-Pareto Module)

## Overview

This repository provides an R-based pipeline for analyzing phenotypic population structure in alfalfa (*Medicago sativa*). The workflow is designed to identify subpopulation patterns in plant architecture traits prior to downstream selection procedures such as GWAS sampling or Pareto optimization.

The pipeline integrates robust data reshaping, missing data handling, multivariate analysis, clustering, and publication-quality visualization.

---

## Key Features

* Robust matrix-to-long data transformation (`matrix_to_long`)
* Automatic handling of irregular CSV matrix formats
* Adaptive missing data strategy:

  * <5% missing → row removal
  * ≥5% missing → KNN imputation
* Z-score standardization with matrix-safe conversion
* Phenotypic correlation analysis
* Principal Component Analysis (PCA)
* K-means clustering with elbow method evaluation
* High-resolution (600 DPI) figure export
* Fully reproducible output structure

---

## Input Data Format

The pipeline expects **matrix-style CSV files** (not tidy format):

### Required Files

* `Height_March.csv`
* `Height_April.csv`
* `Multifoliate_April.csv`

### Expected Structure

|    | C1 | C2 | C3 |
| -- | -- | -- | -- |
| R1 |    |    |    |
| R2 |    |    |    |

* First row: column identifiers
* First column: row identifiers (family or line)
* Remaining cells: trait values

### Output ID Construction

Each observation is converted into:

```
Plant_ID = Row_ID - Col_ID
```

---

## Workflow

### 1. Data Loading

* Convert matrix-format CSV to long format
* Merge multiple traits by `Plant_ID`

### 2. Missing Data Handling

* Compute missing rate per trait
* Apply rule-based cleaning or KNN imputation

### 3. Standardization

* Z-score normalization using `scale()`
* Forced numeric conversion to avoid matrix-column issues

### 4. Correlation Analysis

* Pearson correlation matrix
* Visualization via `corrplot`

### 5. PCA

* Performed on standardized traits
* Biplot visualization using `factoextra`

### 6. Clustering

* Elbow method for optimal K estimation
* K-means clustering (`nstart = 25`)

### 7. Visualization

* PCA-based cluster visualization
* Elliptical grouping (with tolerance for outliers)

### 8. Export

* Full annotated dataset
* PCA variance explained
* Cluster summary statistics

---

## Output Structure

Results are saved in:

```
Pheno_Structure_Analysis_YYYYMMDD/
```

### Data Files

* `Data_01_Pheno_Structure_Full.csv`
* `Data_02_PCA_Variance.csv`
* `Data_03_Cluster_Profiles.csv`

### Figures

* `Plot_00_Optimal_K_Elbow.pdf`
* `Plot_01_Correlation_Matrix.pdf`
* `Plot_02_PCA_Biplot.pdf`
* `Plot_03_Subpopulations.pdf`

---

## Dependencies

```r
dplyr
tidyr
ggplot2
readr
factoextra
corrplot
cluster
stringr
VIM
```

Install with:

```r
install.packages(c("dplyr","tidyr","ggplot2","readr","factoextra","corrplot","cluster","stringr","VIM"))
```

---

## Usage

```r
source("alfalfa_pheno_structure_pca_kmeans_pipeline.R")
```

---

## Methodological Notes

* PCA is performed on pre-scaled data (no internal centering/scaling)
* KNN imputation assumes phenotypic similarity reflects biological proximity
* Ellipse fitting may fail for extreme outliers (does not affect clustering results)
* Matrix-to-long conversion ensures compatibility with field-style data layouts

---

## Applications

* GWAS core population selection
* Phenotypic diversity assessment
* Pre-selection for sequencing
* Trait architecture exploration

---

## License

MIT License

---
