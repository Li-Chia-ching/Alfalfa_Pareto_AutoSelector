---

# 📘 Alfalfa Pareto Selector V3.6: **Robust Multi-Trait Selection Pipeline with Manual Elite Preservation**

---

## 🔬 Overview

This repository provides a **robust, publication-ready R pipeline** for multi-trait selection in alfalfa populations based on:

* March plant height
* April plant height
* Multifoliate score (1–5)

The pipeline integrates:

* **Strict experimental design enforcement (1050 plants)**
* **Robust CSV parsing with artifact filtering**
* **KNN imputation for missing data**
* **Pareto front ranking (multi-objective optimization)**
* **Crowding distance for diversity preservation**
* **Manual elite plant injection**
* **Publication-quality visualization (600 DPI)**

---

## ⚙️ Key Features (V3.6 Update)

### ✅ 1. Strict Master Grid Enforcement

* A **fixed 50 × 21 design (L1–L50 × A–U)** is enforced
* Guarantees **exactly 1050 plants**, regardless of CSV integrity
* Prevents:

  * Missing rows (e.g., L49 deletion)
  * Structural inconsistency across datasets

👉 Even if raw files are corrupted, output remains **structurally valid and analysis-ready**

---

### ✅ 2. Robust CSV Parsing & “Ghost Data” Removal

* Input files are read as **raw character matrices**
* Automatically removes:

  * Empty rows (`NA`)
  * Non-experimental rows (e.g., `"平均值"`, `"Mean"`)
  * Excel artifacts (blank / NULL / malformed values)

✔️ Equivalent to:

```r
filter(!is.na(Row_ID), Row_ID %in% valid_rows)
```

👉 Ensures only **true experimental observations** enter the pipeline

---

### ✅ 3. Safe Left Join onto Master Grid

* Uses **master_grid + left_join()** strategy

**Result:**

* Missing observations → converted to **true NA**
* Structural completeness maintained

👉 Critical advantage:

> Even if an entire row is missing in CSV, all 21 plants are preserved with NA values

---

### ✅ 4. Adaptive Missing Data Strategy

* If missing rate < 5% → **row removal**
* Else → **KNN imputation (k = 5)**

✔️ Guarantees:

```r
stopifnot(no NA remains)
```

---

### ✅ 5. Pareto-Based Multi-Trait Selection

* Traits standardized via **Z-score**
* Non-dominated sorting using:

  * `emoa::nds_rank()`
* Diversity maintained via:

  * **Crowding distance (fixed: no Inf values)**

👉 Final ranking combines:

```
Composite Score = Trait Sum × (1 + Crowding Distance)
```

---

### ✅ 6. ⭐ Manual Elite Preservation 

This pipeline **explicitly preserves manually pre-selected elite plants**.

#### How it works:

* User-defined list (e.g., `"A10", "B30", ...`) is converted to `Plant_ID`
* Automatically matched to dataset
* Injected into final selection

#### Final selection logic:

```
Final 50 = Forced elites + Algorithm-selected fillers
```

#### Guarantees:

* Manual selections are **never lost**
* Algorithm fills remaining slots optimally

👉 This is critical for:

* Field knowledge integration
* Breeding experience retention
* Grant justification & reproducibility

---

## 📊 Visualization (Publication-Ready)

All figures are exported at **600 DPI**, optimized for journals.

### 1. Pareto Front Scatter

* April Height vs. Multifoliate Score
* Color = Pareto rank
* Size = Composite score
* Selected plants highlighted with **red circles**

---

### 2. Growth Dynamics Plot

* March vs. April height
* Diagonal reference line (growth consistency)
* Colorblind-friendly **viridis (cividis)** palette

---

### 3. Parallel Coordinates Plot (NEW)

* Multi-trait trade-offs visualization
* Groups:

  * Rank 1 front
  * Other selected
  * Background population
* Includes:

  * Mean trend lines
  * Category labels directly on plot
  * Clutter reduction via background sampling

---

## 📁 Output Files

All results are saved in a timestamped directory:

```
Alfalfa_Pareto_Select_YYYYMMDD_HHMMSS/
```

| File                          | Description                             |
| ----------------------------- | --------------------------------------- |
| `00_Raw_Merged_Data.csv`      | Master grid merged raw data (with NA)   |
| `01_Full_Processed_Data.csv`  | Cleaned + standardized + ranked dataset |
| `02_Final_Selected_50.csv`    | Final selected plants (核心结果)            |
| `03_Plotting_Data.csv`        | Data used for visualization             |
| `04_Parallel_Coords_Data.csv` | Long-format data for PCP plot           |

---

## 📥 Input Requirements

Three CSV files (matrix format):

* `Height_March.csv`
* `Height_April.csv`
* `Multifoliate_April.csv`

### Format:

|    | A   | B   | C   | ... |
| -- | --- | --- | --- | --- |
| L1 | ... | ... | ... |     |
| L2 | ... | ... | ... |     |

⚠️ Notes:

* First row = column IDs (A–U)
* First column = row IDs (L1–L50)
* Extra rows (e.g., averages) are **automatically removed**

---

## 🧪 Reproducibility

* Random seed fixed:

```r
set.seed(2025)
```

* Deterministic selection pipeline
* Fully reproducible outputs

---

## 🚀 Typical Use Case

This pipeline is designed for:

* GWAS pre-selection
* Core germplasm extraction
* Multi-trait breeding decisions
* Grant proposal justification

---

## 📌 Conceptual Strength

Compared to naive selection approaches, this pipeline:

| Aspect                | Traditional       | This Pipeline                 |
| --------------------- | ----------------- | ----------------------------- |
| Missing data          | Breaks structure  | Reconstructed via grid        |
| CSV robustness        | Fragile           | Artifact-resistant            |
| Multi-trait selection | Weighted sum bias | Pareto optimal                |
| Diversity             | Ignored           | Preserved (crowding distance) |
| Expert knowledge      | Lost              | Explicitly integrated         |

---

### Major updates:

* Strict 1050 master grid enforcement
* Ghost data filtering
* Safe merging strategy
* Stable crowding distance (no Inf/NaN)
* Parallel coordinates visualization
* Colorblind-friendly palettes
* Manual elite preservation mechanism

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
