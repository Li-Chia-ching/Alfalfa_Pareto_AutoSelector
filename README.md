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

# Alfalfa Agronomic Population Structure & Diversity Analysis 🌿
**Pre-Pareto Optimization Module (V2.0)**

![R](https://img.shields.io/badge/R-%23276DC3.svg?style=flat&logo=r&logoColor=white)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-Publication_Ready-brightgreen)

## 📌 Overview
This repository contains the advanced R pipeline for analyzing the **agronomic population structure and plant architecture diversity** of Alfalfa (*Medicago sativa L.*). 

Designed as a rigorous pre-requisite module before Multi-Objective Pareto Selection, this script evaluates the phenotypic landscape of 1,050 individual plants (50 families × 21 replicates). It identifies natural subpopulations, assesses trait redundancies, and ensures that subsequent elite plant selection maintains robust genetic diversity.

## ✨ Key Features & Recent Updates

* 🛡️ **Strict Master Grid Integration (Anti-Phantom-Data):** Implements a definitive $1050$ plant matrix (`L1-L50` × `A-U`). It acts as a bulletproof parser that inherently filters out Excel-generated "ghost" rows/columns (e.g., `NA.xxx`) during CSV loading.
* 📈 **Robust Statistical Upgrades:** Shifted from Pearson to **Spearman rank correlation** to scientifically handle the mixture of continuous variables (Plant Height) and ordinal variables (Multifoliate Score).
* 🎓 **Publication-Ready Terminology:** Upgraded visual and output terminology from generic "Phenotypic Traits" to specific **"Agronomic Traits"** and **"Plant Architecture"**, adhering to high-impact agronomy journal standards.
* 🎨 **Colorblind-Friendly Aesthetics:** Fully integrated with the `viridis` color palette (e.g., *cividis*, *magma*) for all K-means clustering and correlation plots, ensuring maximum accessibility and 600 DPI publication quality.
* 🤖 **Intelligent Imputation:** Automated KNN (k=5) imputation triggered only if the missing data threshold exceeds 5%.

## 📊 Analytical Pipeline

1.  **Data Ingestion & Cleaning:** Reads raw field-layout CSVs and enforces the Master Grid structure.
2.  **Missing Value Handling:** K-Nearest Neighbors (KNN) imputation via the `VIM` package.
3.  **Z-Score Standardization:** Normalizes trait variances.
4.  **Agronomic Correlation (Spearman):** Assesses trade-offs (e.g., Biomass potential vs. Forage quality).
5.  **Principal Component Analysis (PCA):** Reduces dimensionality and visualizes overall plant architecture distribution.
6.  **K-means Clustering:** Evaluates optimal $K$ (Elbow method) and stratifies the population into distinct architecture sub-groups to monitor genetic diversity.

## 📂 Input Data Format

The script expects three CSV files in a spatial field layout (Matrix format) in the working directory:
* `Height_March.csv` (Continuous)
* `Height_April.csv` (Continuous)
* `Multifoliate_April.csv` (Ordinal 1-5 score)

**CSV Structure Example:**
|      | A   | B   | C   | ... | U   |
|------|-----|-----|-----|-----|-----|
| **L1** | 45  | 42  | NA  | ... | 50  |
| **L2** | 39  | 40  | 41  | ... | 38  |

## 🚀 Usage

Simply run the R script in your IDE (e.g., RStudio). The script will automatically check for missing dependencies and install them.

```R
# Required packages automatically handled by the script:
# dplyr, tidyr, ggplot2, readr, factoextra, corrplot, cluster, stringr, VIM, viridis
source("Alfalfa_Agronomic_Structure.R")
