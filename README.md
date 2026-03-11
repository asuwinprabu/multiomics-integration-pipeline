# Multi-Omics Biomarker Discovery for Immunotherapy Response Prediction

## 🧬 Project Overview
Predicting which patients will respond to Immune Checkpoint Blockade (ICB) therapy is one of the most pressing challenges in precision oncology. This project builds an end-to-end Machine Learning pipeline in **R** to address this clinical problem. 

By integrating multi-omics data from The Cancer Genome Atlas (TCGA), including **Tumor Mutational Burden (TMB)** and **Immune Gene Expression Signatures** (e.g., IFN-gamma), we train diagnostic classifiers (Random Forest via `tidymodels`) to predict patient response. The culmination of this project is an interactive **Shiny Dashboard** designed for clinical interpretability.

---

## 🛠️ Technology Stack & Skills Demonstrated
This repository demonstrates an industry-standard Data Science & Bioinformatics workflow optimized for Pharma and Biotech applications:
* **Complex Data Sourcing**: `TCGAbiolinks` to programmatically query and wrangle public clinical and omics cohorts.
* **Biologically-Informed Feature Engineering**: Calculating TMB from MAF files and extracting established RNA-seq signatures.
* **Rigorous Machine Learning**: `tidymodels` framework implemented with strict cross-validation to prevent data leakage.
* **Clinical Interpretability**: `Shiny`, `bslib`, and `vip` (Variable Importance) to deliver complex model results to non-technical stakeholders via an interactive UI.

---

## 📂 Repository Structure
```
├── R/
│   ├── 01_fetch_data.R           # Queries TCGA-SKCM (Melanoma) Clinical & RNA-seq data
│   ├── 02_feature_engineering.R  # Extracts TMB and IFN-gamma signature scores
│   └── 03_train_models.R         # Tidymodels RF classification & hyperparameter tuning
├── app/
│   ├── app.R                     # Interactive Shiny Clinician Dashboard
│   └── models/                   # Contains the trained .rds workflow artifacts
├── data/                         # Ignored by git; raw GDC downloads and processed matrices
├── renv.lock                     # R environment dependency specifications
└── README.md
```

---

## 🚀 How to Run the Pipeline

### 1. Environment Setup
This project uses `renv` to guarantee reproducibility. After cloning the repository, open the project in RStudio and run:
```R
renv::restore()
```

### 2. Execute the Pipeline
Run the scripts sequentially. *(Note: `01_fetch_data.R` downloads real TCGA data which may take some time depending on bandwidth).*

```bash
# Sourcing the data
Rscript R/01_fetch_data.R

# Building the feature matrix
Rscript R/02_feature_engineering.R

# Training the ML Classifier
Rscript R/03_train_models.R
```

### 3. Launch the Dashboard
Once the Random Forest model is saved to `app/models/rf_model.rds`, launch the Shiny proxy:
```R
shiny::runApp("app")
```

---

## 📊 Dashboard Preview
*(Add a screenshot of your Shiny App here before pushing to GitHub)*

The dashboard allow users to:
1. Input synthetic or real patient biomarker profiles.
2. View the probability of response to Immunotherapy (via a Plotly gauge).
3. Inspect the underlying Machine Learning model via Global SHAP/Gini Impurity feature importance plots.
