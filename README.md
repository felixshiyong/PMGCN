# Identifying Omic Biomarkers for Chronic Inflammatory Diseases Using Percolation on Multi-disease Gene Co-expression Networks

## Project Description
PMGCN is a computational biology pipeline designed to identify robust and generalizable biomarker candidates for chronic inflammatory diseases (e.g., Ulcerative Colitis, Crohn's disease). The core methodology involves integrating bulk transcriptomic data from multiple related diseases to construct a unified gene co-expression network. A key network percolation analysis is then applied to this network to pinpoint genes that are crucial for network connectivity, representing potential diagnostic or therapeutic targets.

## Repository Structure & Workflow
The analysis is modular. To reproduce the study, execute the core scripts in the following order:

| File | Language | Description |
| :--- | :--- | :--- |
| **[`1_bulk_process.R`]** | R | **Data Foundation:** Downloads and processes bulk RNA-Seq data (e.g., from GEO). Performs normalization and generates the cross-disease gene co-expression network. |
| **[`2_optimal_P_UC_optimized.ipynb`]** | Python | **Core Analysis:** Executes the controlled percolation process on the multi-disease network to identify key genes critical for connectivity. |
| **[`3_Data Preparation for Machine Learning.R`]** | R | **Validation Prep:** Prepares the feature expression matrix for downstream machine learning analysis (e.g., classifier training) based on percolation results, corresponding to Figures 6 & 7 in the associated manuscript. |
| **[`4_sc_process.R`]** | R | **Deep Validation:** Fetches and processes relevant single-cell RNA-Seq data to validate candidate biomarker expression patterns at cellular resolution. |

**Other Files:** The functions of other files for figure generations are indicated by their respective filenames.

## Getting Started
### Prerequisites
Ensure you have the following environments configured:

*   **R Environment** (version >= 4.4.0 recommended):
    ```r
    # Install core R packages
    install.packages(c("tidyverse", "WGCNA", "BiocManager"))
    BiocManager::install(c("limma", "DESeq2", "Seurat"))
    ```
*   **Python Environment** (version >= 3.10 recommended):
    ```bash
    pip install pandas numpy networkx matplotlib scikit-learn jupyter
    ```

### Execution
Clone the repository and run the scripts sequentially as outlined in the workflow table above.

## Data Sources
All transcriptomic data analyzed in this project are from public repositories:
*   **Bulk RNA-Seq Data:** Primary sources include GEO.
*   **Single-cell RNA-Seq Data:** Sourced from GEO.
*Specific dataset accession numbers (e.g., GSE IDs) are detailed within the corresponding processing scripts or the manuscript's methods section.*

## Citation


## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
