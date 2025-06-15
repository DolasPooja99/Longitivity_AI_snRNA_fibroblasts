
#  snRNA-seq ML Data Preprocessing for Human Optic Nerve and Optic Nerve Head fibroblasts  Cells

This repository includes a Python script `processing.py` that processes `.h5ad` single-nucleus RNA sequencing (snRNA-seq) data and extracts essential features for machine learning (ML) tasks. It generates structured `.parquet` files, including expression matrices, metadata, and principal components (PCA).

Dataset reference:  
[View on Hugging Face](https://huggingface.co/datasets/longevity-db/human-optic-nerve-fibroblasts-snRNAseq)

---

## Dataset Summary

This dataset contains processed snRNA-seq data from **human optic nerve** and **optic nerve head** fibroblasts cells, intended for downstream ML modeling and analysis. The original data is stored in `.h5ad` format using the `AnnData` structure, which is standard in single-cell bioinformatics.

---

## What the Script Does (`processing.py`)

The script processes the `.h5ad` file and outputs the following files in `.parquet` format for efficient loading in Python, R, or big data environments:

### 1. **Creates an Output Directory**
- Based on dataset name (e.g., `snRNA-seq_of_human_optic_nerve_and_optic_nerve_head_fibroblasts_cells_ml_data`)

### 2. **Loads `.h5ad` File**
- Uses `anndata` to read the `.h5ad` file into memory.

### 3. **Exports Key Data Components**
- `expression.parquet`: The expression matrix (genes × cells)
- `feature_metadata.parquet`: Metadata about genes (e.g., gene symbols, IDs)
- `observation_metadata.parquet`: Metadata about cells (e.g., donor ID, cluster, tissue type)

### 4. **Performs PCA**
- Optional standard scaling (`StandardScaler`) applied to gene expression data.
- Reduces dimensionality to `N_PCA_COMPONENTS` (default: 50).
- Outputs:
  - `pca_components.parquet`: Principal component scores for each cell.
  - `pca_explained_variance.parquet`: Variance explained by each principal component.

---

## Requirements

Install required packages:

```bash
pip install pandas anndata scikit-learn pyarrow
```

---

## Usage

1. Replace the path in the script:

```python
H5AD_FILE_PATH = "/path/to/your/file.h5ad"
```

2. Run the script:

```bash
python processing.py
```

3. Output directory structure:

```
snRNA-seq_of_human_optic_nerve_and_optic_nerve_head_fibroblasts_cells_ml_data/
│
├── expression.parquet
├── feature_metadata.parquet
├── observation_metadata.parquet
├── pca_components.parquet
└── pca_explained_variance.parquet
```

---

## Notes

- PCA helps reduce dimensionality and is useful for downstream tasks like clustering, visualization, or model input.
- All files are saved in `.parquet` format for speed and cross-platform compatibility.
- The number of PCA components and whether to apply scaling can be adjusted in the script via:
  ```python
  N_PCA_COMPONENTS = 50
  APPLY_SCALING_BEFORE_PCA = True
  ```

---

## References
- Source Dataset:

## Final Output:
- Final Dataset :[longevity-db/snRNA-seq_of_human_optic_nerve_and_optic_nerve_head_fibroblasts_cells_ml_data](https://huggingface.co/datasets/longevity-db/human-optic-nerve-fibroblasts-snRNAseq)
- File Format: [AnnData Documentation](https://anndata.readthedocs.io/en/latest/)
- PCA: [Scikit-learn PCA](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html)


Contributed by CellVPA Team
Venkatachalam, Pooja, Albert

