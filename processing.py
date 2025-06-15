import pandas as pd
import anndata as ad
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler # For optional scaling before PCA

# --- Configuration ---
# IMPORTANT: Replace with the actual path to your .h5ad file
H5AD_FILE_PATH = "/path/to/your/file.h5ad"

# Name of the directory to store the output Parquet files
OUTPUT_DIR_NAME = "snRNA-seq_of_human_optic_nerve_and_optic_nerve_head_endothelial_cells_ml_data"

# PCA Configuration
N_PCA_COMPONENTS = 50 # Number of principal components to compute. Adjust as needed.
APPLY_SCALING_BEFORE_PCA = True # Set to True to scale data (mean=0, std=1) before PCA

# --- 1. Create Output Directory ---
os.makedirs(OUTPUT_DIR_NAME, exist_ok=True)
print(f"Created output directory: {OUTPUT_DIR_NAME}")

# --- 2. Load the H5AD file into an AnnData object ---
try:
    adata = ad.read_h5ad(H5AD_FILE_PATH)
    print(f"Successfully loaded AnnData object from: {H5AD_FILE_PATH}")
    print(f"AnnData object shape: {adata.shape}")
except FileNotFoundError:
    print(f"Error: H5AD file not found at {H5AD_FILE_PATH}. Please check the path.")
    exit()
except Exception as e:
    print(f"An error occurred while loading the H5AD file: {e}")
    exit()

# --- 3. Save adata.X as expression.parquet ---
expression_parquet_path = os.path.join(OUTPUT_DIR_NAME, "expression.parquet")

# Convert adata.X (which can be a NumPy array or sparse matrix) to a pandas DataFrame
# Ensure it has observation names (index) and variable names (columns)
if isinstance(adata.X, (pd.DataFrame, pd.Series)):
    df_expression = adata.X
elif hasattr(adata.X, 'toarray'): # If it's a sparse matrix (e.g., scipy.sparse.csr_matrix)
    df_expression = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
else: # If it's a dense NumPy array
    df_expression = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)

df_expression.to_parquet(expression_parquet_path, index=False)
print(f"Saved expression data to: {expression_parquet_path}")

# --- 4. Save adata.var as feature_metadata.parquet ---
feature_metadata_parquet_path = os.path.join(OUTPUT_DIR_NAME, "feature_metadata.parquet")
adata.var.to_parquet(feature_metadata_parquet_path, index=False)
print(f"Saved feature metadata to: {feature_metadata_parquet_path}")

# --- 5. Save adata.obs as observation_metadata.parquet ---
# This often contains cell type annotations, patient info, etc., crucial for ML tasks.
observation_metadata_parquet_path = os.path.join(OUTPUT_DIR_NAME, "observation_metadata.parquet")
adata.obs.to_parquet(observation_metadata_parquet_path, index=False)
print(f"Saved observation metadata to: {observation_metadata_parquet_path}")

# --- 6. Perform PCA and save results ---
print(f"\nStarting PCA with {N_PCA_COMPONENTS} components...")

# Get the expression data for PCA. It's best to work with a dense array for PCA.
# If adata.X is a sparse matrix, convert it to a dense array.
X_data = df_expression.values # Get NumPy array from DataFrame

if APPLY_SCALING_BEFORE_PCA:
    print("Scaling expression data before PCA...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_data)
else:
    X_scaled = X_data

pca = PCA(n_components=N_PCA_COMPONENTS)
pca_transformed_data = pca.fit_transform(X_scaled)

# Create a DataFrame for the PCA transformed data
pca_columns = [f"PC{i+1}" for i in range(N_PCA_COMPONENTS)]
df_pca = pd.DataFrame(pca_transformed_data, index=adata.obs_names, columns=pca_columns)

pca_parquet_path = os.path.join(OUTPUT_DIR_NAME, "pca_components.parquet")
df_pca.to_parquet(pca_parquet_path, index=False)
print(f"Saved PCA transformed data to: {pca_parquet_path}")

# Save the explained variance ratio
df_explained_variance = pd.DataFrame({
    'PrincipalComponent': [f"PC{i+1}" for i in range(len(pca.explained_variance_ratio_))],
    'ExplainedVarianceRatio': pca.explained_variance_ratio_
})
explained_variance_parquet_path = os.path.join(OUTPUT_DIR_NAME, "pca_explained_variance.parquet")
df_explained_variance.to_parquet(explained_variance_parquet_path, index=False)
print(f"Saved PCA explained variance ratio to: {explained_variance_parquet_path}")

print(f"\nAll required Parquet files (including PCA) have been created in the '{OUTPUT_DIR_NAME}' directory.")
print("You can now use these files for your submission.")
