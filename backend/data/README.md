# Backend Data Files

This directory contains the data files required for the AI Discovery module's molecular property prediction and clustering features.

## Data Files

### Required Files

1. **bp_data.csv** - Boiling Point data
   - Columns: SMILES, BoilingPoint
   - Used for predicting boiling points of molecules

2. **mp_data.csv** - Melting Point data
   - Columns: SMILES, MeltingPoint
   - Used for predicting melting points of molecules

3. **fp_data.csv** - Flash Point data
   - Columns: SMILES, FlashPoint
   - Used for predicting flash points of molecules

4. **cluster_data.csv** - Clustering data with QM9 properties
   - Columns: SMILES, alpha, mu, gap, homo, lumo, BP, FP, MP, cluster, ...
   - Used for molecular clustering and similarity analysis

## Generating Data Files

If the data files don't exist, they can be generated from the predict_block-master dataset:

```bash
cd /opt/molyte_web_v1.0/backend/data

python3 << 'EOF'
import pandas as pd

# Read the clustering data from predict_block-master
cluster_df = pd.read_csv('/opt/molyte_web_v1.0/predict_block-master/predict_block-master/dm/kmeans_clusters.csv')

# Create BP data
bp_df = cluster_df[['SMILES', 'BP']].dropna()
bp_df.columns = ['SMILES', 'BoilingPoint']
bp_df.to_csv('bp_data.csv', index=False)

# Create MP data
mp_df = cluster_df[['SMILES', 'MP']].dropna()
mp_df.columns = ['SMILES', 'MeltingPoint']
mp_df.to_csv('mp_data.csv', index=False)

# Create FP data
fp_df = cluster_df[['SMILES', 'FP']].dropna()
fp_df.columns = ['SMILES', 'FlashPoint']
fp_df.to_csv('fp_data.csv', index=False)

# Copy cluster data
cluster_df.to_csv('cluster_data.csv', index=False)

print("Data files created successfully!")
EOF
```

## Notes

- These CSV files are large (10-30 MB each) and are excluded from Git via .gitignore
- The backend service will automatically load these files on startup
- If files are missing, the service will attempt to load from predict_block-master as a fallback
- For production deployment, ensure these files are available in the backend/data directory

