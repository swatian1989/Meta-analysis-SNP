# SNP Meta-Analysis of Egyptian and Sham Cohort Genetic Data

## Project Description

### Egyptian Cohort

This project implements a complete, transparent, and reproducible meta-analysis and clustering workflow for SNP (Single Nucleotide Polymorphism) genetic association studies in Egyptian populations. Using robust statistics and unsupervised machine learning, it enables identification of the most significant genetic variants, evaluation of result consistency, and discovery of hidden patterns in the data. All outputs are publication-ready and support open science.

#### Key Steps:

- Data import and cleaning
- Funnel plot (publication bias)
- Heterogeneity analysis (Cochran's Q, I²)
- SNP ranking (by OR and p-value)
- Clustering (KMeans) to find SNP groups
- Meta-regression
- Automatic export to Excel for downstream analysis

**Significance:**

- Highlights top genetic risk loci in Egyptian cohorts
- Reveals heterogeneity and sample-size effects
- Supports research in population genetics and precision medicine
- Open, reproducible pipeline for further studies

---

### Sham Data Cohort

The same workflow is applied to the `sham_data1.xlsx` file, enabling comparative or secondary analysis for another population or cohort ("Sham"). The code is fully reusable—simply update the input file and run.

#### Key Steps:

- Load and clean `sham_data1.xlsx`
- Robust numerical conversion for all relevant columns
- All downstream analyses identical to Egyptian cohort pipeline

**Significance:**

- Enables multi-cohort or multi-population analysis with no extra code changes
- Directly comparable outputs for all cohorts

---

## Python Script (Universal Version for Egypt or Sham Data)

Paste this into `meta_analysis_project.py`. Change only the Excel filename for your cohort:

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import statsmodels.api as sm

# Change this to your Excel file for Egypt or Sham analysis:
FILENAME = 'sham_data1.xlsx'  # For sham cohort; use 'egypt_data.xlsx' for Egypt

df = pd.read_excel(FILENAME)

# These column names work for both Egypt and Sham files if headers are consistent
# If not, run: print(df.columns)
to_numeric_cols = ['OR', 'P-Value', 'No. of Patients', 'No. of controls', 'CI_low', 'CI_high']
for col in to_numeric_cols:
    df[col] = pd.to_numeric(df[col], errors='coerce')

df['logOR'] = np.log(df['OR'])
df['SE_logOR'] = (np.log(df['CI_high']) - np.log(df['CI_low'])) / (2 * 1.96)
plt.figure(figsize=(7,9))
plt.scatter(df['logOR'], 1/df['SE_logOR'], alpha=0.7)
plt.xlabel('Log Odds Ratio')
plt.ylabel('1 / SE of Log(OR)')
plt.title('Funnel Plot')
plt.show()

df['var_logOR'] = df['SE_logOR'] ** 2
df['weight'] = 1 / df['var_logOR']
mean_logOR = np.sum(df['logOR'] * df['weight']) / np.sum(df['weight'])
Q = np.sum(df['weight'] * (df['logOR'] - mean_logOR) ** 2)
df_deg = len(df) - 1
I2 = max(0, 100 * (Q - df_deg) / Q) if Q > df_deg else 0
print(f"Cochran's Q: {Q:.2f}")
print(f"I²: {I2:.1f}%")

top_by_or = df.sort_values('OR', ascending=False)[['Gene', 'SNP', 'OR', 'P-Value']].head(10)
top_by_p = df.sort_values('P-Value')[['Gene', 'SNP', 'OR', 'P-Value']].head(10)
print("\nTop 10 SNPs by OR:\n", top_by_or)
print("\nTop 10 Most Significant SNPs:\n", top_by_p)

features = df[['OR', 'P-Value', 'No. of Patients', 'No. of controls']].copy()
features['P-Value'] = -np.log10(features['P-Value'].replace(0, 1e-12))
scaler = StandardScaler()
features_scaled = scaler.fit_transform(features.fillna(0))
kmeans = KMeans(n_clusters=3, random_state=0)
df['Cluster'] = kmeans.fit_predict(features_scaled)
print("\nCluster assignments (first 10):\n", df[['Gene', 'SNP', 'Cluster']].head(10))

df['Significant (p<0.05)'] = df['P-Value'] < 0.05
summary = df[['Gene', 'SNP', 'Type/position of SNP', 'No. of Patients',
              'No. of controls', 'P-Value', 'OR', 'CI_low', 'CI_high', 'Allele', 'Significant (p<0.05)', 'Cluster']]
print("\nSummary Table (first 10):\n", summary.head(10))

X = df[['No. of Patients', 'No. of controls']].astype(float)
X = sm.add_constant(X)
y = df['OR'].astype(float)
meta_reg = sm.OLS(y, X).fit()
print("\nMeta-regression results:\n", meta_reg.summary())

df.to_excel(FILENAME.replace('.xlsx', '_with_clusters.xlsx'), index=False)
```

---

## How to Use

1. Copy this script to a `.py` file
2. Change the `FILENAME` at the top to your data (Egypt or Sham)
3. Run in your Python environment
4. Find new Excel with results and clusters in the same folder

---

## Research and Clinical Impact

- Accelerates gene discovery in Egyptian and related populations
- Enables robust, side-by-side comparison of different cohorts
- Open code for adaptation and reproducibility

