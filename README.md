# Meta-analysis-SNP
SNP Meta-Analysis of Egyptian and Sham Cohort Genetic Data

Project Description

Egyptian Cohort

This project implements a complete, transparent, and reproducible meta-analysis and clustering workflow for SNP (Single Nucleotide Polymorphism) genetic association studies in Egyptian populations. Using robust statistics and unsupervised machine learning, it enables identification of the most significant genetic variants, evaluation of result consistency, and discovery of hidden patterns in the data. All outputs are publication-ready and support open science.

Key Steps:

Data import and cleaning

Funnel plot (publication bias)

Heterogeneity analysis (Cochran's Q, I²)

SNP ranking (by OR and p-value)

Clustering (KMeans) to find SNP groups

Meta-regression

Automatic export to Excel for downstream analysis

Significance:

Highlights top genetic risk loci in Egyptian cohorts

Reveals heterogeneity and sample-size effects

Supports research in population genetics and precision medicine

Open, reproducible pipeline for further studies

Sham Data Cohort

The same workflow is applied to the sham_data1.xlsx file, enabling comparative or secondary analysis for another population or cohort ("Sham"). The code is fully reusable—simply update the input file and run.

Key Steps:

Load and clean sham_data1.xlsx

Robust numerical conversion for all relevant columns

All downstream analyses identical to Egyptian cohort pipeline

Significance:

Enables multi-cohort or multi-population analysis with no extra code changes

Directly comparable outputs for all cohorts

