# SCHEMA 2.0

Code accompanying the second data freeze of the **SCHEMA (Schizophrenia Exome Meta-Analysis)** project. This repository contains the full analytical pipeline used to process rare variant association study (RVAS) data, run gene- and variant-set-based association analyses, perform enrichment analyses, generate the manuscript figures, and compare SCHEMA results against common-variant (GWAS) findings.

## Overview

SCHEMA aggregates exome sequencing data across schizophrenia case-control cohorts to identify genes and variant classes enriched for rare, damaging variation in cases. This repository covers five stages of the analysis:

1. **RVAS prep and Analysis** — numbered code for annotation, preparation, and CMH analysis of per-gene RVAS
3. **Enrichment analysis** — enrichment of rare variants by functional annotation and within ASD and DD/ID gene sets
4. **Figures** — code to reproduce manuscript figures
5. **GWAS overlap** — integration and comparison of RVAS signal with common-variant GWAS results 

## Data Availability

Individual-level sequencing data are subject to data access agreements and are not included in this repository. Summary statistics and gene-level results from the published freeze are available at https://schema.broadinstitute.org

## Citation

If you use this code, please cite:

> Sealock et al, 2026
> preprint link: 
