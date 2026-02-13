# VariantAnalysis

This repository contains the computational workflows, custom scripts, and analysis code used in the manuscript:

**Integrating Long-Read Structural Variant Analysis with single-nucleus RNA-seq to Elucidate Gene Expression Effects in Disease**

## Overview
This project integrates high-fidelity (HiFi) long-read whole-genome sequencing (WGS) with single-nucleus RNA-sequencing (snRNA-seq) to identify Structural Variants (SVs) that act as expression quantitative trait loci (eQTLs) or drive allele-specific expression (ASE) in Parkinson’s Disease (PD).

The code in this repository covers:
1.  **Individual-specific SV detection:** Ensemble calling using multiple methods.
1.  **Cohort-level SV Processing:** Merging callsets, filtering, and quality control.
2.  **Genotyping:** Re-assign confident genotypes.
3.  **Functional Analysis:** *cis*-eQTL mapping and ASE analysis.
