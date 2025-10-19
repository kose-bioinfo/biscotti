<img src="https://github.com/user-attachments/assets/b33c76f0-0ae7-41b0-aaf4-2ffb58b28be1" alt="italian-cantuccini-biscotti-stacked" width="120">

# BISCOTTI — Crunching Alignments into Perfect BLOSUM Matrices

<p align="left">
  <img src="https://img.shields.io/badge/release-v1.4.0-blue?style=flat-square" alt="Release v1.4.0">
  <img src="https://img.shields.io/badge/language-Python-yellow?style=flat-square" alt="Language: Python">
  <img src="https://img.shields.io/badge/license-MIT-green?style=flat-square" alt="License: MIT">
  <img src="https://img.shields.io/badge/code%20style-black-black?style=flat-square" alt="Code style: black">
  <a href="https://kose-bioinfo.github.io/biscotti/">
    <img src="https://img.shields.io/badge/View%20MathJax%20Docs-Click%20Here-blue?style=flat-square" alt="MathJax Documentation">
  </a>
</p>

BISCOTTI is a high-performance Python tool for generating custom substitution matrices from large-scale multiple sequence alignments (>30k sequences). This tool was developed by S.H.Kose in collaboration with the Tregoning Lab @ Imperial College London. It features redundancy reduction, Henikoff weighting, and BLOSUM-style log-odds scoring.  

## Overview

BISCOTTI implements:

- Redundancy reduction via CD-HIT clustering  
- Henikoff sequence weighting within clusters  
- Weighted amino acid pair counting  
- Computation of log-odds substitution scores (BLOSUM-style)  
- Output of symmetric 20×20 substitution matrices in CSV format  
---

## Documentation


BISCOTTI adapts the widely used **[Henikoff method](https://pmc.ncbi.nlm.nih.gov/articles/PMC50453/)**  for constructing BLOSUM-like matrices by incorporating weighted sequence contributions and custom scoring options. While Henikoff's original approach laid the foundation for position-specific residue weighting, this implementation extends it to handle large datasets efficiently and allows users to customise the scoring scheme, making it more flexible for diverse protein families.

The full mathematical formulation (Henikoff weighting, pair counting, probability calculations, and log-odds scoring) is available on the MathJax-rendered site:

**[View the Biscotti Substitution Matrix MathJax Page](https://kose-bioinfo.github.io/biscotti/)** 

This page provides a complete explanation of the equations and their derivations with properly formatted LaTeX.

---


## Usage

### QC and Align Sequences
Filter and align your sequences using a tool like MAFFT or CLUSTAL W.

### Cluster Sequences
Cluster your FASTA sequences using CD-HIT. For highly conserved proteins, use ~90% identity cutoff:


`cd-hit -i RSV_A.fasta -o clustered_RSVA90.fasta -c 0.90 -n 5 -d 0`


Run BISCOTTI

`python3 biscotti.py --msa_file MSA_FILE --clstr_file CLSTR_FILE --output_file OUTPUT_FILE --nproc 1-10`

**NOTE*: Although BISCOTTI is optimised for efficiency, its performance scales with the number of available CPU cores, with more cores yielding faster computation.


### Output:

`foo.csv `


## Rationale

Biscotti takes aligned sequences, applies **Henikoff weighting** to correct for overrepresentation, counts weighted amino acid pairs, and computes **log-odds substitution scores** with pseudocounts. The output is a robust substitution matrix suitable for downstream bioinformatics analyses.
please see docs above for details. 


## Associated Publication

This repository contains the BLOSUM matrix software developed for the following publication:

Mosscrop, L. G.†, Gerardi, V.†, Kose, S. H., Talts, T., Thomas, C., Paschos, K., Brown, J., Williams, T. C., Skinner, M., Zambon, M., Bravi, B., & Tregoning, J. S.
An integrated in silico and in vitro genotype-to-phenotype pipeline to predict and characterise RSV F site zero escape mutants. (2025, under review)

## Citation 

If you use this code or adapt components in your research, please cite:

**Kose, S. H., Gerardi, V.,Bravi, B.,Tregoning, J.S (2025).** *BISCOTTI: A Scalable Tool for Custom BLOSUM Matrix Construction.* GitHub repository: [https://github.com/kose-bioinfo/biscotti](https://github.com/kose-bioinfo/biscotti)

### Acknowledgements

Special thanks to:  

- [Barbara Bravi](https://github.com/bbravi) – for conceptual guidance, and suggesting the log-odds scoring approach for BLOSUM matrices.  
- Valeria Gerardi – for data preparation and analysis.


## Other References

Henikoff, S., & Henikoff, J. G. (1992). Amino acid substitution matrices from protein blocks. PNAS, 89(22), 10915–10919. PMID: 1438297

Li, W., & Godzik, A. (2006). Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics, 22(13), 1658–1659
