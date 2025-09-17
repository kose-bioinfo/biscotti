<img src="https://github.com/user-attachments/assets/3f0993dc-54ec-4301-bf27-fa961159cf13" alt="BISCOTTI logo" width="100" height="100" align="left" />

# BISCOTTI

_Crunching alignments into your perfect BLOSUM_

<p align="left">
  <!-- Release version badge -->
  <img src="https://img.shields.io/badge/release-v1.4.0-blue?style=flat-square" alt="Release v1.4.0">
  <!-- Development language badge -->
  <img src="https://img.shields.io/badge/language-python-yellow?style=flat-square" alt="Language: Python">
  <!-- MIT License badge -->
  <img src="https://img.shields.io/badge/license-MIT-green?style=flat-square" alt="License: MIT">
  <!-- Code style: black badge -->
  <img src="https://img.shields.io/badge/code%20style-black-black?style=flat-square" alt="Code style: black">
</p>

BISCOTTI or BLOSUM Individualised Substitution Calculator Optimised for Tailored Translation Indices is a Python tool for building custom substitution matrices from multiple sequence alignments (MSAs).

...
# Features

Loads large MSAs (10k+ sequences).

Uses CD-HIT clustering to reduce redundancy.

Applies Henikoff sequence weighting within clusters.

Aggregates weighted amino acid pair counts.

Adds pseudocounts to avoid zero probabilities.

Computes log-odds substitution scores (BLOSUM style).

Outputs a symmetric 20×20 substitution matrix in CSV.

# Usage

QC Filter and align your sequences (e.g with mafft) 

Cluster your FASTA sequences with CD-HIT (e.g., 90% identity for highly conserved proteins):

```
cd-hit -i RSV_A.fasta -o clustered_RSVA90.fasta -c 0.90 -n 5 -d 0
```

Run the script:

```
python3 biscotti.py
```

Output:

RSVA_F_blosum_90.csv — custom 20×20 substitution matrix.

# Algorithm/rationale

```
1. Henikoff weighting

For each sequence i, weight is computed across alignment positions:

w_i = Σ_p ( 1 / ( r_p * n_ap ) )

r_p = number of unique residues at position p

n_ap = count of amino acid a at position p

Normalise weights:

w_i = w_i / Σ_i w_i

2. Weighted pair counting (per cluster)

For each cluster:

Iterate through alignment positions.

Collect weighted amino acid counts.

Count all unique residue pairs (aa1, aa2):

C(aa1, aa2) += 2 * w1 * w2    if aa1 != aa2
C(aa, aa)   += w^2            for aa == aa

3: Combine clusters
C_total(aa1, aa2) = Σ_clusters C_cluster(aa1, aa2)

4: Background frequencies

From symmetric pair counts:

A[a] = 2*C[a,a] + Σ_{b != a} C[a,b]

Normalise:

P_obs(a) = A[a] / Σ_z A[z]

5: Observed vs expected pair probabilities

Observed pair probabilities:

P_obs(a,b) = C_total(a,b) / Σ_{x<=y} C_total(x,y)

Expected probabilities under independence:

P_exp(a,b) = P_obs(a)^2                if a == b
P_exp(a,b) = 2 * P_obs(a) * P_obs(b)   if a != b

Step 6: Log-odds substitution scores
S(a,b) = log2( P_obs(a,b) / P_exp(a,b) )
LOD(a,b) = round( SCALE * S(a,b) )

Scores are symmetric (S(a,b) = S(b,a)).
```
Pseudocounts (α = 0.1) are added to avoid zero probabilities.

Final result is a 20×20 CSV matrix (rows and columns in amino acid order).

# References

[Paper] 

Henikoff, S., & Henikoff, J. G. (1992). Amino acid substitution matrices from protein blocks. PNAS, 89(22), 10915–10919. PMID: 1438297

Li, W., & Godzik, A. (2006). Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics, 22(13), 1658–1659.
