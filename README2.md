BISCOTTI— 

BISCOTTI is a Python tool for building custom substitution matrices from multiple sequence alignments (MSAs).
It is designed for highly conserved viral proteins such as RSV F, but can be applied to any protein alignment.

Features

Loads large MSAs (10k+ sequences).

Uses CD-HIT clustering to reduce redundancy.

Applies Henikoff sequence weighting within clusters.

Aggregates weighted amino acid pair counts.

Adds pseudocounts to avoid zero probabilities.

Computes log-odds substitution scores (BLOSUM style).

Outputs a symmetric 20×20 substitution matrix in CSV.

Usage

Cluster your FASTA sequences with CD-HIT (e.g., 90% identity):

cd-hit -i RSV_A.fasta -o clustered_RSVA90.fasta -c 0.90 -n 5 -d 0


Run the script:

python build_matrix.py


Output:

RSVA_F_blosum_90.csv — custom 20×20 substitution matrix.

Algorithm
Step 1: Henikoff weighting

For each sequence i, weight is computed across alignment positions:

w_i = Σ_p ( 1 / ( r_p * n_ap ) )


r_p = number of unique residues at position p

n_ap = count of amino acid a at position p

Normalize weights:

w_i = w_i / Σ_i w_i

Step 2: Weighted pair counting (per cluster)

For each cluster:

Iterate through alignment positions.

Collect weighted amino acid counts.

Count all unique residue pairs (aa1, aa2):

C(aa1, aa2) += 2 * w1 * w2    if aa1 != aa2
C(aa, aa)   += w^2            for aa == aa

Step 3: Combine clusters
C_total(aa1, aa2) = Σ_clusters C_cluster(aa1, aa2)

Step 4: Background frequencies

From symmetric pair counts:

A[a] = 2*C[a,a] + Σ_{b != a} C[a,b]


Normalize:

P_obs(a) = A[a] / Σ_z A[z]

Step 5: Observed vs expected pair probabilities

Observed pair probabilities:

P_obs(a,b) = C_total(a,b) / Σ_{x<=y} C_total(x,y)


Expected probabilities under independence:

P_exp(a,b) = P_obs(a)^2                if a == b
P_exp(a,b) = 2 * P_obs(a) * P_obs(b)   if a != b

Step 6: Log-odds substitution scores
S(a,b) = log2( P_obs(a,b) / P_exp(a,b) )
LOD(a,b) = round( SCALE * S(a,b) )


Scores are symmetric (S(a,b) = S(b,a)).

Pseudocounts (α = 0.1) are added to avoid zero probabilities.

Final result is a 20×20 CSV matrix (rows and columns in amino acid order).


📚 References

Henikoff, S., & Henikoff, J. G. (1992). Amino acid substitution matrices from protein blocks. PNAS, 89(22), 10915–10919. PMID: 1438297

Li, W., & Godzik, A. (2006). Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics, 22(13), 1658–1659.
