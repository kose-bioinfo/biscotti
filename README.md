# BISCOTTI

Msa2Blosum is a tool to generate custom BLOSUM substitution matrices from large scale (>10,000 sequences) multiple sequence alignments (MSAs). 
It applies Henikoff sequence weighting and pseudocount smoothing to produce robust and unbiased scoring matrices suitable for downstream bioinformatics analyses.

# Features:

Generate custom BLOSUM matrices tailored to your dataset

Henikoff weighting to reduce bias from similar sequences

Pseudocount addition (default α=0.1) to avoid zero probabilities

Easy command-line interface

# Installation

You can install Msa2Blosum by cloning the GitHub repository:

```
git clone https://github.com/skose82/MSA2BLOSUM.git
cd MSA2BLOSUM
```

This package requires the following dependencies:

```
pip install -r requirements.txt
```

# Usage

Input: 

alignment.fasta  

file.clstr (at chosen % identity for your MSA -- e.g. CD-HIT output.) 

Output: 

A Blosum Substitution Matrix in .csv format.

Basic usage example to generate a BLOSUM Matrix from an MSA file (alignment.fasta):

```
python msa2blosum.py -i alignment.fasta -o custom_blosum.mat
```


# Citation

If you use Msa2Blosum in your work, please cite:

```
S.H.Kose. (2025). Msa2Blosum: Custom BLOSUM matrix generator from MSAs [Software]. GitHub repository. https://github.com/skose82/Msa2Blosum
```



