from collections import defaultdict
from Bio import SeqIO
import numpy as np
import pandas as pd
import math
import re
import multiprocessing as mp
from tqdm import tqdm
import time

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
AA_SET = set(AA_LIST)
SCALE = 1
PSEUDOCOUNT = 0.01  # Small pseudocount to avoid zeros in freq calculations


def parse_cd_hit_clstr(file_path):
    clusters = []
    current = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith(">Cluster"):
                if current:
                    clusters.append(current)
                current = []
            else:
                seq_id = line.split('>')[1].split('...')[0]
                current.append(seq_id)
        if current:
            clusters.append(current)
    return clusters


def load_msa(fasta_path):
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}

def henikoff_weights(seqs):
    n = len(seqs)
    L = len(seqs[0])
    weights = np.zeros(n, dtype=float)
    
    arr = np.array([list(seq) for seq in seqs])
    
    for pos in range(L):
        col = arr[:, pos]
        valid_mask = np.isin(col, list(AA_SET))
        valid_aa = col[valid_mask]
        if valid_aa.size == 0:
            continue
        r = len(set(valid_aa))
        
        unique, counts = np.unique(valid_aa, return_counts=True)
        counts_dict = dict(zip(unique, counts))
        
        for i in range(n):
            aa = col[i]
            if aa in counts_dict:
                weights[i] += 1 / (r * counts_dict[aa])
    total_weight = np.sum(weights)
    if total_weight > 0:
        return weights / total_weight
    else:
        return weights

def count_pairs_in_cluster(args):
    cluster, msa_dict = args
    start = time.time()
    
    seqs = [msa_dict[seq_id] for seq_id in cluster if seq_id in msa_dict]
    if len(seqs) < 2 or any(len(seq) != len(seqs[0]) for seq in seqs):
        print(f"Skipped cluster size {len(cluster)} (too small or inconsistent seq lengths)")
        return defaultdict(float)
    
    weights = henikoff_weights(seqs)
    L = len(seqs[0])
    
    arr = np.array([list(seq) for seq in seqs])
    n = len(seqs)
    
    pair_counts = defaultdict(float)
    
    for pos in range(L):
        col = arr[:, pos]
        valid_idx = np.where(np.isin(col, list(AA_SET)))[0]
        if len(valid_idx) < 2:
            continue
        
        valid_aa = col[valid_idx]
        valid_weights = weights[valid_idx]
        
        aa_weight_sum = defaultdict(float)
        for aa, w in zip(valid_aa, valid_weights):
            aa_weight_sum[aa] += w
        
        aa_items = sorted(aa_weight_sum.items())
        
        for i, (aa1, w1) in enumerate(aa_items):
            pair_counts[(aa1, aa1)] += w1 * w1
            for j in range(i + 1, len(aa_items)):
                aa2, w2 = aa_items[j]
                pair_counts[(aa1, aa2)] += 2 * w1 * w2
    
    duration = time.time() - start
    print(f"Processed cluster size {len(cluster)} in {duration:.2f} seconds")
    return pair_counts

def aggregate_pair_counts(pair_counts_list):
    total = defaultdict(float)
    for pc in pair_counts_list:
        for k, v in pc.items():
            total[k] += v
    return total

def compute_log_odds_matrix(pair_counts):
    total_pairs = sum(pair_counts.values())
    if total_pairs == 0:
        raise ValueError("No pairs counted. Check your input data.")

    aa_counts = defaultdict(float)
    for (a1, a2), count in pair_counts.items():
        if a1 == a2:
            aa_counts[a1] += 2 * count
        else:
            aa_counts[a1] += count
            aa_counts[a2] += count
    total_aa = sum(aa_counts.values())

    for aa in AA_LIST:
        aa_counts[aa] += PSEUDOCOUNT
    total_aa += PSEUDOCOUNT * len(AA_LIST)

    aa_freq = {aa: aa_counts[aa] / total_aa for aa in AA_LIST}

    n_pairs_possible = (len(AA_LIST) * (len(AA_LIST) + 1)) // 2
    adjusted_total_pairs = total_pairs + PSEUDOCOUNT * n_pairs_possible

    matrix = pd.DataFrame(index=AA_LIST, columns=AA_LIST, dtype=int)

    for a1 in AA_LIST:
        for a2 in AA_LIST:
            pair = tuple(sorted((a1, a2)))
            obs_count = pair_counts.get(pair, 0) + PSEUDOCOUNT
            obs_freq = obs_count / adjusted_total_pairs

            if a1 == a2:
                exp_freq = aa_freq[a1] ** 2
            else:
                exp_freq = 2 * aa_freq[a1] * aa_freq[a2]

            if exp_freq > 0 and obs_freq > 0:
                score = math.log2(obs_freq / exp_freq)
            else:
                score = -10.0

            lod = int(round(SCALE * score))
            matrix.at[a1, a2] = matrix.at[a2, a1] = lod

    return matrix

def build_blosum_parallel(msa_dict, clusters, nproc=8):
    pair_counts_list = []
    with mp.Pool(nproc) as pool:
        for result in tqdm(pool.imap_unordered(count_pairs_in_cluster, [(cluster, msa_dict) for cluster in clusters], chunksize=1), total=len(clusters), desc="Clusters processed"):
            pair_counts_list.append(result)
    pair_counts = aggregate_pair_counts(pair_counts_list)
    return compute_log_odds_matrix(pair_counts)

# --- Run ---
if __name__ == "__main__":
    msa_file = "RSV_A.fasta"
    clstr_file = "clustered_RSVA90.clean.clstr"
    output_file = "RSVA_F_blosum_90__test___v5_scale0_psuedo_0.01_Clean2.csv"
    nproc = 6  # Adjust for your CPU

    print("Loading MSA...")
    msa = load_msa(msa_file)
    print(f"Loaded {len(msa)} sequences")

    print("Loading clusters...")
    clusters = parse_cd_hit_clstr(clstr_file)
    print(f"Loaded {len(clusters)} clusters")
    #print("Cluster sizes:", [len(c) for c in clusters])

    print("Computing matrix (parallel with progress)...")
    start = time.time()
    matrix = build_blosum_parallel(msa, clusters, nproc=nproc)
    print(f"Done in {time.time() - start:.2f} seconds")

    matrix.to_csv(output_file)
    print(f"Matrix saved to {output_file}")

