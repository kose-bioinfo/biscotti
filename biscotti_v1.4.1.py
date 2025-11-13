# BISCOTTI: Crunching Alignments into Perfect BLOSUM Matrices
# Copyright (C) 2025 S.H.Kose
# Licensed under the GNU General Public License v3.0
# See LICENSE file for details

#!/usr/bin/env python3
from collections import defaultdict
from Bio import SeqIO
import numpy as np
import pandas as pd
import math
import multiprocessing as mp
from tqdm import tqdm
import time
import argparse


# Config

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
AA_SET = set(AA_LIST)
SCALE = 1
PSEUDOCOUNT = 0.1             # pseudocount to avoid zeros
ADD_MSA_SINGLETONS = True    # True => add singleton clusters for base IDs not in .clstr
NPROC_DEFAULT = 4


# Parsing & loading


def parse_cd_hit_clstr_first_id(file_path: str):
  
    clusters, cur = [], []
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">Cluster"):
                if cur:
                    seen = set()
                    cur = [x for x in cur if not (x in seen or seen.add(x))]
                    clusters.append(cur)
                cur = []
                continue
            if ">" in line:
                sid = line.split(">", 1)[1].split()[0].rstrip(",")
                if sid.endswith("..."):
                    sid = sid[:-3]
                cur.append(sid)
    if cur:
        seen = set()
        cur = [x for x in cur if not (x in seen or seen.add(x))]
        clusters.append(cur)
    return clusters


def load_msa_instances(fasta_path: str):
    seq_by_key = {}
    id_to_keys = defaultdict(list)
    seen = defaultdict(int)

    for rec in SeqIO.parse(fasta_path, "fasta"):
        base = rec.id
        if base in id_to_keys:
            seen[base] += 1
            ukey = f"{base}_{seen[base]}"
        else:
            ukey = base
        seq_by_key[ukey] = str(rec.seq)
        id_to_keys[base].append(ukey)

    return seq_by_key, id_to_keys


def expand_clusters_with_instances(clusters_base_ids, id_to_keys, msa_keys_set=None):
    
    expanded = []
    for cl in clusters_base_ids:
        members = []
        for base in cl:
            for ukey in id_to_keys.get(base, []):
                if msa_keys_set is None or ukey in msa_keys_set:
                    members.append(ukey)
        seen = set()
        members = [x for x in members if not (x in seen or seen.add(x))]
        if members:
            expanded.append(members)
    return expanded


def add_missing_as_singletons(expanded_clusters, seq_by_key, id_to_keys, clusters_base_ids):
 
    clustered_bases = set(b for cl in clusters_base_ids for b in cl)
    all_bases = set(id_to_keys.keys())
    missing_bases = sorted(all_bases - clustered_bases)
    if not missing_bases:
        return expanded_clusters, 0

    for base in missing_bases:
        for ukey in id_to_keys[base]:
            if ukey in seq_by_key:
                expanded_clusters.append([ukey])
    return expanded_clusters, sum(len(id_to_keys[b]) for b in missing_bases)


# Henikoff weights & counts


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
                weights[i] += 1.0 / (r * counts_dict[aa])

    tot = np.sum(weights)
    return weights / tot if tot > 0 else weights


def count_pairs_in_cluster(args):
    cluster_keys, seq_by_key = args
    start = time.time()
    seqs = [seq_by_key[k] for k in cluster_keys if k in seq_by_key]
    if len(seqs) < 2:
        return defaultdict(float)
    L = len(seqs[0])
    if any(len(s) != L for s in seqs):
        return defaultdict(float)

    weights = henikoff_weights(seqs)
    arr = np.array([list(seq) for seq in seqs])
    pair_counts = defaultdict(float)

    for pos in range(L):
        col = arr[:, pos]
        valid_idx = np.where(np.isin(col, list(AA_SET)))[0]
        if len(valid_idx) < 2:
            continue
        valid_aa = col[valid_idx]
        valid_w = weights[valid_idx]

        aa_weight_sum = defaultdict(float)
        for aa, w in zip(valid_aa, valid_w):
            aa_weight_sum[aa] += w

        aa_items = sorted(aa_weight_sum.items())
        for i, (aa1, w1) in enumerate(aa_items):
            pair_counts[(aa1, aa1)] += w1 * w1
            for j in range(i + 1, len(aa_items)):
                aa2, w2 = aa_items[j]
                pair_counts[(aa1, aa2)] += 2.0 * w1 * w2

    print(f"Processed cluster size {len(cluster_keys)} in {time.time() - start:.2f} seconds")
    return pair_counts


def aggregate_pair_counts(pair_counts_list):
    total = defaultdict(float)
    for pc in pair_counts_list:
        for k, v in pc.items():
            total[k] += v
    return total


# Log-odds matrix


def compute_log_odds_matrix(pair_counts):
    total_pairs = sum(pair_counts.values())
    if total_pairs == 0:
        raise ValueError("No pairs counted. Check input data.")

    aa_counts = defaultdict(float)
    for (a1, a2), count in pair_counts.items():
        if a1 == a2:
            aa_counts[a1] += 2.0 * count
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
            obs_count = pair_counts.get(pair, 0.0) + PSEUDOCOUNT
            obs_freq = obs_count / adjusted_total_pairs
            if a1 == a2:
                exp_freq = aa_freq[a1] ** 2
            else:
                exp_freq = 2.0 * aa_freq[a1] * aa_freq[a2]
            score = math.log2(obs_freq / exp_freq) if (exp_freq > 0 and obs_freq > 0) else -10.0
            lod = int(round(SCALE * score))
            matrix.at[a1, a2] = matrix.at[a2, a1] = lod
    return matrix


# Parallel build

def build_blosum_parallel(seq_by_key, clusters_keys, nproc=NPROC_DEFAULT):
    pair_counts_list = []
    with mp.Pool(nproc) as pool:
        iterable = [(cluster, seq_by_key) for cluster in clusters_keys]
        for result in tqdm(pool.imap_unordered(count_pairs_in_cluster, iterable, chunksize=1),
                           total=len(clusters_keys), desc="Clusters processed"):
            pair_counts_list.append(result)
    pair_counts = aggregate_pair_counts(pair_counts_list)
    return compute_log_odds_matrix(pair_counts)


# Main

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build a BLOSUM matrix from MSA and clusters")

    # All arguments are required
    parser.add_argument("--msa_file", required=True, help="Input MSA file")
    parser.add_argument("--clstr_file", required=True, help="Cluster file")
    parser.add_argument("--output_file", required=True, help="Output CSV file")
    parser.add_argument("--nproc", type=int, required=False, help="Number of processes for parallel computation, default = 4")
    parser.add_argument("--pseudocount", type=float, required=True, help="Input pseudocount to account for zero values")
    args = parser.parse_args()

    msa_file = args.msa_file
    clstr_file = args.clstr_file
    output_file = args.output_file
    nproc = args.nproc
    pseudocount = args.pseudocount

    print("Loading MSA (ALL instances)...")
    seq_by_key, id_to_keys = load_msa_instances(msa_file)
    print(f"Loaded {len(seq_by_key)} records. Unique base IDs: {len(id_to_keys)}")

    print("Loading clusters...")
    clusters_base = parse_cd_hit_clstr_first_id(clstr_file)
    print(f"Loaded {len(clusters_base)} clusters (base IDs)")

    clusters_keys = expand_clusters_with_instances(clusters_base, id_to_keys, msa_keys_set=set(seq_by_key.keys()))
    if ADD_MSA_SINGLETONS:
        clusters_keys, added = add_missing_as_singletons(clusters_keys, seq_by_key, id_to_keys, clusters_base)
        if added:
            print(f"Added {added} singleton records for base IDs missing from clusters.")

    covered_keys = len({k for cl in clusters_keys for k in cl})
    print(f"Clusters ready: {len(clusters_keys)} | Coverage: {covered_keys} / {len(seq_by_key)} records")
    print("Cluster sizes (first 10):", [len(c) for c in clusters_keys[:10]])

    print("Computing matrix...")
    start = time.time()
    matrix = build_blosum_parallel(seq_by_key, clusters_keys, nproc=nproc)
    print(f"Done in {time.time() - start:.2f} seconds")

    matrix.to_csv(output_file)
    print(f"Matrix saved to {output_file}")
