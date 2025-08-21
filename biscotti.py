from collections import defaultdict
from Bio import SeqIO
import numpy as np
import pandas as pd
import math
import multiprocessing as mp
from tqdm import tqdm
import time

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
AA_SET = set(AA_LIST)
SCALE = 1
PSEUDOCOUNT = 0.1          # small pseudocount to avoid zeros
ADD_SINGLETONS = True     # set True to add MSA-only sequences as singleton clusters

# ---------------------------
# .clstr parsing & sanitizing
# ---------------------------

def parse_cd_hit_clstr_first_id(file_path):
    """
    Parse CD-HIT .clstr:
      - start new cluster at lines beginning with '>Cluster'
      - on member lines, extract ONLY the first ID after the first '>'
      - strip any trailing '...' (defensive)
      - de-duplicate IDs within a cluster (order-preserving)
    Returns: list[list[str]]
    """
    clusters, cur = [], []
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">Cluster"):
                if cur:
                    # dedupe while preserving order
                    seen = set()
                    cur = [x for x in cur if not (x in seen or seen.add(x))]
                    clusters.append(cur)
                cur = []
                continue
            if ">" in line:
                # take only the first ID token after first '>'
                sid = line.split(">", 1)[1].split()[0].rstrip(",")
                if sid.endswith("..."):  # defensive
                    sid = sid[:-3]
                cur.append(sid)
    if cur:
        seen = set()
        cur = [x for x in cur if not (x in seen or seen.add(x))]
        clusters.append(cur)
    return clusters

def load_msa(fasta_path):
    # Biopython record.id is the first whitespace-delimited token (good for CD-HIT matching)
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}

def sanitize_clusters(clusters_raw, msa_keys):
    """
    Intersect cluster members with the MSA keys and drop empty clusters.
    Ensures cluster sizes never exceed the MSA size and removes stray IDs.
    """
    clean = []
    for cl in clusters_raw:
        members = [sid for sid in cl if sid in msa_keys]
        if members:
            clean.append(members)
    return clean

# ---------------------------
# Henikoff weighting & counts
# ---------------------------

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
                # Henikoff position weight contribution
                weights[i] += 1.0 / (r * counts_dict[aa])

    total_weight = np.sum(weights)
    return weights / total_weight if total_weight > 0 else weights

def count_pairs_in_cluster(args):
    cluster, msa_dict = args
    start = time.time()

    # Pull sequences for members present in MSA (cluster already sanitized, but keep defensive)
    seqs = [msa_dict[sid] for sid in cluster if sid in msa_dict]
    if len(seqs) < 2:
        # too small to contribute pairs
        return defaultdict(float)
    # require consistent lengths
    L = len(seqs[0])
    if any(len(s) != L for s in seqs):
        # skip inconsistent clusters
        return defaultdict(float)

    weights = henikoff_weights(seqs)
    arr = np.array([list(seq) for seq in seqs])

    pair_counts = defaultdict(float)

    # per-column weighted AA co-occurrence (ignore gaps / non-AA)
    for pos in range(L):
        col = arr[:, pos]
        valid_idx = np.where(np.isin(col, list(AA_SET)))[0]
        if len(valid_idx) < 2:
            continue

        valid_aa = col[valid_idx]
        valid_w  = weights[valid_idx]

        # aggregate weights by residue at this column
        aa_weight_sum = defaultdict(float)
        for aa, w in zip(valid_aa, valid_w):
            aa_weight_sum[aa] += w

        aa_items = sorted(aa_weight_sum.items())  # sort for deterministic ordering

        # add diagonal and off-diagonal contributions; keep pair keys sorted
        for i, (aa1, w1) in enumerate(aa_items):
            pair_counts[(aa1, aa1)] += w1 * w1
            for j in range(i + 1, len(aa_items)):
                aa2, w2 = aa_items[j]
                pair_counts[(aa1, aa2)] += 2.0 * w1 * w2

    duration = time.time() - start
    print(f"Processed cluster size {len(cluster)} in {duration:.2f} seconds")
    return pair_counts

def aggregate_pair_counts(pair_counts_list):
    total = defaultdict(float)
    for pc in pair_counts_list:
        for k, v in pc.items():
            total[k] += v
    return total

# ---------------------------
# Log-odds matrix
# ---------------------------

def compute_log_odds_matrix(pair_counts):
    total_pairs = sum(pair_counts.values())
    if total_pairs == 0:
        raise ValueError("No pairs counted. Check your input data (ID matching, cluster parsing).")

    # background AA counts from pair counts (symmetric)
    aa_counts = defaultdict(float)
    for (a1, a2), count in pair_counts.items():
        if a1 == a2:
            aa_counts[a1] += 2.0 * count
        else:
            aa_counts[a1] += count
            aa_counts[a2] += count
    total_aa = sum(aa_counts.values())

    # pseudocounts on single-AA background
    for aa in AA_LIST:
        aa_counts[aa] += PSEUDOCOUNT
    total_aa += PSEUDOCOUNT * len(AA_LIST)

    aa_freq = {aa: aa_counts[aa] / total_aa for aa in AA_LIST}

    # pseudocounts on pair space
    n_pairs_possible = (len(AA_LIST) * (len(AA_LIST) + 1)) // 2
    adjusted_total_pairs = total_pairs + PSEUDOCOUNT * n_pairs_possible

    matrix = pd.DataFrame(index=AA_LIST, columns=AA_LIST, dtype=int)

    for a1 in AA_LIST:
        for a2 in AA_LIST:
            pair = tuple(sorted((a1, a2)))
            obs_count = pair_counts.get(pair, 0.0) + PSEUDOCOUNT
            obs_freq = obs_count / adjusted_total_pairs

            # expected under independence; off-diagonal has factor 2
            if a1 == a2:
                exp_freq = aa_freq[a1] ** 2
            else:
                exp_freq = 2.0 * aa_freq[a1] * aa_freq[a2]

            score = math.log2(obs_freq / exp_freq) if (exp_freq > 0 and obs_freq > 0) else -10.0
            lod = int(round(SCALE * score))
            matrix.at[a1, a2] = matrix.at[a2, a1] = lod

    return matrix

# ---------------------------
# Parallel build
# ---------------------------

def build_blosum_parallel(msa_dict, clusters, nproc=8):
    # intersect clusters with MSA keys (safety)
    msa_keys = set(msa_dict.keys())
    clusters = sanitize_clusters(clusters, msa_keys)

    # optional: add singleton clusters for any MSA sequences not covered
    if ADD_SINGLETONS:
        covered = {sid for cl in clusters for sid in cl}
        missing = sorted(msa_keys - covered)
        clusters += [[sid] for sid in missing]

    pair_counts_list = []
    with mp.Pool(nproc) as pool:
        iterable = [(cluster, msa_dict) for cluster in clusters]
        for result in tqdm(pool.imap_unordered(count_pairs_in_cluster, iterable, chunksize=1),
                           total=len(clusters), desc="Clusters processed"):
            pair_counts_list.append(result)
    pair_counts = aggregate_pair_counts(pair_counts_list)
    return compute_log_odds_matrix(pair_counts)

# ---------------------------
# Run
# ---------------------------

if __name__ == "__main__":
    msa_file = "RSV_B.fasta"
    clstr_file = "clustered_RSVB90_remade.fasta.clstr"
    output_file = "RSVB_F_blosum_90_Sc1_pseudo0.1_singltons_v7.csv"
    nproc = 6  # Adjust for your CPU

    print("Loading MSA...")
    msa = load_msa(msa_file)
    print(f"Loaded {len(msa)} sequences (unique first-token IDs)")

    print("Loading clusters...")
    clusters_raw = parse_cd_hit_clstr_first_id(clstr_file)
    msa_keys = set(msa.keys())

    # sanitize + report coverage
    clusters = sanitize_clusters(clusters_raw, msa_keys)
    covered = len({sid for cl in clusters for sid in cl})
    missing = len(msa_keys - {sid for cl in clusters for sid in cl})
    print(f"Loaded {len(clusters)} clusters covering {covered} of {len(msa)} MSA sequences "
          f"({missing} not in clusters)")
    print("Cluster sizes (sanitized):", [len(c) for c in clusters])

    if ADD_SINGLETONS and missing:
        print(f"Added {missing} singleton clusters (full coverage).")

    print("Computing matrix (parallel with progress)...")
    start = time.time()
    matrix = build_blosum_parallel(msa, clusters, nproc=nproc)
    print(f"Done in {time.time() - start:.2f} seconds")

    matrix.to_csv(output_file)
    print(f"Matrix saved to {output_file}")
