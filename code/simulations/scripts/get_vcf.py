#!/usr/bin/env python
import random
import numpy as np
import pandas as pd
import tskit
import logging
import argparse
import sys
import atexit
import warnings
import msprime
import time
import pyslim
from multiprocessing import Process, Pipe

def equal_stratified_host_sample(df, n, rng, group_col="subpop", host_col="host", fitness_col="fitness"):
    groups = sorted(df[group_col].unique())
    if len(groups) == 1:
        valid = df[df[fitness_col] > 0]
        if len(valid) < n:
            raise ValueError(f"Only {len(valid)} eligible samples found in sole group '{groups[0]}', need n={n}.")
        w = valid[fitness_col].to_numpy(dtype=float)
        p = w / w.sum()
        idx = rng.choice(valid.index.values, size=n, replace=False, p=p)
        sampled_df = valid.loc[idx]
    elif len(groups) == 2:
        k = n // 2
        if n % 2 != 0:
            raise ValueError(f"Sample size n={n} not divisible by 2 for stratified sampling.")
        groupA, groupB = groups
        subA_valid = df[(df[group_col] == groupA) & (df[fitness_col] > 0)]
        subB_valid = df[(df[group_col] == groupB) & (df[fitness_col] > 0)]
        nA = len(subA_valid)
        nB = len(subB_valid)
        # Both have at least k
        if nA >= k and nB >= k:
            wA = subA_valid[fitness_col].to_numpy(dtype=float)
            pA = wA / wA.sum()
            idxA = rng.choice(subA_valid.index.values, size=k, replace=False, p=pA)
            sampledA = subA_valid.loc[idxA]
            wB = subB_valid[fitness_col].to_numpy(dtype=float)
            pB = wB / wB.sum()
            idxB = rng.choice(subB_valid.index.values, size=k, replace=False, p=pB)
            sampledB = subB_valid.loc[idxB]
            sampled_df = pd.concat([sampledA, sampledB])
        # A is short, B fills in
        elif nA < k and nB >= n - nA and nA > 0:
            sampledA = subA_valid
            wB = subB_valid[fitness_col].to_numpy(dtype=float)
            pB = wB / wB.sum()
            idxB = rng.choice(subB_valid.index.values, size=n - nA, replace=False, p=pB)
            sampledB = subB_valid.loc[idxB]
            sampled_df = pd.concat([sampledA, sampledB])
        # B is short, A fills in
        elif nB < k and nA >= n - nB and nB > 0:
            sampledB = subB_valid
            wA = subA_valid[fitness_col].to_numpy(dtype=float)
            pA = wA / wA.sum()
            idxA = rng.choice(subA_valid.index.values, size=n - nB, replace=False, p=pA)
            sampledA = subA_valid.loc[idxA]
            sampled_df = pd.concat([sampledA, sampledB])
        else:
            raise ValueError(f"Not enough eligible hosts after fitness filter: {nA} + {nB} < n={n}")
    else:
        raise ValueError(f"Function only supports 1 or 2 groups, found {len(groups)}.")
    return sampled_df.sample(frac=1, random_state=rng).reset_index(drop=True)[host_col].tolist()

def equal_stratified_host_sample_unique(df, n, rng, group_col, host_col, filter, max_tries=10000):
    num_groups = df[group_col].nunique()
    if num_groups != 1:
        raise ValueError(f"Cannot filter for multiple groups")
    k = n // num_groups
    if filter == "parent":
        cols = ["parent1", "parent2"]
    elif filter == "grandparent":
        cols = ["grandparent1", "grandparent2", "grandparent3", "grandparent4"]
    else:
        raise ValueError(f"Unknown filter: {filter}")

    sets = {idx: set(row) for idx, *row in df[cols].itertuples(index=True, name=None)}
    chosen_rows = []
    seen = set()
    idxs_all = df.index.tolist()
    success = False
    for _ in range(max_tries):
        rng.shuffle(idxs_all)
        picks = []
        local_seen = set()
        for idx in idxs_all:
            s = sets[idx]
            if s.isdisjoint(local_seen):
                picks.append(idx)
                local_seen.update(s)
                if len(picks) == n:
                    break
        if len(picks) == n:
            chosen_rows = picks
            success = True
            break
    if not success:
        raise ValueError("Could not select enough hosts globally without shared %s." % filter)
    chosen_rows = rng.permutation(chosen_rows)
    return df.loc[chosen_rows, host_col].tolist()

def process_subsample_iteration(
    seed, host_tree_path, matches_path,
    n, rng, stratify_col, filter, recapitate, ancestral_ne
):
    if recapitate:
        uncapped_host_ts = pyslim.update(tskit.load(host_tree_path))
        host_ts = pyslim.recapitate(uncapped_host_ts, ancestral_Ne=ancestral_ne)
    else:
        host_ts = tskit.load(host_tree_path)
    matches_df = pd.read_csv(matches_path, sep="\t")

    # Sampling logic
    if filter == "none":
        selected_hosts = equal_stratified_host_sample(
            matches_df, n, rng=rng,
            group_col=stratify_col, host_col="host", fitness_col="fitness")
    else:
        selected_hosts = equal_stratified_host_sample_unique(
            matches_df, n, rng=rng, group_col=stratify_col, host_col="host", filter=filter)

    individual_node_sets = [h_ind.nodes for h_ind in host_ts.individuals() if f"h{h_ind.metadata['pedigree_id']}" in selected_hosts]
    all_host_nodes = np.concatenate(individual_node_sets)
    filtered_host_ts = host_ts.simplify(samples=all_host_nodes)
    return filtered_host_ts

def main():
    parser = argparse.ArgumentParser(description="Unified tree permutation test pipeline.")
    parser.add_argument("--output_file", type=str, required=True)
    parser.add_argument("--host_tree", type=str, required=True)
    parser.add_argument("--pairs", type=str, required=True)
    parser.add_argument("--sample_size", type=int, default=200)
    parser.add_argument("--stratify_col", type=str, choices=["subpop", "allele"])
    parser.add_argument("--filter", type=str, choices=["none", "parent", "grandparent"], default="none")
    parser.add_argument("--recapitate", action="store_true", default=False)
    parser.add_argument("--ancestral_ne", type=int, default=10000)
    parser.add_argument("--mut", type=float, default=1e-7)
    args = parser.parse_args()

    # Seed parsing
    seed_part = args.host_tree.split("_seed_")[-1]
    if "_ST_" in seed_part:
        seed = int(seed_part.split("_ST_")[0])
    elif "." in seed_part:
        seed = int(seed_part.split(".")[0])
    rng = np.random.RandomState(seed)
    ts = process_subsample_iteration(
            seed,
            args.host_tree,
            args.pairs,
            n=args.sample_size,
            rng=rng,
            stratify_col=args.stratify_col,
            filter=args.filter,
            recapitate=args.recapitate,
            ancestral_ne=args.ancestral_ne
        )
    mts = msprime.sim_mutations(ts, rate=args.mut, keep=False)
    with open(args.output_file, "w") as f:
        mts.write_vcf(f, allow_position_zero = True)
    print(f"VCF written to {args.output_file}")

if __name__ == "__main__":
    main()