#!/usr/bin/env python
import random
import numpy as np
import pandas as pd
import tskit
from scipy.stats import norminvgauss
from rpy2.robjects import r, pandas2ri, numpy2ri, globalenv
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector, FloatVector
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
import logging
import argparse
import sys
import atexit
import warnings
import msprime
import time
import pyslim
from multiprocessing import Process, Pipe

warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)
rpy2_logger.setLevel(logging.ERROR)

worker_conn = None
worker_proc = None
td_call_counter = 0
TD_RESET_INTERVAL = 100000

def start_worker():
    global worker_conn, worker_proc
    parent_conn, child_conn = Pipe()
    p = Process(target=worker_loop, args=(child_conn,))
    p.daemon = True
    p.start()
    if parent_conn.poll(20.0):  # 20s timeout per user request
        msg = parent_conn.recv()
        if msg != "ready":
            raise RuntimeError("Worker failed to initialize correctly.")
    else:
        p.terminate()
        raise TimeoutError("Worker did not signal readiness in time.")
    worker_conn = parent_conn
    worker_proc = p

def stop_worker(force_kill=True):
    global worker_conn, worker_proc
    if worker_conn is not None:
        try:
            worker_conn.send(None)
        except Exception:
            pass
    if worker_proc is not None:
        worker_proc.join(timeout=1)
        if worker_proc.is_alive() and force_kill:
            print("Worker timeout — forcing kill", file=sys.stderr)
            worker_proc.terminate()
            worker_proc.join()
    worker_conn = None
    worker_proc = None

def treedist_worker(host_newick, symb_newick, timeout=4.0):
    global worker_conn, worker_proc, td_call_counter
    if worker_conn is None or worker_proc is None or not worker_proc.is_alive():
        stop_worker()
        start_worker()
    try:
        worker_conn.send((host_newick, symb_newick))
    except Exception:
        stop_worker()
        start_worker()
        return np.nan
    if worker_conn.poll(timeout):
        try:
            result = float(worker_conn.recv())
        except Exception:
            result = np.nan
    else:
        print("Timeout — restarting worker", file=sys.stderr)
        stop_worker()
        start_worker()
        return np.nan
    td_call_counter += 1
    if td_call_counter % TD_RESET_INTERVAL == 0:
        print("Resetting worker after 1000 calls", file=sys.stderr)
        stop_worker()
        start_worker()
    return result

def worker_loop(conn):
    from rpy2.robjects import r
    from rpy2.robjects.packages import importr
    try:
        ape = importr('ape')
        TreeDist = importr('TreeDist')
        conn.send("ready")
    except Exception:
        conn.send("error")
        return
    while True:
        msg = conn.recv()
        if msg is None:
            break
        host_newick, symb_newick = msg
        try:
            r.assign("host_newick", host_newick)
            r.assign("symb_newick", symb_newick)
            dist = float(
                r('TreeDist::TreeDistance(ape::read.tree(text = host_newick), ape::read.tree(text = symb_newick))[1]')[0])
            conn.send(dist)
        except Exception:
            conn.send(float('nan'))

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

def compute_pval(values, original):
    if not np.isfinite(original):
        return np.nan
    values = np.asarray(values)
    values = values[np.isfinite(values)]
    if len(values) < 990 or np.std(values) == 0:
        print(f"Too few finite permutations ({len(values)} < 990) or all identical — returning NaN", file=sys.stderr)
        return np.nan
    nig_params = norminvgauss.fit(values)
    p_nig = norminvgauss.cdf(original, *nig_params)
    return -np.log10(p_nig)

def compute_zscore(true_val, permuted_vals):
    permuted_vals = np.asarray(permuted_vals)
    permuted_vals = permuted_vals[np.isfinite(permuted_vals)]
    if len(permuted_vals) < 990 or not np.isfinite(true_val):
        print(f"Too few finite permutations < 990 or true value is nan", file=sys.stderr)
        return np.nan
    mean_perm = np.mean(permuted_vals)
    std_perm = np.std(permuted_vals)
    if std_perm == 0:
        print(f"permutations all the same", file=sys.stderr)
        return np.nan
    return (true_val - mean_perm) / std_perm

def process_subsample_iteration(
    seed, host_tree_path, symbiont_tree_path, matches_path,
    n, permutations, rng, stratify_col, filter, recapitate, ancestral_ne
):
    if recapitate:
        uncapped_host_ts = pyslim.update(tskit.load(host_tree_path))
        host_ts = pyslim.recapitate(uncapped_host_ts, ancestral_Ne=ancestral_ne)
        uncapped_symbiont_ts = pyslim.update(tskit.load(symbiont_tree_path))
        symbiont_ts = pyslim.recapitate(uncapped_symbiont_ts, ancestral_Ne=ancestral_ne)
    else:
        host_ts = tskit.load(host_tree_path)
        symbiont_ts = tskit.load(symbiont_tree_path)
    matches_df = pd.read_csv(matches_path, sep="\t")
    symbiont_to_host = matches_df.set_index("symb")["host"].to_dict()
    host_to_symbiont = matches_df.set_index("host")["symb"].to_dict()

    if filter == "none":
        selected_hosts = equal_stratified_host_sample(
            matches_df, n, rng=rng,
            group_col=stratify_col, host_col="host", fitness_col="fitness")
    else:
        selected_hosts = equal_stratified_host_sample_unique(
            matches_df, n, rng=rng, group_col=stratify_col, host_col="host", filter=filter)

    individual_node_sets = [h_ind.nodes for h_ind in host_ts.individuals() if f"h{h_ind.metadata['pedigree_id']}" in selected_hosts]
    all_host_nodes = np.concatenate(individual_node_sets)
    symbiont_names = [host_to_symbiont[name] for name in selected_hosts]
    symbiont_nodes = [s_ind.nodes[0] for s_ind in symbiont_ts.individuals() if f"s{s_ind.metadata['pedigree_id']}" in symbiont_names]
    filtered_symbiont_ts = symbiont_ts.simplify(samples=symbiont_nodes)
    node_labels_dict = {s_ind.nodes[0]: symbiont_to_host[f"s{s_ind.metadata['pedigree_id']}"] for s_ind in filtered_symbiont_ts.individuals()}
    symbiont_newick_str = str(filtered_symbiont_ts.first().as_newick(node_labels=node_labels_dict))
    globalenv["symbiont_newick"] = symbiont_newick_str
    r("shuffled_tree <- read.tree(text = symbiont_newick)")
    permuted_symb_newicks = []
    symb_labels = list(node_labels_dict.values())
    for i in range(permutations):
        shuffled_labels = StrVector(rng.permutation(symb_labels))
        r("shuffled_tree$tip.label <- {}".format(shuffled_labels.r_repr()))
        newick_i = str(r("ape::write.tree(shuffled_tree)")[0])
        permuted_symb_newicks.append(newick_i)
    filtered_host_ts = host_ts.simplify(samples=all_host_nodes)
    host_node_sets = [h_ind.nodes for h_ind in filtered_host_ts.individuals()]
    host_labels = [f"h{h_ind.metadata['pedigree_id']}" for h_ind in filtered_host_ts.individuals()]

    td_values = []
    all_shuffled_td_values = []
    start_pos = []
    for idx, tree in enumerate(filtered_host_ts.trees()):
        if tree.num_roots == 1:
            start, stop = tree.interval
            start_pos.append(start)
            avg_tmrca_matrix = np.zeros((n, n))
            for i in range(n):
                for j in range(i+1, n):
                    sample_i_1, sample_i_2 = host_node_sets[i]
                    sample_j_1, sample_j_2 = host_node_sets[j]
                    avg_tmrca_matrix[i, j] = (
                        tree.tmrca(sample_i_1, sample_j_1) +
                        tree.tmrca(sample_i_1, sample_j_2) +
                        tree.tmrca(sample_i_2, sample_j_1) +
                        tree.tmrca(sample_i_2, sample_j_2)
                    ) / 2
                    avg_tmrca_matrix[j, i] = avg_tmrca_matrix[i, j]
            r_mat_tmrca = r['matrix'](FloatVector(avg_tmrca_matrix.flatten()), nrow=n, ncol=n)
            r_df_tmrca = r['as.data.frame'](r_mat_tmrca)
            r_df_tmrca.colnames = host_labels
            r_df_tmrca.rownames = host_labels
            r.assign("D", r_df_tmrca)
            r('host_tree <- phangorn::upgma(D)')
            host_newick = str(r("ape::write.tree(host_tree)")[0])
            td = treedist_worker(host_newick, symbiont_newick_str, timeout=1)
            td_values.append(td)
            shuffled_td_values = np.zeros(permutations)
            for i in range(permutations):
                shuffled_td = treedist_worker(host_newick, permuted_symb_newicks[i], timeout=1)
                shuffled_td_values[i] = shuffled_td
            all_shuffled_td_values.append(shuffled_td_values)
    td_values = np.asarray(td_values)
    mean_td = np.nanmean(td_values)
    mean_shuffled = np.nanmean(all_shuffled_td_values, axis=0)
    median_shuffled = np.nanmedian(all_shuffled_td_values, axis=0)
    global_pval_mean = np.mean(mean_shuffled <= mean_td) if np.isfinite(mean_td) and np.isfinite(mean_shuffled).all() else np.nan
    global_pval_median = np.mean(median_shuffled <= mean_td) if np.isfinite(mean_td) and np.isfinite(median_shuffled).all() else np.nan
    global_mean_zscore = compute_zscore(mean_td, mean_shuffled)
    nan_shuffled = np.isnan(all_shuffled_td_values).sum()
    nan_tds = np.isnan(td_values).sum()
    p_vals = []
    for i in range(len(td_values)):
        log10_pval_td = compute_pval(all_shuffled_td_values[i] if i < len(all_shuffled_td_values) else [], td_values[i])
        p_vals.append(log10_pval_td)
    nan_p_vals = np.isnan(p_vals).sum()
    pval_str = ";".join([f"{val}" for val in p_vals])
    start_str = ";".join([f"{val}" for val in start_pos])
    return [
        seed, filtered_host_ts.num_trees, sum([t.num_roots > 1 for t in filtered_host_ts.trees()]),
        global_pval_mean, global_pval_median, global_mean_zscore,
        nan_shuffled, nan_tds, nan_p_vals, pval_str, start_str
    ]

atexit.register(stop_worker)

def main():
    parser = argparse.ArgumentParser(description="Unified tree permutation test pipeline.")
    parser.add_argument("--output_file", type=str, required=True)
    parser.add_argument("--symb_tree", type=str, required=True)
    parser.add_argument("--host_tree", type=str, required=True)
    parser.add_argument("--pairs", type=str, required=True)
    parser.add_argument("--sample_size", type=int, default=200)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--stratify_col", type=str, choices=["subpop", "allele"])
    parser.add_argument("--filter", type=str, choices=["none", "parent", "grandparent"], default="none")
    parser.add_argument("--recapitate", action="store_true", default=False)
    parser.add_argument("--ancestral_ne", type=int, default=10000)
    args = parser.parse_args()

    pandas2ri.activate()
    numpy2ri.activate()
    ape = importr('ape')
    phangorn = importr('phangorn')
    start_worker()
    # Seed parsing
    seed_part = args.host_tree.split("_seed_")[-1]
    if "_ST_" in seed_part:
        seed = int(seed_part.split("_ST_")[0])
    elif "." in seed_part:
        seed = int(seed_part.split(".")[0])
    rng = np.random.RandomState(seed)
    with open(args.output_file, "w") as out_f:
        if seed == 0:
            header = "seed\tnum_trees\tfiltered_trees\tglobal_pval_mean\tglobal_pval_median\tglobal_mean_zscore\tnan_shuffled\tnan_tds\tnan_p_vals\tp_vals\tstart_vals\n"
            out_f.write(header)
        result = process_subsample_iteration(
            seed,
            args.host_tree,
            args.symb_tree,
            args.pairs,
            n=args.sample_size,
            permutations=args.permutations,
            rng=rng,
            stratify_col=args.stratify_col,
            filter=args.filter,
            recapitate=args.recapitate,
            ancestral_ne=args.ancestral_ne
        )
        result_line = "\t".join([str(x) if not isinstance(x, str) else x for x in result]) + "\n"
        out_f.write(result_line)
    stop_worker()

if __name__ == "__main__":
    main()