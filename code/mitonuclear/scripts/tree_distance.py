import sys
import numpy as np
import pickle
import os
import time
import tskit
import argparse
from rpy2.robjects import r, pandas2ri, numpy2ri
from rpy2.robjects.vectors import FloatVector
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
import logging
from scipy.stats import norminvgauss
from multiprocessing import Process, Pipe
import atexit

pandas2ri.activate()
numpy2ri.activate()
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
    if parent_conn.poll(10.0):
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

def treedist_worker(host_newick, symb_newick, timeout=1.0):
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
        print("Resetting worker after 100000 calls", file=sys.stderr)
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
    except Exception as e:
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
            
def compute_pval(values, original):
    if not np.isfinite(original): return np.nan
    values = np.asarray(values)
    values = values[np.isfinite(values)]
    if len(values) < 990:
        print(f"Too few finite permutations ({len(values)} < 990) — returning NaN", file=sys.stderr)
        return np.nan
    nig_params = norminvgauss.fit(values)
    pval = norminvgauss.cdf(original, *nig_params)
    return -np.log10(pval)

atexit.register(stop_worker)

def main():
    parser = argparse.ArgumentParser(description="Measure cophylogeny between population and mitochondrial genome.")
    parser.add_argument("--sample_list_path", type=str, required=True)
    parser.add_argument("--tree_sequence_file", type=str, required=True)
    parser.add_argument("--mito_tree", type=str, required=True)
    parser.add_argument("--shuffled_trees", type=str, required=True)
    parser.add_argument("--output_file", type=str, required=True)
    parser.add_argument("--total_para", type=int, required=True)
    parser.add_argument("--parallele_idx", type=int, required=True)
    parser.add_argument("--permutations", type=int, required=True)
    args = parser.parse_args()

    ts = tskit.load(args.tree_sequence_file)
    with open(args.sample_list_path, 'r') as f: idsToExtract = [line.strip() for line in f]
    with open(args.mito_tree, 'r') as f: mito_newick = f.read().strip()
    with open(args.shuffled_trees, "rb") as f: shuffled_newick_list = pickle.load(f)

    num_samples = ts.num_samples
    num_individuals = num_samples // 2
    if num_individuals != len(idsToExtract):
        print("Mismatch between individuals and sample list")
        sys.exit(1)

    total_trees = ts.num_trees

    partitions = np.linspace(0, total_trees, args.total_para + 1, dtype=int)
    first_tree_idx = 0 if args.parallele_idx == 0 else partitions[args.parallele_idx] + 1
    last_tree_idx = partitions[args.parallele_idx + 1]

    output_dir = os.path.dirname(args.output_file)
    os.makedirs(output_dir, exist_ok=True)
    start_worker()

    with open(args.output_file, 'w') as out_file:
        if first_tree_idx == 0:
            out_file.write(f"start_pos\tend_pos\tspan\ttd\tlog10_pval\tnan_shuffled\tshuffled_td_values\n")

        for idx, tree in enumerate(ts.trees()):
            if idx < first_tree_idx or idx > last_tree_idx:
                continue
            print(f"working on {idx}")
            start, stop = tree.interval
            out_file.write(f"{start}\t{stop}\t{stop-start}\t")

            avg_tmrca_matrix = np.zeros((num_individuals, num_individuals))
            for i in range(num_individuals):
                for j in range(i + 1, num_individuals):
                    s1, s2 = 2 * i, 2 * i + 1
                    t1, t2 = 2 * j, 2 * j + 1
                    avg_tmrca_matrix[i, j] = (tree.tmrca(s1, t1) + tree.tmrca(s1, t2) + tree.tmrca(s2, t1) + tree.tmrca(s2, t2)) / 2
                    avg_tmrca_matrix[j, i] = avg_tmrca_matrix[i, j]

            r_mat = r['matrix'](FloatVector(avg_tmrca_matrix.flatten()), nrow=num_individuals, ncol=num_individuals)
            r_df = r['as.data.frame'](r_mat)
            r_df.colnames = idsToExtract
            r_df.rownames = idsToExtract
            r.assign("D", r_df)
            r('host_tree <- phangorn::upgma(D)')
            host_newick = str(r("ape::write.tree(host_tree)")[0])
            td = treedist_worker(host_newick, mito_newick)

            shuffled_td_values = []
            for i in range(args.permutations):
                shuffled_td = treedist_worker(host_newick, shuffled_newick_list[i])
                shuffled_td_values.append(shuffled_td)
            log10_pval = compute_pval(shuffled_td_values, td)
            nan_count = np.isnan(shuffled_td_values).sum()
            td_strs = ";".join([f"{v:.8f}" for v in shuffled_td_values])
            out_file.write(f"{td}\t{log10_pval}\t{nan_count}\t{td_strs}\n")

    stop_worker()
    print(f"Processed {total_trees} trees in {time.time() - start_time:.2f} seconds.")

if __name__ == "__main__":
    start_time = time.time()
    main()