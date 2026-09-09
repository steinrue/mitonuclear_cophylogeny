import sys
import numpy as np
import random
import os
import time
import tskit
import io
from rpy2.robjects import r, pandas2ri, numpy2ri, globalenv
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector, FloatVector
from tqdm import tqdm
import argparse
import pickle


def reorder_tree_vector(ref_labels, tree2_labels, t2_vector):
    t1_argsort = np.argsort(np.array(ref_labels))
    t2_argsort = np.argsort(np.array(tree2_labels))
    return reorder_tree_vector_internal(t1_argsort, t2_argsort, t2_vector)

def reorder_tree_vector_internal(t1_idxs, t2_idxs, t2_vector):
    n = len(t1_idxs)
    t1_pos_in_t2 = t2_idxs[np.argsort(t1_idxs)]

    vector_idxs = np.zeros_like(t2_vector, dtype=int)
    vector_is = np.zeros_like(vector_idxs, dtype=int)
    vector_js = np.zeros_like(vector_idxs, dtype=int)
    overall_idx = 0
    full_js_range = np.arange(1,t1_pos_in_t2.shape[0])
    for i in range(t1_pos_in_t2.shape[0]):
        num_js = full_js_range.shape[0]-i
        vector_is[overall_idx:overall_idx+num_js] = t1_pos_in_t2[i]
        vector_js[overall_idx:overall_idx+num_js] = t1_pos_in_t2[full_js_range[i:]]
        overall_idx += num_js
    super_true_is = np.minimum(vector_is, vector_js)
    super_true_js = np.maximum(vector_is, vector_js)
    vector_idxs = super_true_js - (super_true_is+1) + super_true_is * n - super_true_is * (super_true_is+1) // 2
    return t2_vector[vector_idxs]

def calc_single_kc_distance(tip_labels1, edge_matrix1, tip_labels2, edge_matrix2):
    tv1 = compute_tree_vector(tip_labels1, edge_matrix1)
    tv2 = compute_tree_vector(tip_labels2, edge_matrix2)
    tv2_sorted = reorder_tree_vector(tip_labels1, tip_labels2, tv2)
    return np.sqrt(np.sum((tv1 - tv2_sorted)**2))


def calc_one_to_many_kc_distance(ref_tip_labels, ref_edge_matrix, tip_labels_list, edge_matrix_list):
    assert len(tip_labels_list) == len(edge_matrix_list)
    tree_dists = []
    tv1 = compute_tree_vector(ref_tip_labels, ref_edge_matrix)
    for t_i in range(len(tip_labels_list)):
        tv_temp = compute_tree_vector(tip_labels_list[t_i], edge_matrix_list[t_i])
        tv_reordered = reorder_tree_vector(ref_tip_labels, tip_labels_list[t_i], tv_temp)
        tree_dists.append(np.sqrt(np.sum((tv1 - tv_reordered)**2)))
    return tree_dists

def compute_tree_vector(tip_labels, edge_matrix):
    n_tips = len(tip_labels)
    vector = np.zeros(((n_tips*(n_tips-1))//2,),dtype=int)


    tips_sorted = np.argsort(edge_matrix[:, 1])
    ancestor_grid = np.zeros((1, n_tips), dtype=int)
    ancestor_grid[0, :] = np.arange(1, n_tips+1)
    while np.any(ancestor_grid[-1, :] != n_tips+1):
        tip_pos = np.searchsorted(edge_matrix[:, 1][tips_sorted], ancestor_grid[-1, :])
        new_indices = tips_sorted[tip_pos]
        ancestor_grid = np.vstack((ancestor_grid, edge_matrix[new_indices, 0]))

    non_tip_nodes = np.arange(n_tips+2, np.max(edge_matrix)+1)
    non_tip_tf = np.zeros((non_tip_nodes.shape[0], ancestor_grid.shape[1]), dtype=bool)
    for nt_i, nt_node in enumerate(non_tip_nodes):
        non_tip_tf[nt_i, np.where(ancestor_grid[1:]==nt_node)[1]] = True
    non_tip_tf = non_tip_tf[::-1, :]
    non_tip_dists = np.zeros_like(non_tip_nodes)
    internal_edges = edge_matrix - (n_tips+1)
    internal_edges = internal_edges[internal_edges[:, 1] > 0, :]
    iter_i = 0
    while internal_edges.shape[0]:
        next_coords = internal_edges[internal_edges[:, 0]==0, 1]
        if iter_i:
            non_tip_dists[next_coords+iter_i-1] = non_tip_dists[iter_i - 1] + 1
        else:
            non_tip_dists[next_coords+iter_i-1] = 1
        internal_edges -= 1
        internal_edges = internal_edges[internal_edges[:, 1] > 0, :]
        iter_i += 1

    non_tip_dists = non_tip_dists[::-1]
    non_tip_dists = np.hstack((non_tip_dists, [0]))

    giant_tf_i = np.zeros((non_tip_nodes.shape[0], vector.shape[0]), dtype=bool)
    giant_tf_j = np.zeros((non_tip_nodes.shape[0], vector.shape[0]), dtype=bool)

    rc2 = 0
    for i in np.arange(ancestor_grid.shape[1]-1):
        giant_tf_i[:, rc2:rc2+ancestor_grid.shape[1]-i-1] = non_tip_tf[:, i][:, np.newaxis]
        giant_tf_j[:, rc2:rc2+ancestor_grid.shape[1]-i-1] = non_tip_tf[:, np.arange(i+1, ancestor_grid.shape[1])]
        rc2 += ancestor_grid.shape[1]-i-1

    giant_tf_ij = giant_tf_i & giant_tf_j
    giant_tf_ij = np.vstack((giant_tf_ij, np.zeros(giant_tf_ij.shape[1], dtype=bool)+True))
    vector2 = non_tip_dists[np.argmax(giant_tf_ij, axis=0)]
    return vector2


# Activate pandas <-> R dataframes conversion
pandas2ri.activate()
numpy2ri.activate()

# Import necessary R packages
ape = importr('ape')
phangorn = importr('phangorn')
TreeDist = importr('TreeDist')
Quartet = importr('Quartet')
readr = importr('readr')
arrowr = importr('arrow')

def main():
    rng = np.random.default_rng(42)
    parser = argparse.ArgumentParser(description="Compare Newick Trees.")

    # Define required positional arguments
    parser.add_argument("-tree_list", type=str, required=True, help="Path to the pickle file of newick trees")
    parser.add_argument("-focal_tree_idx", type=int, required=True, help="Index of the focal tree")
    parser.add_argument("-num_perms", type=int, required=True, help="Number of permutations")
    parser.add_argument("-o", type=str, required=True, help="Output file for focal tree results.")

    # Parse the arguments
    args = parser.parse_args()

    # Assign arguments to variables
    tree_list_path = args.tree_list
    focal_idx = args.focal_tree_idx
    num_perms = args.num_perms
    output_file = args.o

    # Load the Newick trees list from the pickle file
    with open(tree_list_path, 'rb') as f:
        tree_list = pickle.load(f)

    # Read the two Newick strings into an R phylo object
    globalenv["focal_tree"] = ape.read_tree(text=tree_list[focal_idx])
    globalenv["test_tree"] = ape.read_tree(text=tree_list[focal_idx + 100])

    # Generate new tip labels in the format "n<original_label>"
    r("focal_tree$tip.label <- paste0('n', focal_tree$tip.label)")
    r("test_tree$tip.label <- paste0('n', test_tree$tip.label)")
    focal_tree = r("focal_tree")
    test_tree = r("test_tree")

    # Extract tip labels from the phylo object in R
    tip_labels = r("test_tree$tip.label")

    # Start timing the process
    start_time = time.time()

    # Ensure the directory for the output file exists
    output_dir = os.path.dirname(output_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Initialize a list to store all shuffled_qd values
    labelset = []
    shuffled_qd_values = []
    shuffled_td_values = []
    shuffled_rf_values = []
    shuffled_kc_values = []

    # Open the main output file for writing
    with open(output_file, 'w') as out_file:
        if focal_idx == 0:
            out_file.write("td\tkc\trf\tqd\n")

        for _ in tqdm(range(num_perms), desc="Processing Trees"):
            # Shuffle test tree labels
            shuffled_labels = StrVector(rng.permutation(tip_labels))
            labelset.append(shuffled_labels)
            r("perm_tree <- test_tree")
            r("perm_tree$tip.label <- {}".format(shuffled_labels.r_repr()))
            permuted_tree = r("perm_tree")

            # Calculate qd distances
            statuses = r('QuartetStatus(focal_tree, perm_tree)')
            shuffled_qd = float(Quartet.QuartetDivergence(statuses, similarity = False)[0])
            shuffled_qd_values.append(shuffled_qd)

            # Calculate td Distances
            shuffled_td = float(TreeDist.TreeDistance(focal_tree, permuted_tree)[0])
            shuffled_td_values.append(shuffled_td)

            # Calculate Robinson Folds Distance
            shuffled_rf = float(TreeDist.RobinsonFoulds(focal_tree, permuted_tree)[0])
            shuffled_rf_values.append(shuffled_rf)

        # Generate shuffled kc distances
        focal_tree_tip_labels = focal_tree.rx2("tip.label")
        focal_tree_edge_matrix = np.array(focal_tree.rx2("edge"))
        tv_focal = compute_tree_vector(focal_tree_tip_labels, focal_tree_edge_matrix)
        test_tree_tip_labels = test_tree.rx2("tip.label")
        test_tree_edge_matrix = np.array(test_tree.rx2("edge"))
        tv_test = compute_tree_vector(test_tree_tip_labels, test_tree_edge_matrix)
        for shffl in labelset:
            tv_shuffled = reorder_tree_vector(test_tree_tip_labels, shffl, tv_focal)
            shuffled_kc_values.append(np.sqrt(np.sum((tv_test-tv_shuffled)**2)))

        # Sort the Values and turn them into strings
        shuffled_td_str = ";".join([f"{val:.8f}" for val in np.sort(shuffled_td_values)])
        shuffled_qd_str = ";".join([f"{val:.8f}" for val in np.sort(shuffled_qd_values)])
        shuffled_rf_str = ";".join([f"{val:.8f}" for val in np.sort(shuffled_rf_values)])
        shuffled_kc_str = ";".join([f"{val:.8f}" for val in np.sort(shuffled_kc_values)])
        # Write the results to the output file
        out_file.write(f"{shuffled_td_str}\t{shuffled_kc_str}\t{shuffled_rf_str}\t{shuffled_qd_str}\n")

    # End timing the process
    end_time = time.time()
    progress_time = end_time - start_time
    print(f"Processed {len(tree_list)} trees in {progress_time:.2f} seconds.")

if __name__ == "__main__":
    main()