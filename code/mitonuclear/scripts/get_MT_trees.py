import numpy as np
import pickle
import argparse
from rpy2.robjects import r, pandas2ri, numpy2ri, globalenv
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector

# Activate pandas <-> R dataframes conversion
pandas2ri.activate()
numpy2ri.activate()

# Import necessary R package
ape = importr('ape')
phangorn = importr('phangorn')

def main():
    # Sets random generator up with seed 42
    rng = np.random.default_rng(42)

    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Create mitochondrial trees.")

    # Define required positional arguments
    parser.add_argument("-mt_msa", type=str, help="Path to the mitochondrial fasta")
    parser.add_argument("-true_tree", type=str, help="Path to the output file for the true tree")
    parser.add_argument("-shuffled_trees", type=str, help="Path to the output pickle file for shuffled trees")
    parser.add_argument("--permutations", type=int, required=True, help="Number of permutations per tree comparison")

    # Parse the arguments
    args = parser.parse_args()

    # Assign arguments to variables
    mt_msa = args.mt_msa
    true_tree = args.true_tree
    shuffled_trees = args.shuffled_trees
    permutations = args.permutations

    # Calculate distance matrix using the K80 model
    globalenv['mt_msa'] = mt_msa
    r('mito_alignment <- read.FASTA(mt_msa, type = "DNA")')
    r('mt_dm <- as.matrix(dist.dna(mito_alignment, model="K80"))')
    r('mt_upgma <- upgma(mt_dm)')
    mt_upgma = r['mt_upgma']

    # Write the true phylogenetic tree to a Newick file
    ape.write_tree(mt_upgma, file=true_tree)

    # Generate permuted trees and collect as Newick strings
    col_names = r('colnames(mt_dm)')
    shuffled_tree_list = []

    for i in range(permutations):
        shuffled_labels = StrVector(rng.permutation(col_names))
        r("mt_upgma$tip.label <- {}".format(shuffled_labels.r_repr()))
        newick_str = r('write.tree(mt_upgma)')
        shuffled_tree_list.append(str(newick_str[0]))

    # Save shuffled trees as a pickle
    with open(shuffled_trees, "wb") as f:
        pickle.dump(shuffled_tree_list, f)

if __name__ == "__main__":
    main()
