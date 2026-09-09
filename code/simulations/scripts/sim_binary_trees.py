import msprime
import tskit
import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import argparse

def generate_uncorrelated_binary_trees(num_trees, num_tips, mu=1e-8, rec=1e-8, Ne=10000, model='hudson'):
    
    newick_trees = []

    rand_seeds = np.random.randint(1, 2**32, num_trees)
    
    for _ in range(num_trees):
        ts = msprime.simulate(sample_size=num_tips, Ne=Ne, mutation_rate=mu, random_seed=rand_seeds[_], model=model)
        
        # only one tree in the sequence
        for tree in ts.trees():
            newick_trees.append(tree.newick())

    return newick_trees

def parse_arguments():
    parser = argparse.ArgumentParser(description="Simulate binary trees and save them to a file.")
    parser.add_argument("-o", "--output_file", type=str, required=True, help="Path to the output file for saving the binary trees (pickle format).")
    parser.add_argument("--num_trees", type=int, required=True, help="Number of trees to simulate.")
    parser.add_argument("--num_leaves", type=int, required=True, help="Number of leaves (tips) per tree.")
    return parser.parse_args()

def main():
    args = parse_arguments()
    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    newick_trees = generate_uncorrelated_binary_trees(args.num_trees, args.num_leaves)
    with open(args.output_file, 'wb') as f:
        pickle.dump(newick_trees, f)

    print(f"Newick trees saved to '{args.output_file}'")

if __name__ == "__main__":
    main()