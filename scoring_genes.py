from network import *
import numpy as np

def gene_substitution(locus_index, itr_gene, solutions, sol_index):
    # pop solutions structure: {1: {0: [genes], 1: [genes], ...}, ..., Q: {0: [genes], 1: [genes], ...}}
    # solution structure: {locus_index0: [genes]}
    # get the all chosen genes in all the solutions
    sol_chosen_genes = solutions.chosen_genes[sol_index]
    # get the chosen gene at the specified index of the specified solution
    locus_chosen_gene = sol_chosen_genes[locus_index]
    chosen_replaced = list(map(lambda x: x.replace(locus_chosen_gene, itr_gene), sol_chosen_genes))

    return chosen_replaced

def compute_density(chosen_replaced, network):
    # density formula: edge counts
    # count edges connected among genes in chosen_replaced
    # get the connections from network (direct neighbors)
    gene_pairs = [(a, b) for idx, a in enumerate(chosen_replaced) for b in chosen_replaced[idx + 1:]]
    density = 0
    for gene_pair in gene_pairs:
        try:
            edge_count = network.find_edge(gene_pair[0], gene_pair[1])
        except: 
            edge_count = 0
        density = density + edge_count

    return density

def empty_locus_case(locus_index, chosen_replaced, network):
    # remove the indicated locus completely from the solution
    del chosen_replaced[locus_index]
    gene_pairs = [(a, b) for idx, a in enumerate(chosen_replaced) for b in chosen_replaced[idx + 1:]]
    empty_locus_density = 0
    for gene_pair in gene_pairs:
        try:
            edge_count = network.find_edge(gene_pair[0], gene_pair[1])
        except: 
            edge_count = 0
        empty_locus_density  = empty_locus_density + edge_count

    return empty_locus_density 

def score_gene(locus_index, chosen_replaced, network):
    density = compute_density(chosen_replaced, network)
    empty_locus_density = empty_locus_case(locus_index, chosen_replaced, network)
    gene_score = np.abs(density - empty_locus_density)

    return gene_score

