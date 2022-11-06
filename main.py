from scoring_genes import *
from utilities_module3 import * 
from PF_solutions import *

# By default: if no 2 column txt file is provided
# get loci candidate genes from the gmt file for FA
loci_candidate_dict = get_loci_candidate_genes()
print(loci_candidate_dict)
# create Prix fixe solutions object
solutions = PFsolutions(loci_candidate_dict)
#print(solutions.gene_scores)
# create network object
network = Network("STRING_network.txt")
network.get_network_interactions_degrees()
# randomly generate chosen genes for each solution (one gene per locus)
solutions.generate_chosen_genes()
# get number of solutions = number of all possible combinations of chosen genes
num_sols = len(solutions.chosen_genes)
print("Number of generated solutions: ",num_sols)
# get the generated chosen genes
chosen_genes = solutions.chosen_genes
for sol_index in range(1,num_sols):
    for locus in solutions.loci_set.keys():
        locus_gene_candidates = solutions.loci_set[locus] 
        locus_chosen_gene = solutions.get_sol_chosen_genes(sol_index)
        locus_itr_genes = list(set(locus_gene_candidates) - set(locus_chosen_gene))
        for itr_gene in locus_itr_genes:
            chosen_replaced = gene_substitution(locus, itr_gene, solutions, sol_index)
            gene_score = score_gene(locus, chosen_replaced, network)
            solutions.update_gene_scores(itr_gene, gene_score)
for key, value in solutions.gene_scores.items():
        solutions.gene_scores[key] = sum(value)/len(value)
final_gene_scores = solutions.gene_scores
print(final_gene_scores)




