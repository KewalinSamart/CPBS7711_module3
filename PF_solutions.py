import itertools
import random
import pandas as pd 
import networkx as nx 
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, text
import numpy as np


class PFsolutions():
    def __init__(self, loci_candidate_dict, annotated_candidate_dict, chosen_genes = []):
        self.loci_set = loci_candidate_dict 
        self.annotated_candidate_dict = annotated_candidate_dict
        self.chosen_genes = chosen_genes # empty list by default
        # combine all genes from different loci
        candidate_genes = list(itertools.chain.from_iterable(loci_candidate_dict.values()))
        self.gene_scores = {gene: [0] for gene in candidate_genes} # dict = {gene1:[score_sol1, score_sol2,...],gene2:...}
        self.final_gene_scores = self.gene_scores
        self.final_sol_df = pd.DataFrame()
        
    def generate_chosen_genes(self,iteration=5000):
        # this function take a long time to compute the products!!!
        chosen_genes = []
        for itr in range(iteration):
            random_genes = []
            for gene_candidates in self.loci_set.values():
                random_gene = random.choice(gene_candidates)
                random_genes.append(random_gene)
            tp_random_genes = tuple(random_genes)
            chosen_genes.append(tp_random_genes)
        self.chosen_genes = chosen_genes

    def get_sol_chosen_genes(self, sol_index):
        # get the corresponding tuple in chosen_genes list 
        # (at position sol_index-1 as sol_index starts with 1)
        sol_chosen_genes = self.chosen_genes[sol_index-1]

        return list(sol_chosen_genes)
    
    def update_gene_scores(self, gene, score):
        if gene in self.gene_scores.keys():
            self.gene_scores[gene].append(score)
        else:
            self.gene_scores[gene] = [score]
    
    def compute_final_scores(self):
        # scored_solutions is a PEsolutions object
        gene_scores_dict = self.gene_scores
        for key, value in gene_scores_dict.items():
            gene_scores_dict[key] = sum(value)/len(value)
        self.final_gene_scores = gene_scores_dict

    def finalize_final_sol(self):
        final_scores_df = pd.DataFrame(self.final_gene_scores.items(), columns=['gene','score'])
        loci_candidates_df = pd.DataFrame(self.annotated_candidate_dict.items(), columns=['gene','locus'])
        self.final_sol_df = pd.merge(final_scores_df, loci_candidates_df, on ='gene')
    
    def output_final_sol(self, output_dir='example_result/final_solution.txt'):
        self.final_sol_df.to_csv(output_dir, header=['gene','score','locus'], index=None, sep=' ', mode='a')
        print("The final solution was saved at ", output_dir)

    def visualize_scored_sol(self):
        final_sol = self.final_sol_df
        num_loci = len(self.loci_set)
        genes = list(self.gene_scores.keys())
         
        pass 

