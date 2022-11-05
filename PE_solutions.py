import itertools

class PEsolutions():
    def __init__(self, loci_candidate_dict):
        self.loci_set = loci_candidate_dict 
        self.chosen_genes = []
        self.pop_solutions = {}
        self.gene_scores = {} # dict = {gene1:[score_sol1, score_sol2,...],gene2:...}
        self.final_gene_scores = {}
        
    def generate_chosen_genes(self):
        # this function take a long time to compute the products!!!
        candidate_gene_sets = list(self.loci_set.values())
        all_combinations = list(itertools.product(*candidate_gene_sets))
        self.chosen_genes = all_combinations
        
    def generate_pop_solutions(self):
        num_solutions = len(self.chosen_genes)
        # generate a dict of loci set (num_solutions(Q) copies)
        # {1: {0: [genes], 1: [genes], ...}, ..., Q: {0: [genes], 1: [genes], ...}}
        pop_solutions = {}
        pop_sol_keys = [i for i in range(1,num_solutions+1)]
        for sol_num in pop_sol_keys:
            pop_solutions[sol_num] = self.loci_set
        self.pop_solutions = pop_solutions

    def get_single_solution(self, sol_index):
        # this method gets a solution from the multi-solutions object
        pop_solutions = self.pop_solutions
        selected_solution = pop_solutions[sol_index]

        return selected_solution

    def get_sol_chosen_genes(self, sol_index):
        # get the corresponding tuple in chosen_genes list 
        # (at position sol_index-1 as sol_index starts with 1)
        sol_chosen_genes = self.chosen_genes[sol_index-1]

        return sol_chosen_genes
    
    def update_gene_scores(self, gene, score):
        prev_score_list = self.gene_scores[gene]
        self.gene_scores[gene] = prev_score_list.append(score)
    
    def final_gene_scores(self):
        # scored_solutions is a PEsolutions object
        gene_scores_dict = self.gene_scores
        for key, value in gene_scores_dict.items():
            gene_scores_dict[key] = sum(value)/len(value)
        self.final_gene_scores = gene_scores_dict

