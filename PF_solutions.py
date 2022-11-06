import itertools

class PFsolutions():
    def __init__(self, loci_candidate_dict):
        self.loci_set = loci_candidate_dict 
        self.chosen_genes = []
        # combine all genes from different loci
        candidate_genes = list(itertools.chain.from_iterable(loci_candidate_dict.values()))
        self.gene_scores = {gene: [0] for gene in candidate_genes} # dict = {gene1:[score_sol1, score_sol2,...],gene2:...}
        self.final_gene_scores = self.gene_scores
        
    def generate_chosen_genes(self):
        # this function take a long time to compute the products!!!
        candidate_gene_sets = list(self.loci_set.values())
        all_combinations = list(itertools.product(*candidate_gene_sets))
        self.chosen_genes = all_combinations

    def get_sol_chosen_genes(self, sol_index):
        # get the corresponding tuple in chosen_genes list 
        # (at position sol_index-1 as sol_index starts with 1)
        sol_chosen_genes = self.chosen_genes[sol_index-1]

        return list(sol_chosen_genes)
    
    def update_gene_scores(self, gene, score):
        #print(self.gene_scores)
        if gene in self.gene_scores.keys():
            self.gene_scores[gene].append(score)
        else:
            self.gene_scores[gene] = [score]
        #print(self.gene_scores)
    
    def final_gene_scores(self):
        # scored_solutions is a PEsolutions object
        gene_scores_dict = self.gene_scores
        for key, value in gene_scores_dict.items():
            gene_scores_dict[key] = sum(value)/len(value)
        self.final_gene_scores = gene_scores_dict

