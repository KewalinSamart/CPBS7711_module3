import numpy as np
import pandas as pd
import random 
from utilities_module3 import *

# move this to PEsoultions class
def get_chosen_candidate_genes(gmt_filename, delimiter='\t'):
    myfile = open(gmt_filename)
    loci_genes_dict = {}
    loci_candidate_dict = {}
    loci_index = 0
    for line in myfile :
        locus_gene = line.split(delimiter)[1].split(" ")[2]
        loci_genes_dict[loci_index] = locus_gene
        candidate_genes = line.split("\t")[2:]
        candidate_genes = [elem.replace('\n', '') for elem in candidate_genes]
        loci_candidate_dict[loci_index] =  candidate_genes
        loci_index = loci_index + 1
        
    return loci_genes_dict, loci_candidate_dict

def gene_substitution(itr_gene, chosen_gene):
    pass

def compute_density(itr_gene, solution):
    pass

def empty_locus_case(locus_index, solution):
    pass

def score_gene(scored_solutions, gene):
    pass 