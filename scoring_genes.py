import numpy as np
import pandas as pd
import random 
from utilities_module3 import *
from shortest_path_module2 import *

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

def generate_random_solsubnet(loci_candidate_dict):
    random_loci_candidates = {}
    for key in loci_candidate_dict.keys():
        locus_candidates = loci_candidate_dict[key]
        num_list = [2,3,4]
        ramdom_num = random.choice(num_list)
        random_locus_genes = random.choices(locus_candidates, k=ramdom_num)
        random_loci_candidates[key] = random_locus_genes
    return random_loci_candidates

# default random loci candidates
random_loci_candidates = {0: ['LCMT1', 'SULT1A2', 'CTB-134H23.2', 'SULT1A2'], 1: ['C11orf91', 'HIPK3', 'NAT10'], 2: ['AXIN2', 'CLTC', 'CTD-2535L24.2', 'RP11-15E18.4'], 3: ['NCBP1', 'KRT8P11', 'HIATL2', 'XPA'], 4: ['DEF8', 'DEF8'], 5: ['KISS1', 'OPTC', 'KISS1', 'CNTN2'], 6: ['TAMM41', 'RAF1'], 7: ['KIF6', 'CLPSL2', 'C6orf222'], 8: ['KBTBD6', 'POSTN', 'MRPS31'], 9: ['CCP110', 'C16orf45'], 10: ['AC087477.1', 'VPS33B', 'TICRR', 'AC016251.1'], 11: ['ZNF500', 'TFAP4', 'TMEM186', 'DNASE1']}

def gene_candidate_subnetwork(network, candidate_gene, locus_index, loci_genes_dict = random_loci_candidates):
    '''
    Replace 'chosen' gene at locus i with given candidate gene.
    Compute new subnetwork using all other chosen genes and the given candidate gene
    Return recomputed subnetwork
    '''
    # replace 'chosen' gene at locus i with given candidate gene.
    chosen_genes = list(loci_genes_dict.values())
    new_query_genes = list(map(lambda x: x.replace(loci_genes_dict[locus_index], candidate_gene), chosen_genes))
    # compute new subnetwork using all other chosen genes and the given candidate gene
    recomputed_subnetwork = get_final_subnetwork(network, query_nodes = new_query_genes)
    return recomputed_subnetwork

def empty_locus_case(loci_genes_dict, locus_index):
    '''
    Remove the chosen gene at the given locus index from the list of chosen genes.
    Compute new subnetwork using all other choden genes.
    Return empty locus subnetwork
    '''
    gene_to_remove = locus_index[locus_index]
    loci_genes = list(loci_genes_dict.values())
    loci_genes.remove(gene_to_remove)
    empty_locus_genes = loci_genes
    return empty_locus_genes 

def score_gene(recomputed_subnetwork, empty_locus_subnetwork):
    '''
    Calculate gene score: score = |new density - empty locus density| 
    recomputed_subnetwork (a list of paths: tuple(gene1,gene2,weight)) is subnetwork built on replacing the chosen gene i with a candidate geneat locus i
    empty_locus_subnetwork (a list of paths: tuple(gene1,gene2,weight)) is subnetwork with locus i removed 
    Return gene score
    '''
    # compute new subnetwork density
    recomputed_density = calculate_density(recomputed_subnetwork)
    # compute subnetwork density with an empty locus
    empty_locus_density = calculate_density(empty_locus_subnetwork)
    gene_score = np.abs(recomputed_density - empty_locus_density)
    return gene_score