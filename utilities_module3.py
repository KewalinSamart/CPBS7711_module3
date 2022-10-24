import pandas as pd
from collections import Counter

def get_genes_gmt(filename, delimiter = "\t"):
    subset_info = pd.read_csv(filename ,sep=delimiter, header=None, engine='python')
    subgenes_info = subset_info[0].str.split(pat = delimiter, n=-1, expand=False)
    subgenes_info_list = subgenes_info.tolist()
    subgenes_list = []
    for i in range(len(subgenes_info_list)):
        representative_gene = subgenes_info_list[i][1].split()[2]
        representative_gene = representative_gene.replace(",","")
        subgenes_list.append(representative_gene)
    return sorted(subgenes_list)

def intersection(list1, list2):
    intersection_list = [element for element in list1 if element in list2]
    return intersection_list

def paths_to_df(indiv_shortest_paths):
    indiv_path_df = pd.DataFrame(indiv_shortest_paths, columns =['gene1', 'gene2', 'weight'])
    return indiv_path_df

def get_network_interactions_degrees(network_df):
    # convert network_df to list of tuples
    network_list = list(network_df.itertuples(index=False, name=None))
    network_interactions = {} 
    node_list = []
    for element in network_list:
        strs = list(filter(lambda x : type(x) == str, element))
        weights = list(filter(lambda x : type(x) == float, element))[0]
        network_interactions[tuple(sorted(strs))] = weights
        node_list.append(strs[0])
    network_degrees = dict(Counter(node_list))
    return network_interactions, network_degrees

def calculate_density(subnetwork_paths):
    # this function computes subnetwork density 
    subnetwork_df = paths_to_df(subnetwork_paths)
    subnetwork_df['interaction_strength'] = 1 - subnetwork_df.weight
    interaction_strength = list(subnetwork_df.interaction_strength)
    sum_weight = sum(interaction_strength)
    # get number of vertices
    num_vertices = len(list(set(list(subnetwork_df.gene1)+list(subnetwork_df.gene2))))
    density = 2*sum_weight/(num_vertices*(num_vertices-1))
    return density