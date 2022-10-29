from multiprocessing.spawn import is_forking
import pandas as pd
import numpy as np
import networkx as nx
import ipycytoscape as cy
from ipycytoscape import *
from collections import Counter

# class to keep track of network's activities, check network's properties, and find a subnetwork for a given gene set
class Network():
    def __init__(self, network_name, flip="no", sep = "\t"):
        # if a network data frame is provided, assign it to network var.
        if isinstance(network_name, pd.DataFrame) == True:
            network = network_name
        # read in a network file if file name is provided
        elif isinstance(network_name, str) == True:
            # network in 3-columns format
            network = pd.read_csv(network_name, sep = sep, names = ["gene1","gene2", "weight"],index_col=False)
        
        # flip the weights so a higher weight = a stronger biological interaction
        if flip == "yes":
            network.weight = 1-network.weight
        elif flip == "no":
            pass
        # initialize attributes: network data frame, list of visited nodes
        self.network_df = network
        self.network_interactions = {}
        self.network_degrees = {}

    def get_network_interactions_degrees(self):
        # this function gets and returns all interactions and degree information of the network (dicts)
        # convert network_df to list of tuples
        network_df = self.network_df
        network_list = list(network_df.itertuples(index=False, name=None))
        network_interactions = {} 
        node_list = []
        for element in network_list:
            strs = list(filter(lambda x : type(x) == str, element))
            weights = list(filter(lambda x : type(x) == float, element))[0]
            network_interactions[tuple(sorted(strs))] = weights
            try:
                str_ = strs[0]
            except:
                print(element)
                pass 
            node_list.append(str_)
        network_degrees = dict(Counter(node_list))
        self.network_interactions = network_interactions
        self.network_degrees = network_degrees
        
    def get_interaction(self, gene1, gene2):
        # this function get weight of the edge connected given genes
        gene_pair = tuple(sorted([gene1, gene2]))
        edge_weight = self.network_interactions[gene_pair] 
        return edge_weight 
        