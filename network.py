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
        self.visited_nodes = []

    def update_visited_nodes(self, visited_node):
        # this method updates visited nodes in the network
        self.visited_nodes.append(visited_node)

    def clear_visited_nodes(self):
        # this method clears out all the visited nodes recorded
        self.visited_nodes = []
    
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
        return network_interactions, network_degrees
    
    def get_interactions(self, path_list):
        # this function gets and returns interactions on the given path
        interacting_pairs = []
        for pos in range(len(path_list)-1):
            gene_pair = tuple(sorted(list((path_list[pos],path_list[pos+1]))))
            weight = self.get_network_interactions_degrees()[0][gene_pair]
            interacting_pairs.append(gene_pair+(weight,))
        return interacting_pairs

    def get_neighbors(self, node):
        # this functions gets and returns neighbors of the given node
        interactions = self.get_network_interactions_degrees()[0]
        neighbor_dict = {}
        degrees = self.get_network_interactions_degrees()[1]
        one_degree = [k for k,v in degrees.items() if int(v) == 1]
        for tub, weight in interactions.items():
            if node in tub:
                index = tub.index(node)
                neighbor = tub[abs(1-index)] 
                #if (neighbor not in one_degree) and (neighbor != dest_node):
                if neighbor not in self.visited_nodes: 
                    neighbor_dict[neighbor] = weight
                #else:
                #    self.update_visited_nodes(neighbor)
        return neighbor_dict
        
    def visualize_network(self):
        # this function visualizes network object on jupyter notebook interface
        edge_data = {'source': list(self.network_df.gene1),
                     'target': list(self.network_df.gene2),
                     'weight': list(self.network_df.weight)}
        link_df = pd.DataFrame.from_dict(edge_data)
        node_data = {'id': list(set(self.network_df.gene1))}

        node_df = pd.DataFrame.from_dict(node_data)
        G = nx.from_pandas_edgelist(link_df, edge_attr=True)
        node_df = pd.DataFrame(G.degree(), columns=['id', 'degree'])
        cytoscapeobj = cy.CytoscapeWidget()
        cytoscapeobj.graph.add_graph_from_networkx(G)

        cytoscapeobj.set_style(
        [
            {
                 'selector': 'node',
                 'style': {
                     'font-family': 'helvetica',
                     'font-size': '20px',
                     'label': 'data(id)'
                }
            },   
            { # 'BRCA2', 'ERCC4', 'FANCA', 'FANCC', 
            # 'FANCD2', 'FANCE', 'FANCF', 'FANCI', 'PALB2', 'RAD51C', 'SLX4', 'UBE2T'
                'selector': 'node[id = "BRCA2"]',
                'style': {
                'font-family': 'helvetica',
                'font-size': '20px',
                'label': 'data(id)',
                'background-color': 'red'}
            },
            {
                'selector': 'node[id = "ERCC4"]',
                'style': {
                'font-family': 'helvetica',
                'font-size': '20px',
                'label': 'data(id)',
                'background-color': 'red'}
            },
            {
                'selector': 'node[id = "FANCA"]',
                'style': {
                'font-family': 'helvetica',
                'font-size': '20px',
                'label': 'data(id)',
                'background-color': 'red'}
            },
            {
                'selector': 'node[id = "FANCC"]',
                'style': {
                'font-family': 'helvetica',
                'font-size': '20px',
                'label': 'data(id)',
                'background-color': 'red'}
            },
            {
                'selector': 'node[id = "FANCD2"]',
                'style': {
                'font-family': 'helvetica',
                'font-size': '20px',
                'label': 'data(id)',
                'background-color': 'red'}
            },
            {
                'selector': 'node[id = "FANCE"]',
                'style': {
                'font-family': 'helvetica',
                'font-size': '20px',
                'label': 'data(id)',
                'background-color': 'red'}
            },
            {
                'selector': 'node[id = "FANCF"]',
                'style': {
                'font-family': 'helvetica',
                'font-size': '20px',
                'label': 'data(id)',
                'background-color': 'red'}
            },
            {
                'selector': 'node[id = "FANCI"]',
                'style': {
                'font-family': 'helvetica',
                'font-size': '20px',
                'label': 'data(id)',
                'background-color': 'red'}
            },
            {
                'selector': 'node[id = "PALB2"]',
                'style': {
                'font-family': 'helvetica',
                'font-size': '20px',
                'label': 'data(id)',
                'background-color': 'red'}
            },
            {
                'selector': 'node[id = "RAD51C"]',
                'style': {
                'font-family': 'helvetica',
                'font-size': '20px',
                'label': 'data(id)',
                'background-color': 'red'}
            },
            {
                'selector': 'node[id = "SLX4"]',
                'style': {
                'font-family': 'helvetica',
                'font-size': '20px',
                'label': 'data(id)',
                'background-color': 'red'}
            },
            {
                'selector': 'node[id = "UBE2T"]',
                'style': {
                'font-family': 'helvetica',
                'font-size': '20px',
                'label': 'data(id)',
                'background-color': 'red'}
            },
            {
                 'selector': 'node[degree>0]',
                 'style': {
                     'width': '100px',
                     'height': '100px'
                 }
            },
            {
                 'selector': 'node[degree>1]',
                 'style': {
                     'width': '150px',
                     'height': '150px'
                 }
            },
            {
                 'selector': 'node[degree>2]',
                 'style': {
                     'width': '200px',
                     'height': '200px'
                 }
            }
        ]
        )
        return cytoscapeobj
    