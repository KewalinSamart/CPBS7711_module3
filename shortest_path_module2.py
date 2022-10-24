import numpy as np
from network import *
from priority_queue import *
from utilities_module3 import *

#network = Network(filename, flip="no")
def shortest_path(network, start_node, query_nodes):
    # This function computes shortest path from the start_node to all the query_nodes
    # returns subnetwork interactions and cumulative weight(distance) to all query nodes from the start node
    shortest_paths = {}
    gene_list = list(set(list(network.network_df.gene1)+list(network.network_df.gene2)))
    distance_list = [np.inf]* len(gene_list)
    shortest_paths = dict(zip(gene_list, distance_list))
    shortest_paths[start_node] = 0
    # update visited_nodes
    network.update_visited_nodes(start_node)
    # get its neighbors
    node_distance_dict = network.get_neighbors(start_node)
    indiv_neighbors_dict = network.get_neighbors(start_node)
    # record individual paths for subnetwork output and visualization
    ##############################################
    # transform node_distance_dict to list = [(start_node, dest_node, weight),...]
    indiv_paths = [(start_node,k, v) for k, v in indiv_neighbors_dict.items()]
    indiv_shortest_paths = []
    ##############################################
    # record paths found # first itr -- no update-- just transform to {start_node: (dest_node, weight), ... }
    new_node_distance_dict = {}
    new_node_distance_dict[start_node] = [(k, v) for k, v in node_distance_dict.items()]

    # add it to priority queue
    prio_q = Priority_queue()
    prio_q.add_queue(start_node, new_node_distance_dict)
    priority = prio_q.get_priority() 
    prio_q.remove(priority)

    # update shortest_paths dict if the priority passes QC conditions
    # ****QC conditions:
    prio_dest = priority[1]
    prio_weight = priority[2]

    if (shortest_paths[prio_dest] == np.inf) or (shortest_paths[prio_dest] > prio_weight) :
        shortest_paths[prio_dest] = prio_weight
        network.update_visited_nodes(prio_dest)
        # set start_node = the destination node of the selected priority queue
        start_node = prio_dest
        ##########################
        # get index corresponding to the priority queue
        index = [tup[1] for tup in indiv_paths].index(prio_dest)
        indiv_shortest_paths.append(indiv_paths[index])
        #########################
    record_cum_paths = []
    while len(prio_q.queue) != 0:
        if (all(gene in network.visited_nodes for gene in query_nodes)):
            prio_q.clear_queue()
            network.clear_visited_nodes()
            break
        # start new iteration
        # look for its neighbors
        node_distance_dict = network.get_neighbors(start_node)
        indiv_neighbors_dict = network.get_neighbors(start_node)
        # record individual paths for subnetwork output and visualization
        ##############################################
        # transform node_distance_dict to list = [(start_node, dest_node, weight),...]
        # update individual paths list (include the previous one)
        indiv_paths = [(start_node,k, v) for k, v in indiv_neighbors_dict.items()] + indiv_paths
        ##############################################
        # update the weight by adding the shortest path weight of the current
        # use dict container to update: # make weights in node_distance_dict cumulative
        for node in node_distance_dict.keys():
            node_distance_dict[node] = node_distance_dict[node] + shortest_paths[start_node]
    
        # then transform it to {start_node: (dest_node, weight), ... }
        new_node_distance_dict = {}
        new_node_distance_dict[start_node] = [(k, v) for k, v in node_distance_dict.items()]
        
        # add it to priority queue
        prio_q.add_queue(start_node, new_node_distance_dict)
        priority = prio_q.get_priority() # (through, dest, cum_weight) -- same initial node throughout
        prio_q.remove(priority)
     
        prio_dest = priority[1]
        prio_weight = priority[2]
    
        if (shortest_paths[prio_dest] == np.inf) or (shortest_paths[prio_dest] > prio_weight) :
            shortest_paths[prio_dest] = prio_weight
            network.update_visited_nodes(prio_dest)
            record_cum_paths.append(priority)
            #prev_start_node = start_node
            start_node = prio_dest
            ##########################
            # get index corresponding to the priority queue
            index = [tup[1] for tup in indiv_paths].index(prio_dest)
            indiv_shortest_paths.append(indiv_paths[index])
            #########################
    
    return list(set(indiv_shortest_paths)), shortest_paths


# read in network separately and put it to the function
# network = Network(network_filename, flip = "yes") # yes to flip the weights coming with STRING network
def get_path_cost(network, start_node, dest_node):
    # This function computes path cost info and stores it in a dictionary for backtracking the path
    # get start node's neighbors
    # add neighbors to queue
    priority_q = Priority_queue()
    node_distance_dict = network.get_neighbors(start_node)
    new_node_distance_dict = {}
    new_node_distance_dict[start_node] = [(k, v) for k, v in node_distance_dict.items()]
    priority_q.add_queue(start_node, new_node_distance_dict)
    # get the priority queue
    priority = priority_q.get_priority()
    # priority = (prev, current, weight)
    priority_q.remove(priority)
    # initialize path-cost dictionary 
    gene_list = sorted(list(set(list(network.network_df.gene1)+list(network.network_df.gene2))))
    path_cost_list = [None]* len(gene_list)
    gene_path_cost_dict = dict(zip(gene_list, path_cost_list))
    gene_path_cost_dict[priority[1]] = (priority[0],priority[2])
    network.update_visited_nodes(priority[0])
     
    while len(priority_q.queue) != 0:
        next_start_node = priority[1]
       
        # get neighbors
        node_distance_dict = network.get_neighbors(next_start_node)
        new_node_distance_dict = {}
        new_node_distance_dict[next_start_node] = [(k, v) for k, v in node_distance_dict.items()]
        # add neighbors to the queue
        priority_q.add_queue(next_start_node, new_node_distance_dict)
        # get the priority queue
        priority = priority_q.get_priority()
        # remove the priority from queue
        priority_q.remove(priority)
        ### conditions to update weight
        # if previous weight in the dict exists (continue from an existing path) -> update the weight 
        if (gene_path_cost_dict[priority[0]] != None):
            weight = priority[2] + gene_path_cost_dict[priority[0]][1]
            chosen_index = 1
            if gene_path_cost_dict[priority[1]] != None: 
                #check which way is shorter
                weight1 = priority[2] + gene_path_cost_dict[priority[0]][1]
                weight2 = priority[2] + gene_path_cost_dict[priority[1]][1]
                # if from priority[0] is shorter then update the gene_path_cost_dict(priority[1])
                if weight1 < weight2:
                    weight = weight1
                    chosen_index = 1
                # if from priority[1] is shorter then update the gene_path_cost_dict(priority[0])
                if weight1 > weight2:
                    weight = weight2
                    chosen_index = 0  
            
        else: # if not, just add the current weight
            weight = priority[2]
            chosen_index = 1     # priority[0]
        ### QC conditions to update a shortest path
        if chosen_index == 1:
            if (gene_path_cost_dict[priority[chosen_index]] == None) or (gene_path_cost_dict[priority[chosen_index]][1] > weight):
                gene_path_cost_dict[priority[chosen_index]] = (priority[0], weight)
                network.update_visited_nodes(priority[chosen_index]) # priority[0]
            else:
                pass
        elif chosen_index == 0:
            if (gene_path_cost_dict[priority[chosen_index]] == None) or (gene_path_cost_dict[priority[chosen_index]][1] > weight):
                gene_path_cost_dict[priority[chosen_index]] = (priority[1], weight)
                network.update_visited_nodes(priority[chosen_index]) # priority[1]
            else:
                pass
        if next_start_node == dest_node:
            network.clear_visited_nodes()
            priority_q.clear_queue()
            break
    gene_path_cost_dict[start_node] = (start_node, 0)

    return gene_path_cost_dict

def get_path(gene_path_cost_dict, start_node, dest_node):
    # This function takes the out put from get_path_cost func  
    # returns nodes in order visited on the shortest path from start node to destination node
    # remove dict elements with value = None
    gene_path_cost_dict = {key: value for key, value in gene_path_cost_dict.items() if value is not None}
    # initialize path_list to store path
    path_list = [dest_node]
    # backtrack to get previous nodes starting from the destination node to the start node
    while dest_node != start_node:
        prev_node = gene_path_cost_dict[dest_node][0] 
        path_list.insert(0,prev_node)
        dest_node = prev_node     
    return path_list    

def get_interactions_subnetwork(network_df, path_list):
    # This function all get gene-gene interactions in a given network data frame of the genes 
    # in the given path list (output of get_path func)
    interactions_subnetwork = []
    for pos in range(len(path_list)-1):
        gene_pair = tuple(sorted(list((path_list[pos],path_list[pos+1]))))
        weight = network_df.get_network_interactions_degrees()[0][gene_pair]
        interactions_subnetwork.append(gene_pair+(weight,))
    return interactions_subnetwork

def get_final_subnetwork(network, query_nodes):
     # This function generate a subnetwork based on the given query nodes
    subnet_interactions = []
    for start_node in query_nodes:
        indiv_shortest_path = shortest_path(network=network,start_node = start_node, query_nodes = query_nodes)[0]
        primary_subnetwork = Network(paths_to_df(indiv_shortest_path))
        for dest_node in query_nodes:
            gene_path_cost_dict = get_path_cost(primary_subnetwork,start_node, dest_node )
            try:
                paths = get_path(gene_path_cost_dict, start_node=start_node, dest_node=dest_node )
                subnet_interactions.extend(get_interactions_subnetwork(primary_subnetwork,paths))
            except:
                pass
            primary_subnetwork.clear_visited_nodes()
    return list(set(subnet_interactions))