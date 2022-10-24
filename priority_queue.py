import pandas as pd
import numpy as np

class Priority_queue:
    def __init__(self):
        self.queue = {}
        
    def add_queue(self, start_node, node_distance_dict):
        # add element to the priority queue
        neighbor_list = node_distance_dict[start_node]
        for neighbor in neighbor_list:
            self.queue[(start_node, neighbor[0])] = neighbor[1]
     
    def get_priority(self):
        # get key of the min(weight) 
        key_min = min(self.queue, key=self.queue.get)
        val_min = self.queue[key_min]
        # set to priority queue to (gene1, gene2, min(weight))
        priority_key = key_min + (val_min,)
        
        return priority_key
    
    def remove(self, priority):
        # remove the element that has the most priority
        del self.queue[priority[:2]]

    def clear_queue(self):
        # delete all the queue
        self.queue.clear()    
