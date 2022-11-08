from scoring_genes import *
from visualize_finalsol import *
from network import *
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Given a population of solutions for a given set of loci, this script is going to score genes 
    on the loci using the method in Tasan et al. 2015 and visualize the final solution.
    """)
    parser.add_argument("-solutions", help="(n+1)-clumns (n+1 loci; where n = 0,1,..) containing solutions (a chosen gene per locus; set to example_solutions.txt by default)")
    parser.add_argument("-num_loci",help="number of loci in the solutions; set to 12 by default")
    parser.add_argument("-network", help="txt file containing gene-gene interaction network (undirected; can be weighted/unweighted, but weights will not be used in gene scoring); set to STRING_network.txt by default")
    parser.add_argument("-output_dir", help="path to store final output; set to example_result/final_solution.txt")
    parser.add_argument("-score_cutoff", help="score cutoff for circular-layout (via Networkx) solution visualization; set to 0.25 by default")


    args = parser.parse_args()
    solutions_filename = args.solutions
    num_loci = args.num_loci
    network_filename = args.network
    output_dir = args.output_dir
    score_cutoff = args.score_cutoff

    final_sol = get_final_solution(solutions_filename, num_loci, network_filename, output_dir)    
    finalsol_viz(final_sol, Network(network_filename).network_df, num_loci, score_cutoff)
   
    

   


