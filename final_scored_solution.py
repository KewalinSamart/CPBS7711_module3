from scoring_genes import *
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Given a population of solutions for a given set of loci, this script is going to score genes 
    on the loci using the method in Tasan et al. 2015 and visualize the final solution.
    """)
    parser.add_argument("-solutions", help="(n+1)-clumns (n+1 loci; where n = 0,1,..) containing solutions (a chosen gene per locus; set to example_solutions.txt by default)")
    parser.add_argument("-network", help="txt file containing gene-gene interaction network (undirected; can be weighted/unweighted); set to STRING_network.txt by default")
    parser.add_argument("-flip_weight", help="a string either 'yes' or 'no' indicating whether or not to flip weights in the network; set to 'no' by default")
    parser.add_argument("-output_dir", help="path to store final output; set to example_result/final_solution.txt")

    args = parser.parse_args()
    solutions_filename = args.solutions
    network_filename = args.network
    flip_weight = args.flip_weight
    output_dir = args.output_dir

    get_final_solution(solutions_filename, network_filename, flip_weight, output_dir=output_dir)

   


