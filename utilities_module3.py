# by default, this program would generat random prix fixe solutions using the gmt file for FA 
def get_loci_candidate_genes(loci_genes_filename = "Module1_FA_Day3_input.gmt.txt"):
    '''
    This function read in txt file with loci and their candidate genes
    input file format:  
    first colunm: loci index e.g. 0, 1, 11
    second column and so on: gene candidate at the specific locus index e.g. AC092291.2
    Returns a dict with keys: loci indices, values: lists of gene candidates
    '''
    myfile = open(loci_genes_filename)
    loci_candidate_dict = {}
    loci_index = 0
    for line in myfile :
        if loci_genes_filename == "Module1_FA_Day3_input.gmt.txt":
            candidate_genes = line.split("\t")[2:]
            candidate_genes = [elem.replace('\n', '') for elem in candidate_genes]
        else:
            candidate_genes = line.split("\t")[1:]
        loci_candidate_dict[loci_index] =  candidate_genes
        loci_index = loci_index + 1
    return loci_candidate_dict

