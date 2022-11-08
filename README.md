# CPBS7711_module3
From a population of Prix Fixe (PF) solutions for a given set of loci, score the genes on the loci using the method in Tasan et al. and visualize the final solution

## Workflow
### Inputs
1. Undirected un/weighted network file in `.txt` format (`\t` separated) containing gene-gene interactions.
  - `column 1&2`: gene names
  - `column 3`: interaction strengths range from 0 to 1.
  - Example input `STRING_network.txt` –– weighted network:
```{r}
BLOC1S6 BLOC1S3	  0.24
RAB3D   CHML      0.847
MYL7    MYO15A    0.842
...
```
**NOTE:** For unweighted network, network file would be exact the same format as above but without the weight column.

2. A population of PF solutions; `.txt` (`\t` separated) file containing n+1 columns where each column is named 0,1,...,n indicating locus index and contains their corresponding chosen genes (one gene per locus for a single solution (i.e. row)); number of rows = number of PF solutions 
  - Example input `toy_loci_set.txt`
```{r}
0 1 2 3 4 5 6 7 8 9 10 11
XPO6 FBXO3 C17orf64 CORO2A DEF8 CNTN2 RPL32 FANCE FREM2 PARN ANPEP CDIP1
KDM8 CAT DDX5 NR4A3 GAS8 PPFIA4 BRK1 TULP1 SERTM1 RPS15A MESP1 RP11-127I20.4
PLK1 SVIP GDPD1 NANS AC133919.6 SOX13 FANCD2 SLC26A8 SERTM1 NPIPP1 SEMA4B HMOX2
...
```

### Output
1. A `.txt` (`\t` separated) file containing genes, their final scores, and the locus that each gene belongs to
- `column 1`: gene names
- `column 2`: final gene score; higher scores indicate better contributions 
- `column 3`: associated locus (locus index)

- Example output `final_solution.txt`
```{r}
gene     score                 locus
PALB2    0.13280434328669868   0
NUPR1    0.02125051082958725   0
SLC5A11  0.016724454415663878  0
CCDC73   0.0                   1
CAPRIN1  0.02042483660130719   1
RCN1     0.05896805896805897   1
...
```
2. Solution visualization in `png` format

Visualization details:
  - Each circle represents a gene
  - Different colors indicate different loci that the genes belong to
  - Circle's size represents their final gene score
  - Edges are gene-gene interaction based on the input network

    2.1 full final solution with no score cutoff applied –– kamada kawai layout
      - Example output `example_result/example_finalsol_kkviz.png`:
      ![alt text](https://github.com/KewalinSamart/CPBS7711_module3/blob/main/example_result/example_finalsol_kkviz.png?raw=true)
    2.2 highest-scores solution with a score cutoff (0.25 by default) –– circular layout
      - Example output `example_result/example_finalsol_ccviz.png`:
      ![alt text](https://github.com/KewalinSamart/CPBS7711_module3/blob/main/example_result/example_finalsol_ccviz.png?raw=true)

3. Final solution subnetwork in `.json` for user's customization in [`Cytoscape`](https://cytoscape.org/) 
  - Example output `example_finalsol_json.js`:
 ```{r}
 {"data": [], "directed": false, "multigraph": false, "elements": {"nodes": [{"data": {"locus": 8, "score": 0.0358671285918076, "color": "#64678B", "id": "AKAP11", "value": "AKAP11", "name": "AKAP11"}}, {"data": {"locus": 11, "score": 0.1770076779414816, "color": "#D5A612", "id": "SLX4", "value": "SLX4", "name": "SLX4"}}
 ```
 
### Command
```{r}
python3 "final_scored_solution.py" -solutions [solutions file name (string)] -num_loci [loci number (int)] -network [network file name (string)] -output_dir [path to output directory (string)] -score_cutoff [score cutoff (float)]    
```

### Arguments
- sys.argv[0]: 
- sys.argv[1]: 
- sys.argv[2]: 
- sys.argv[3]: 
- sys.argv[4]: 

### Example command
```{r}
python3 
```

## Installation and Dependencies
- `Python 3.8.3`
- `pandas 1.5.0`
- `numpy 1.23.3`
