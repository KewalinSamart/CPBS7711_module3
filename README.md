# CPBS7711_module3
From a population of solutions for a given set of loci, score the genes on the loci using the method in Tasan et al. and visualize the final solution

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

2. A population of Prix Fixe (PF) solutions; `.txt` (`\t` separated) file containing n+1 columns where n is in {0,1,...} with their corresponding chosen gene (one gene per locus for a single solution); number of rows = number of PF solutions 

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

- Example output 
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


### Command
```{r}
python3  
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
