# MD_network_map
**Project Description**

These scripts accompany and form the basis of the study of metabolic diseases (MD) network 
titled 'Degree adjusted large scale network analysis reveals novel putative metabolic disease genes'
The study uses disease-gene data from the Comparative Toxicogenomics Database (CTD)
and two protein-protein interaction networks (PPINs) - HPRD and BioGRID - to predict novel 
MD genes based on network centrality. 


**Platform/Requirements**

* Python v3.7
* Networkx library v2.3


**Script usage (In order of requirements)**

*Data Pre-processing*
1. Extract_genes_fromCtd_v4.py
2. Sel_Human_Int.py (only for BioGRID PPIN)
3. HGNC_conversion_O2019_v3.py
4. PPI_to_dict.py

*Network generation and centrality calculations*

5. MD_analysis.py (calls Bet_cent_parallel_2.py)
6. Calls.py (calls Binning_optimization_v3.py)
7. HPRD_random_runs_O2019.py (calls Binning_optimization_v3.py)

*Analysis* 

8. Modify_csv.py
9. Combine_csv_files_v2.py
10. Analysis.py
11. Adj_P_values.py


**Data Requirements**

 * CTD disease-genes raw data 
 * HPRD raw data 
 * BioGRID MV raw data

**Notes**

To follow this method, we need betweenness centrality analysis for two kinds 
of networks: the Metabolic Diseases network (MD) and the Randomly constructed 
networks (Random networks)

*Data Pre-processing*
For the MD network, we downloaded the disease-gene lists from CTD.
For the 4 categories considered under MD, we extract the gene symbols, 
(Extract_genes_fromCtd_v4.py). (After obtaining the gene list for each category,
they need to be combined manually into a single file. This file can be 
processed using Extract_genes_fromCtd_v4.py to ensure non-redundant MD gene list)

We use HPRD and BioGRID MV datasets for the protein-protein interaction network
(PPIN) data.

For BioGRID, we use Sel_Human_Int.py to filter out human interactions. The output
from this script is used for converting to HGNC based symbols (HGNC_conversion_O2019_v3.py)

For HPRD, HGNC_conversion_O2019_v3.py has a function to remove single interactions.

HGNC_conversion_O2019_v3.py takes in as input HGNC file with approved symbols, old and alternate 
symbols. It first creates a reference file, in a two column format, with approved symbols in one column
and old/alternate symbol in the second.
This reference file is then used to convert the MD gene list as well as PPINs to standard HGNC symbols.

After HGNC conversion, the script PPI_to_dict.py is used to convert PPI (which is a two column format 
with a list of interactions) to a dictionary format to speeden up network construction. 

*Network generation and centrality calculations*

For the networks, the MD analysis script can be used to obtain betweenness centralities for MD network
genes. It calls the 'Bet_cent_parallel_2.py' The output file '03betweenness.txt' will be required for the analysis.

Following the MD network centralities, one would need a set of random network centralities to compare them
with and calculate p-values.
The scripts in the Random Network Analysis should be used in the following order:
1. First use Calls.py to generate stratification files. This script calls Binning_optimization_v3.py
2. Then use the HPRD_random_runs_O2019.py to construct and analyse a user defined number of random networks.
This script also calls Binning_optimization_v3.py for generating a random gene list for each run.
The output of this script are two .csv files, one containing all the genelists for each run and the other 
containing all the betweenness centralities. 

*Analysis*

Random networks for this manuscript were generated in chunks, and hence, to combine 
the output of different runs, two scripts were used:
Modify_csv.py - Removes single lines from random network output .csv files and renames them
Combine_csv_files_v2.py - Combines several random run files into a single file.
Need to run this file twice, once for genelist files and once for files with betweenness 
centralities. 
 
The output of this processing is used for analysis.
 This script, which requires the betweenness centrality values for the MD network:
1. Selects centralities based on presence in random list. Outputs a file 
called 'SelectedData_'
2. Finally, it outputs a file containing gene name, number of hits, total number
of occurences, p-values, Average centrality in random networks, Centrality in MD
and Relative centrality (MD centrality/Average)

This file with p-values can be then used as input for calculating the adj-p values
This script removes all the entries where the average betweenness centrality is
0 and calculates adjusted-p values.  

**Licence**

This code is freely available for non-commercial users under the GPLv3 licence


**Maintainence**

This repository is maintained by the Systems Biology Group, Department of Life Sciences and Medicine, University of Luxembourg
For questions or comments, please contact:

apurva.badkas@uni.lu
