Workflow

Extraction of genes from CTD raw data

PPI conversion:
For HPRD, run HGNC conversion directly
For BioGRID, run Select_Human_Int, then run HGNC conversion using the output file.

Select_Human_Int : Dependencies: BioGRID raw datafile
Output: BioGRID interactions file 2 column format

HGNC conversion

Dependencies: Gene list (1 column gene symbols)
PPI network files:
HPRD - raw datafile
BioGRID - Select_Human_Int output file (2 column format)


