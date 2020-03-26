# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 22:17:31 2018

@author: Admin
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 14:32:30 2018
@author: Admin

On searching under 'Batch query' in the CTD database, for specific MeSH Disease IDs, 4 files are obtained.
This code takes these files for processing further.
Part 1: From the entire file, only obtain the gene symbols
Part 2: Combine all the four different lists into a single file (To be done manually)
Part 3: Ensure non-redundancy

Update: File corrected on 27.10.2018. Check comments below 

"""

def readfile(filetoread, p=None):    
    if p is None:
        listoflines = []
    openfile = open(filetoread)
    outputpar = openfile.readlines()
    listoflines = [linein for linein in outputpar if linein.split() != [] and linein.split() != '' ]    
    openfile.close()
    listoflines[-1] = listoflines[-1] + '\n'    
    return listoflines

def writefilewospace(filenametowrite, datatowrite):
    with open("%s.txt" %filenametowrite, "w+") as f: f.writelines("".join(datatowrite))


#categories = ['Liver','Metabolic','Malnutrition','Overnutrition'] #To run through all 4 files
categories = ['Liver_Diseases','Malnutrition'] #To run through all 4 files

filelist = []
glist = []
# From query files, obtain only the gene symbols
for disease in categories:
    Ctd_file = open("%s_Test.txt" %disease,'r') #Metabolicdiseases_CTDquery, Liverdiseases_CTDquery,Overnutrition_CTDquery
    Ctd_list = Ctd_file.readlines()
    del Ctd_list[0]                 # Delete header. Each query file has a header
# Change in code: 
    genelistR = []                  # Empty list initialization 
    for item in Ctd_list:
        item.strip()                # Remove trailing spaces if present   
        linebits = item.split('\t') # Split tab seperated line
        genesymbol = linebits[3]    # Gene symbols occur in the 4th column
        g3 = genesymbol.split()[0]     #Remove trailing spaces if any
        genelistR.append(g3)        # append to genelist
    l3 = set(genelistR)             # Remove duplicates if any
    l4 = list(l3)                   # l4 is a non-redundant list
    glist.append(l4)                # glist is collecting all lists from 4 files
    l42 = [item +'\n' for item in l4] # Add newline character for writing out the list
    filetoname = "%s_list_extracted.txt" %disease
    with open(filetoname, "w+") as f: f.writelines(l42)
    Ctd_file.close()
    genelistR = [] # Change: Reset list after writing to file, so genes from one list do not get carried to next
    filelist.append(filetoname)    # Collecting all filenames to concatenate data
    time.sleep(2)

# Ensure list is non-redundant

# Combine the independent gene lists into a single file. Use this single file as the input for below.
MDentfileread = readfile('Combined_Red_Test.txt')
MDentfile = [item23.split()[0] for item23 in MDentfileread]
setMDentfile = set(MDentfile)
lstMDentfile = list(setMDentfile)
l44 = [item +'\n' for item in lstMDentfile]
writefilewospace('Combined_Red_Test_NR',l44)
