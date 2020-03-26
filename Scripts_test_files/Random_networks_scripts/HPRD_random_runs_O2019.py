# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 18:34:21 2019

@author: Admin
File name: HPRD_random_runs.py
Replaces all previous versions of Masterscripts
Node removal with random genelist based filtering

Inputs: 1. Disease-gene list
        2. PPI file
        3. Bet_cent_parallel to run the betweenness centrality algorithm
WARNING: This script needs to be run in an external window. Spyder has issues with multiprocessing module    

Notes: For AppendRandlisttable, main1() returns the randlist, which is collected by the variable f1
and f1 is then supplied as an argument to main3(f1) for the table to get appendend with the current
genelist

\
"""

import sys
import os
import glob
from collections import OrderedDict
import networkx as nx
import time
import csv
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import Bet_cent_parallel_2
import logging
import importlib
import json

logging.basicConfig(level = logging.DEBUG, filename ='logfile.log',format='%(asctime)s - %(message)s')



Option = 3
jobID = '04' #for local desktop runs
#jobID =sys.argv[1]   # for HPC runs
print ("Saving to file " + jobID)
dbname = "HPRD" #Option: HPRD, BioGRID,MINT
basefilename = jobID + "CV_O2019_%s" %dbname      #Used to name all the files generated
betweenName = jobID + "betweeness"      #Naming convention, required while running multiple jobs               
degreeName = jobID + "degree"
Setno = 'S1'


"""
#Reads a file into a list, closes the original file. Lists used for further computation. Input: File to read
#Removes empty lines with the linein.split() conditions
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

#Writing files: without and with space
def writefilewospace(filenametowrite, datatowrite):
    with open("{0:s}.txt".format(filenametowrite), "w+") as f: f.writelines("".join(datatowrite))

def writefile(filenametowrite, datatowrite):
    with open("{0:s}.txt".format(filenametowrite), "w+") as f: f.writelines("\n".join(datatowrite))

def writefilewithspace(filenametowrite, datatowrite):
#    datatowritewithspace = [item +'\n' for item in datatowrite]
    with open("{0:s}.txt".format(filenametowrite), "w+") as f: f.writelines("\n".join(str(item) for item in datatowrite))
    

"""
#PPIextraction: Extracts disease specfic network based on the input diseasegene list from the PPI file
Inputs: diseasegenelist - Initial seed list. Input to the function as list after reading the file
        ppifile - PPI file read
        filetosave - basefilename
        returns NR PPI list for Node extraction for adding Connected components
"""
def PPIextraction(diseasegenelist,filetosave,l=None):
    if l is None:
        genelist2 = []  #Clears the variables for each run
        genelist = [] 
        ppilist = [] 
        result = []
        result3 = []
        result2 = []
    genelist2 = diseasegenelist  #Genelist
    genelist = [item15 for item15 in genelist2 if item15 != '\n'] #Removes empty line if present
    with open('PPI_dict_HPRD_GC_NR.txt', 'rb') as handle:
        b = pickle.loads(handle.read())

    print('\tPPI Extraction')
    start = time.time()
    for item3 in genelist: 
        if item3.split() != []: 
            item4 = item3.split()[0]    #Genename from genelist
        for k, vals in b.items():
                if item4 == k:          #k: key of the dictionary
                    ppilist.append(vals[:]) #if gene matches key, collect all interactions
    result = sum(ppilist, [])      #merge separate lists into a single PPI list           
    result3 = list(set(result))
    print("\t\tTime: %.4F" % (time.time() - start))  

    filenamebydb = "%s_PPI" %filetosave #"17082018_%s" %dbname
    writefile(filenamebydb,result)
    return result3

"""
AddConnectedComponents extracts all nodes from the above extracted PPI, removes the MD nodes. From
among the first neighbours letf, collects interactions between any two first neighbours
"""    
def AddConnectedComponents(res3,diseasegenelist,ppifile,filetosave, lr = None):
    if lr is None:
        MDgenelist = []
        PPIfullfile = []
        MDlist = []
        Pfile = []
        GL,G1,G2,G3,G4,G5,FNeighbors,PPIfilt,PPIsel,TotPPI_red, TotPPIset = ([] for i in range(11))
       
    MDgenelist = diseasegenelist  
    PPIfullfile = ppifile
    MDlist = [item01.split()[0] for item01 in MDgenelist]  # MD nodes
    Pfile = res3

    GL = [item.split() for item in Pfile]   # Split two interactions
    G1 = [item2[0] for item2 in GL]         # Collect first interactor
    G2 = [item3[1] for item3 in GL]         # Collect second interactor
    G5 =list(set(G1+G2))

    FNeighbors = [item for item in G5 if item not in MDlist]      # collect all non-MD nodes from above list

    PPIfilt = [item03 for item03 in PPIfullfile for item04 in FNeighbors if item04 in item03] #Rough filtering
    # PPI sel: Add interactions only if both interactors are first neighbours
    PPIsel = [item02 for item02 in PPIfullfile if item02.split()[0] in FNeighbors and item02.split()[1] in FNeighbors]
    TotPPI =list(set(res3 + PPIsel)) # Combine the previous PPI and added components 
    Logstr = "Number of nodes:{0},Number of Interactions:{1}".format(str(len(G5)),str(len(TotPPI)))
    logging.info(f'{Logstr}')
    filenamebydb = "%s_TotPPI" %filetosave      # This file then used to create giant component
    writefilewospace(filenamebydb,TotPPI)
    writefile('FNnodes',FNeighbors)
    writefilewospace('CC',PPIsel)

def createGiantComponent(filenametoread,filenametowrite):
    try:
        import matplotlib.pyplot as plt
    except:
        raise
    try:
        from networkx import graphviz_layout
        layout=nx.graphviz_layout
    except ImportError:
        print("PyGraphviz not found; drawing with spring layout; will be slow.")
        layout=nx.spring_layout
    G = nx.read_edgelist("%s_TotPPI.txt" %filenametoread)   # Using the complete PPI from above 
    Gcc=sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)
    G0=Gcc[0]
    nx.write_edgelist(G0, "%s_GiantComponent.txt" %filenametowrite)



"""
AppendTable collects centrality values
"""

def AppendTable(filename, lj = None):
    if lj is None:
        gene = []
        xyz = []
        xwz = []
    print('\tAppendtable')
    start = time.time()
    ReadBetfile = readfile('%s.txt'%filename)
    print (ReadBetfile[0:5])
    gene = [item1.split(':') for item1 in ReadBetfile]
    xyz = [item.append(gene2[1].split()[0]) for item in LofGenes for gene2 in gene if item[0] == gene2[0]]
    xwz = [item4.append('Ab') for item4 in LofGenes if len(item4) !=(i+2)]
    print("\t\tTime: %.4F" % (time.time() - start))
    return LofGenes

def AppendRandListTable(ijk, lj = None):
    if lj is None:
        gene = []
        abc = []
        ghi = []
    print('\tAppendRandlist')
    start = time.time()
    ReadBetfile = ijk  
    print (ReadBetfile[0:5])
    genetrim = [item14 for item14 in ReadBetfile if item14 != '\n']
    gene = [item1.split()[0] for item1 in genetrim if item1.split() != []]
    abc = [item5.append(gene3) for item5 in LofRandList for gene3 in gene if item5[0] == gene3]
    ghi = [item6.append('Ab') for item6 in LofRandList if len(item6) !=(i+2)]
    print("\t\tTime: %.4F" % (time.time() - start))
    return LofRandList

def writecentvalues(datafromappendtable,loopno):
    print('\tWriteCentVals')
    start = time.time()
    testdf = pd.DataFrame(data = datafromappendtable) 
    fname =  "{0:s}_BetweenessMeasures_{1:s}.csv".format(str(Option),jobID)
    testdf.to_csv(fname)
    print("\t\tTime: %.4F" % (time.time() - start))

def writerandlistable(datafromappendrandlisttable,loopno):
    print('\tWriteRandlistTable')
    start = time.time()
    testdf = pd.DataFrame(data = datafromappendrandlisttable) 
    fname =  "{0:s}_RandListTable_{1:s}.csv".format(str(Option),jobID)
    testdf.to_csv(fname)
    print("\t\tTime: %.4F" % (time.time() - start))

def calculateDegree(filenametoread,filenametowrite):
    G = nx.read_edgelist("%s_GiantComponent.txt" %filenametoread)
    G.remove_edges_from(G.selfloop_edges())
    nx.write_edgelist(G, "%s_processed_GC.txt"%filenametowrite)
    d = dict(G.degree)
    n = OrderedDict(sorted(d.items(), key=lambda t: t[1], reverse=True))
    with open ('%s_degree.txt'%filenametowrite, 'w') as fp:
        for p in n.items():
            fp.write("%s:%s\n" % p)

"""
processDegreefile_c and NewGiantComponent_c applicable for case of removing single degree
nodes except the ones in rand gene lists
"""

def processDegreefile_c(diseasegenelist, filenametoread,filenametowrite,l=None):
    if l is None:
        linestoread = []
        rejectedlist = []
        Rejlist = []
        Rejgenelist =[]
        genelistwSpace = []
        genelist = []
        rejectedlist_1 = []
    genelst = list(set(diseasegenelist))
    genelist = [item10 for item10 in genelst if item10 != '\n']
    print (genelist[0:5])
    filename = '%s_degree.txt'%filenametoread
    linestoread = readfile(filename)
    rejectedlist_1 = [item3.split(':')[0] for item3 in linestoread if item3.split(':')[1].split()[0]=='1' or item3.split(':')[1].split()[0]=='0']
    print(rejectedlist_1[0:5])
#    rejlist = [item12 for item12 in genelist if item12 =! '']
    rejectedlist = [item9 for item9 in rejectedlist_1 for item10 in genelist if item10.split() != [] if item9 == item10.split()[0] ]
    print("Secondlist",rejectedlist[0:5])
    trim = [rejectedlist_1.remove(item11) for item11 in rejectedlist]
    Rejgenelist = list(set(rejectedlist_1))
    Removegenes = "GeneswDeglessthan1_%s" %filenametowrite
    writefile(Removegenes,Rejgenelist)


def NewGiantComponent_c(filenametoread,filenametowrite, p2=None):
    if p2 is None:
        removedGC1 = []
        removedGC2 = []
        readGCfile = []
        dxy = []
        removedGC = []
        readGCfile = []
        filtGCfile = []
        readRejgenefile = []
        remGCset = []
    print('\tNew GC')
    start = time.time()  
    readGCfile = readfile("%s_processed_GC.txt" %filenametoread)
    readRejgenefile = readfile("GeneswDeglessthan1_%s.txt" %filenametoread) 
    filtGCfile = [iteminGC for iteminGC in readGCfile for remgene in readRejgenefile if remgene.split()[0] in iteminGC]
    removedGC1 = [item4 for item4 in filtGCfile for item in readRejgenefile if item.split()[0] == item4.split()[0]]
    removedGC2 = [item6 for item6 in filtGCfile for item7 in readRejgenefile if item7.split()[0] == item6.split()[1]]
    
    remGClst = list(set(removedGC1 + removedGC2))
    dxy = [readGCfile.remove(item) for item in remGClst]
    print("\t\tTime: %.4F" % (time.time() - start))
    print('\tNew GC-write')
    start = time.time()
    writefilewospace("RevisedGC2_%s_GiantComponent"%filenametowrite,readGCfile)
    print("\t\tTime: %.4F" % (time.time() - start))

"""
processDegreefile_c and NewGiantComponent_c applicable for case of removing single degree
nodes except the ones in rand gene lists
"""

def main1c():
    from Binning_optimization_v3 import RandomGeneList
    print("Step - 1")
    Diseasegenelist = RandomGeneList(c,i)
    print("Step - 2")
    PPIext = PPIextraction(Diseasegenelist,basefilename)  #Geneslist_NR2,Geneslist_NR2_16082018,1920_unique_genes 
    CC = AddConnectedComponents(PPIext,Diseasegenelist,PPIfile,basefilename)
    print("Step - 3")
    FirstGC = createGiantComponent(basefilename,basefilename)
    Firstdeg = calculateDegree(basefilename,basefilename)
    FilteredDeg = processDegreefile_c(Diseasegenelist,basefilename,basefilename)
    SecondGC = NewGiantComponent_c(basefilename,basefilename)
    print("Step - 4")
    return Diseasegenelist

def main2b():
    import Bet_cent_parallel_2
    print("Step - 5")
    reffilename = "RevisedGC2_%s" %basefilename
    print("Step - 6")
    if __name__ == "__main__":
        print("Step - 7")
        Bet_cent_parallel_2.runcent(reffilename,jobID)

def main3(ijk):
    tableinput = AppendTable(betweenName)
    randlistinput = AppendRandListTable(ijk)
    time.sleep(2)


#def runloop():
for c in list(range(1)):
    Nodeslist = readfile('All_Nodes_HPRD.txt')
    SFulllist = sorted(list(set(Nodeslist)))
    timestamp = time.time()
    LofGenes = [item.split() for item in SFulllist]
    LofRandList = [item2.split() for item2 in SFulllist]
    PPIfile = readfile('HPRD_GC_PPI.txt')

#    PPIfile = readfile(ppifilename)
    for i in list(range(1)):
        if __name__ == "__main__":
            f1= main1c()
            main2b()
            main3(f1)
            if i == 1:
                writecentvalues(LofGenes,timestamp)
                writerandlistable(LofRandList,timestamp)                          

