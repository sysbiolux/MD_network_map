# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 16:12:32 2019

@author: apurva.badkas
"""

# -*- coding: utf-8 -*-
"""
This script does the following:
    Binning: Take a file with nodes and degrees and output files into different stratified lists based on the user defined min separation
    value
    It also outputs a file containing information on the number of genes in each of the stratified lists
    And a dictionary file containing the information of the degrees included in each of the stratified lists
    
    MDStratification: Reads in the dictionary file with degree information, also MDgenelist
                        Outputs a degree distribution list for MD genes
"""

# Min no of genes in a bin = 3
# Min separation between degrees = 5

import time
import csv
import pickle
import json


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

#Degree
def calculateDegree(filenametoread,filenametowrite): 
    import networkx as nx
    from collections import OrderedDict
    G = nx.read_edgelist("%s.txt" %filenametoread)
    G.remove_edges_from(G.selfloop_edges())
    nx.write_edgelist(G, "%s_processed_GC.txt"%filenametowrite)
    d = dict(G.degree)
    n = OrderedDict(sorted(d.items(), key=lambda t: t[1], reverse=True))
    with open ('%s_degree.txt'%filenametowrite, 'w') as fp:
        for p in n.items():
            fp.write("%s:%s\n" % p)

#Binning is for the entire PPI
#Binning requires a degree file input
def Binning(Degfile, min_sep):
    from Binning_optimization_v3 import readfile
    DegFile2 = readfile(Degfile)
    RangeDegs = [int(item01.split(':')[1]) for item01 in DegFile2 if item01.split()!=''] # Collect all degrees
    MaxDeg = max(RangeDegs)
    LstRange = list(set(RangeDegs)) # Non-redundant range of degrees
    LofRange = [[item04] for item04 in LstRange if item04 != 0]  #Used to build the dictionary below
    #DictofDegs = [item02.append(NofGene.split(':')[0]) for item02 in LofRange for NofGene in DegFile if int(NofGene.split(':')[1])==item02[0]]
    
    #collect all genes with the same degree
    for item08 in LofRange:
        for item09 in DegFile2:
            if item08[0] == int(item09.split(':')[1]):
                item08.append(item09.split(':')[0])
    bins = list()  #Contains all the stratified sets
    current_bin = list() #Container for single stratified layer
    min_len = min_sep  #min no of degrees per bin
    count = 00
    strcounter = 1
    totalcount = 00
#    DegsInStrat = []
    DegListToWrite = []
    DegsInStrat = {}
    BigDict = {}
    ListItemsForEachCategory = []
    sorLofRange = sorted(LofRange,reverse=True) # Starting from the highest degree

    for items in sorLofRange:
        if count < min_len:
            lofitems = len(items)-1  #excluding degree (int)
            current_bin.append(items) 
            c = list(sum(current_bin[:], []))
            count = count + 1
            totalcount = totalcount + lofitems
            if (count >= min_len and len(c)>= 2*min_len) or (sorLofRange.index(items) == len(sorLofRange)-1 ):  #min 3 items in a bin
                GeneNames = [item10+'\n' for item10 in c if isinstance(item10, str)]
                GeneNames2 = [item12 for item12 in c if isinstance(item12, str)]
                bins.append(GeneNames)
                ListNo = 'List_{}'.format(strcounter)
                DegsInStrat[ListNo] = [item11 for item11 in c if isinstance(item11, int)]
                BigDict[ListNo] = GeneNames2
                ListItemsForEachCategory.extend([ListNo +':' +str(totalcount)])
#                DegListToWrite.extend([ListNo +':' +str(DegsInStrat[strcounter-1])])
#                writefilewospace(ListNo, GeneNames)
                current_bin = list()
                count = 00 
                totalcount = 00
                strcounter = strcounter + 1
    writefilewithspace('LengthOfLists', ListItemsForEachCategory)
    json.dump(DegsInStrat, open('DictofListsandDegs.txt','w'))
    json.dump(BigDict, open('AllGenesDict.txt','w'))
#    writefilewithspace('DegsDistPPI', DegListToWrite)
    return DegsInStrat #bins

def MDStratification(DegFileimp,MDgenes):
    #ReadMDGenelist
#    from Calls import DegFileimport, MDgenes, Stratfile
    from collections import Counter
    DegFile3 = readfile(DegFileimp)
    MDgenelist = readfile(MDgenes)
    Stratification = json.load(open('DictofListsandDegs.txt'))
    MDgeneDegrees = [item20 for item20 in DegFile3 for item19 in MDgenelist if item19.split()[0] ==item20.split(':')[0]]
    MDdegs = [int(item22.split(':')[1]) for item22 in MDgeneDegrees]
    
    new_dict = [a for a, b in Stratification.items() for i2 in MDdegs if i2 in b]
    occurrences = Counter(new_dict).most_common()
    json.dump(occurrences, open('MDgeneDegDist.txt','w'))
#    writefile('MDgeneDegDist', occurrences)
    return occurrences

def RandomGeneList(Loops,loop2, p = None):
    import os
    import glob
    import random
    if p is None:
        RandomGL = []
    MDgenedist = json.load(open('MDgeneDegDist.txt'))    
    F3convert = json.load(open('AllGenesDict.txt'))
    for item46 in MDgenedist:
        for key,value in list(F3convert.items()):
            if key == item46[0]:
                print(key)
                randomgenes = random.sample(value, int(item46[1]))  
                RandomGL.append(randomgenes)
    RandGenes = list(sum(RandomGL[:], []))
    Fname = 'RandomGL_{0}_{1}'.format(Loops,loop2)
    writefilewithspace(Fname, RandGenes)
    return RandGenes
    #Obtain degree distribution
    #Obtain lengths of each 

