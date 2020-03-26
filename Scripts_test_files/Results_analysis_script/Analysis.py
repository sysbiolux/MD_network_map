# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 23:34:46 2019

@author: Admin
"""

import os
import csv
import itertools
import sys
import math as mt
import operator

def writefilewospace(filenametowrite, datatowrite):
    with open("%s.txt" %filenametowrite, "w") as f: f.writelines("\n".join(datatowrite))

def readfile(filetoread, p=None):
    if p is None:
        listoflines = []
    openfile = open(filetoread)
    outputpar = openfile.readlines()
    listoflines = [linein for linein in outputpar if linein.split() != [] and linein.split() != '' ]
    openfile.close()
    listoflines[-1] = listoflines[-1] + '\n'
    return listoflines

def readcsvfile(filetoread):
    listoflines = []

    with open(filetoread) as csv_file:
        openfile = csv.reader(csv_file)
        for linein in openfile:
            if linein != [] :
                listoflines.append(linein)
    return listoflines


def writecentvalues(csvftowrite, datafromappendtable):
    with open('%s.csv'%csvftowrite,"w") as f: #Comb_File_Set4_12112018  "SelectedData_%s_%s" %(RType,No)
        wr = csv.writer(f)
        wr.writerows(datafromappendtable)
#End of definitions

SetNo = '6'
ONo = '03'
Betno = '03'
RunType = 'TEST'



def RandRunsDataSelection(Betno,SetNo,No0,RType, pqr = None):
    if pqr is None:
        indlist, indlistUncomon, LofGenes, LofGenes2, ListofGenes, ListofGenes2  = ([] for i in range(6))
#Read in the Masterlist of Nodes. Serves as reference to add selected values

    redlist = readfile('Allnodes.txt')  #Read Masterfile with all genes
    full = set(redlist)
    fulllist = list(full)
    SFulllist = sorted(fulllist) #Sorted Master genelist
    
    #Split and create a list as reference to add to
    ListofGenes = [item2.split() for item2 in SFulllist] # Remove newline character
    ListofGenes2 = [item3[0] for item3 in ListofGenes] # List of gene symbols
    
    #Data is appended to list of lists
    LofGenes = [[item] for item in ListofGenes2]  # List of lists for appending values
    LofGenes2 = [[item2] for item2 in ListofGenes2]
    
    #Read in input files: MD genelist and MD betweenness. Some MD genes have no betweenness, and 
    #the betweenness list has one step interactors, thus potential non-MD genes
    MDlist = readfile('MDgenes.txt')# Read in MD genelist
    MDBetfile = readfile('Betweenness.txt')
    MDBetsplt = [item2.split(':')[0] for item2 in MDBetfile]
    MDlistsplt = [item4.split()[0] for item4 in MDlist]
    MDUncommon = [item6 for item6 in MDBetsplt if item6 not in MDlistsplt]
    
    #Datafiles to select the data from. Randlists and Centrality values
    Randlist = readcsvfile('RandGenes.csv') # Read Randomlist table file
    Centralitytable = readcsvfile('BetCent.csv') # Read Centrality table file    
    SMDlist = sorted(MDlist)
    SMDUncommon = sorted(MDUncommon)
    indlist = [] # indlist is collecting all the indices of the MD genes in the ListofGenes
    indlistUncommon = []
    #Select for genes in the MD genelist. Collect indices of MD genes in Ref List
    if len(ListofGenes2) == len(LofGenes):
        for MDgene in SMDlist:
            MDgene2 = MDgene.split()
            indlist.append([l for l,val in enumerate(ListofGenes2) if val == MDgene2[0]])
    else: 
        print ('Lists mismatch. Unequal sizes')
        
    #Select for genes in the betweenness list, not in MD list. Collect indices of Non MD genes in bet file    
    if len(ListofGenes2) == len(LofGenes):
        for MDUncommongene in SMDUncommon:
            MDUncommongene2 = MDUncommongene.split()
            indlistUncommon.append([i2 for i2,val2 in enumerate(ListofGenes2) if val2 == MDUncommongene2[0]])
    else: 
        print ('Lists mismatch. Unequal sizes')
    
    #Remove empty lists in the list of indices
    Indlist = indlist
    Nonemptylist = [item5 for item5 in Indlist if item5 != []]
    #Merge all elements into a single list
    c = list(sum(Nonemptylist[:], []))
    d = list(sum(indlistUncommon[:], [])) #Merge items into a single list to loop through
    
    
    #Select data based on the indices
    for ij in c:
        if Randlist[ij][0].split() != '' and Centralitytable[ij][0].split() != '':
          del Randlist[ij][0]  # del gene name
          del Centralitytable[ij][0]   #del gene name
          for ind, item1 in enumerate(Randlist[ij]):
              if item1 != 'Ab' :
                  LofGenes[ij].append(Centralitytable[ij][ind])
              else:
                  LofGenes2[ij].append(Centralitytable[ij][ind])
              
    for ik in d:
        if Randlist[ik][0].split() != '' and Centralitytable[ik][0].split() != '':
          del Randlist[ik][0]  # del gene name
          del Centralitytable[ik][0]   #del gene name
          for ind3, item3 in enumerate(Randlist[ik]):
              if item3 == 'Ab' :
                  LofGenes[ik].append(Centralitytable[ik][ind3])
              else:
                  LofGenes2[ik].append(Centralitytable[ik][ind3])
        
    
    writecentvalues("SelectedData_%s_%s" %(RType,ONo),LofGenes)


def InputProcessing(RType, Setno, ONo, pt = None): #'NonEmptylist_%s_%s' %(RType,No)
    if pt is None:
        item2, Actlist, lstimproved, Nonemptylist = ([] for i in range(4))
    
    Readfile = readcsvfile('SelectedData_%s_%s.csv'%(RType,ONo)) # Read input file. 9545 entries 'SelectedData_%s_%s.csv'%(RType,No)

#First stage filtering: Collect entries with data. Returns 1896 entries. The MD genes
    item2 = [item1 for item1 in Readfile if len(item1)>1]

# Filter out genes which have data, discard lines with only genenames

#Second stage filtering: Remove all entries of 'Absent'.
    Actlist = []
    for lst in item2:
        lstimproved = [item4 for item4 in lst if item4 != 'Ab'] #Append data that isnt 'Absent'
        Actlist.append(lstimproved)    #List with all genes with some data

#Third stage filtering: Remove items with no values. ## No values bec earlier the only entries were Absent which were removed    
    Nonemptylist = [item6 for item6 in Actlist if len(item6)>1]

    writecentvalues('Actlist', Actlist)
    writecentvalues('NonEmptylist_%s_%s' %(RType,ONo), Nonemptylist)


def PValues(RType,SNo,No,betno, lt = None):
    if lt is None:
        MDgenelist, MDSpecGenes, MDgenename2, Selgenes, Selcentvals, l1, l2,l3, l4, PTable, = ([] for k in range(10))
    CentvaluesTable = readcsvfile('NonEmptylist_%s_%s.csv'%(RType,No))
    MDBetweenness = readfile('Betweenness.txt')
    MDOriglist = readfile('Allnodes.txt')
    MDspltList = [item.split()[0] for item in MDOriglist]

    MDgenelist = []
    TempVals = []
    MDSpecGenes = []

    MDgenename2 = [lineinMD3.split(':')[0]+','+ lineinMD3.split(':')[1] for lineinMD3 in MDBetweenness]

    Selgenes = [item for item in CentvaluesTable for item2 in MDgenename2 if item[0]== item2.split(',')[0]]
    Selcentvals = [item3 for item3 in MDgenename2 for item4 in Selgenes if item3.split(',')[0] == item4[0]]
    l1 = sorted(Selgenes)
    l2 = sorted(Selcentvals)

    l4 = []
    for x,y in zip(l2,l1):
        l3 = [x.split(',')[0],x.split(',')[1][:-1],y[1:]]
        l4.append(l3)


    PTable = []
    for value in l4:
        TempVals = []
        L= []
        Av = []
        Relbet = []
        TempVals = [value1 for value1 in value[2] if float(value1) >= float(value[1])]
        Numofhits = len(TempVals)
        TotalNo = len(value[2])
        pvalue = (float(Numofhits) + 1)/(float(TotalNo) + 1)
        L = [float(n) for n in value[2]]
        Av = sum(L)/float(TotalNo)
        if Av != 0:
            Relbet = float(value[1])/float(Av)
    
        PTable.append([value[0],str(Numofhits),str(TotalNo),str(pvalue),Av,value[1],Relbet,])    
    writecentvalues('Pvals_%s_1_%s'%(RType,ONo),PTable)
    

def main_runs():
    RandRunsDataSelection(Betno, SetNo, ONo, RunType)
    InputProcessing(RunType, SetNo, ONo) #RType, Setno, ONo
    PValues(RunType, SetNo,ONo, Betno)

main_runs()