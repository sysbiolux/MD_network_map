# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 18:34:21 2019

@author: Admin
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor
File name: Masterscript_v6
Replaces all previous versions of Masterscripts

Extraction of disease-specific network based on a gene-list.
Changes:
    1. Many of the for loops replaced by list comprehensions
    2. Parallelized version of betweenness centrality algorithm

Inputs: 1. Disease-gene list
        2. PPI file
        3. Bet_cent_parallel to run the betweenness centrality algorithm
WARNING: This script needs to be run in an external window. Spyder has issues with multiprocessing module    

Notes: For AppendRandlisttable, main1() returns the randlist, which is collected by the variable f1
and f1 is then supplied as an argument to main3(f1) for the table to get appendend with the current
genelist

Options for runs:
    1. No node removal
    2. Node removal with no filtering
    3. Node removal with Rand filtering
    4. Node removal with MD filtering

"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 20:19:28 2019
@author: Admin
"""


from collections import OrderedDict
import networkx as nx
import time
import pickle

import matplotlib.pyplot as plt
import Bet_cent_parallel_2


ts = time.gmtime()
print(time.strftime("%Y-%m-%d %H:%M:%S", ts))




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

def writefile(filenametowrite, datatowrite):
    with open("%s.txt" %filenametowrite, "w") as f: f.writelines("\n".join(datatowrite))

def writefilewospace(filenametowrite, datatowrite):
    with open("%s.txt" %filenametowrite, "w") as f: f.writelines("".join(datatowrite)) 
    

"""
#PPIextraction: Extracts disease specfic network based on the input diseasegene list from the PPI file
Inputs: diseasegenelist - Initial seed list. Input to the function as list after reading the file
        ppifile - PPI file read
        filetosave - basefilename
"""
def PPIextraction(diseasegenelist,filetosave,l=None):
    if l is None:
        genelist2 = []
        genelist = [] 
        ppilist = [] 
        result = []
        result3 = []
        result2 = []
    genelist2 = diseasegenelist 
    genelist = [item15 for item15 in genelist2 if item15 != '\n']
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
    result2 = set(result)           
    result3 = list(result2)                 
    print("\t\tTime: %.4F" % (time.time() - start))                
    filenamebydb = "%s_PPI" %filetosave #"17082018_%s" %dbname
    writefile(filenamebydb,result2)
    return result3
    
def AddConnectedComponents(res3,diseasegenelist,ppifile,filetosave, lr = None):
    if lr is None:
        MDgenelist = []
        PPIfullfile = []
        MDlist = []
        Pfile = []
        GL,G1,G2,G3,G4,G5,FNeighbors,PPIfilt,PPIsel,TotPPI_red, TotPPIset = ([] for i in range(11))
       
    MDgenelist = diseasegenelist  
    PPIfullfile = ppifile
    MDlist = [item01.split()[0] for item01 in MDgenelist]
    Pfile = res3

    GL = [item.split() for item in Pfile]
    G1 = [item2[0] for item2 in GL]
    G2 = [item3[1] for item3 in GL]
    G3 = G1 + G2
    G4 = set(G3)
    G5 = list(G4)

    FNeighbors = [item for item in G5 if item not in MDlist]

    PPIfilt = [item03 for item03 in PPIfullfile for item04 in FNeighbors if item04 in item03]
    PPIsel = [item02 for item02 in PPIfullfile if item02.split()[0] in FNeighbors and item02.split()[1] in FNeighbors]    
    TotPPI_red = res3 + PPIsel
    TotPPIset = set(TotPPI_red)
    TotPPI = list(TotPPIset)    
    filenamebydb = "%s_TotPPI" %filetosave
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
    G = nx.read_edgelist("%s_TotPPI.txt" %filenametoread)   
    Gcc=sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)
    G0=Gcc[0]
    nx.write_edgelist(G0, "%s_GiantComponent.txt" %filenametowrite)

def RemoveSingleInteractions(filenametoread,l=None):
    if l is None:
        removedInteraction = []
        intlist = []
    GCfileopen = open('%s_GiantComponent.txt'%filenametoread)
    GCfileread = GCfileopen.readlines()

    for line11 in GCfileread:
        inp2 = line11.split()
        if '-' in inp2[0] and len(inp2[0])==1:
            removedInteraction.append(line11)
        elif '-' in inp2[1] and len(inp2[1])==1:
            removedInteraction.append(line11)
        else:
            intlist.append(line11)
    with open("%s_FilteredGiantComponent.txt" %filenametoread, "w") as f2: f2.write(''.join(intlist))
    GCfileopen.close()  
 
               
def calculateDegree(filenametoread,filenametowrite):    
    G = nx.read_edgelist("%s_FilteredGiantComponent.txt" %filenametoread)
    G.remove_edges_from(G.selfloop_edges())
    nx.write_edgelist(G, "%s_processed_GC.txt"%filenametowrite)
    d = dict(G.degree)
    n = OrderedDict(sorted(d.items(), key=lambda t: t[1], reverse=True))
    with open ('%s_degree.txt'%filenametowrite, 'w') as fp:
        for p in n.items():
            fp.write("%s:%s\n" % p)
            
"""
processDegreefile_b and NewGiantComponent_b applicable for case of node removal without filtering
"""            

def processDegreefile_b(diseasegenelist, filenametoread,filenametowrite,l=None):
    if l is None:
        linestoread = []
        rejectedlist = []
        Rejlist = []
        Rejgenelist =[]
        genelistwSpace = []
        genelist = []
        rejectedlist_1 = []
    genelistwSpace = diseasegenelist 
    genelist = [item10 for item10 in genelistwSpace if item10 != '\n']
    print (genelist[0:5])
    filename = '%s_degree.txt'%filenametoread
    linestoread = readfile(filename)
    rejectedlist_1 = [item3.split(':')[0] for item3 in linestoread if item3.split(':')[1].split()[0]=='1' or item3.split(':')[1].split()[0]=='0' ]
    print(rejectedlist_1[0:5])
    Rejlist = set(rejectedlist_1)
    Rejgenelist = list(Rejlist)    
    Removegenes = "GeneswDeglessthan1_%s" %filenametowrite
    writefile(Removegenes,Rejgenelist)



def NewGiantComponent_b(filenametoread,filenametowrite, p2=None):
    if p2 is None:
        removedGC1 = []
        removedGC2 = []
        readGCfile = []
        dxy = []
        removedGC = []
        filtGCfile = []
        readRejgenefile = []
        remGCset = []
    print('\tNew GC')
    start = time.time()  
    readGCfile = readfile("%s_processed_GC.txt" %filenametoread)
    readRejgenefile = readfile("GeneswDeglessthan1_%s.txt" %filenametoread) 
    filtGCfile = [iteminGC for iteminGC in readGCfile for remgene in readRejgenefile if remgene.split()[0] in iteminGC ]
    removedGC1 = [item4 for item4 in filtGCfile for item in readRejgenefile if item.split()[0] == item4.split()[0]]
    removedGC2 = [item6 for item6 in filtGCfile for item7 in readRejgenefile if item7.split()[0] == item6.split()[1]]
    removedGC = removedGC1 + removedGC2
    remGCset = set(removedGC)
    remGClst = list(remGCset)
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

def processDegreefile_c(diseasegenelist, filenametoread,filenametowrite,l=None):
    if l is None:
        linestoread = []
        rejectedlist = []
        Rejlist = []
        Rejgenelist =[]
        genelistwSpace = []
        genelist = []
        rejectedlist_1 = []
    genelistwSpace = set(diseasegenelist)
    genelst = list(genelistwSpace)
    genelist = [item10 for item10 in genelst if item10 != '\n']
    print (genelist[0:5])
    filename = '%s_degree.txt'%filenametoread
    linestoread = readfile(filename)
    rejectedlist_1 = [item3.split(':')[0] for item3 in linestoread if item3.split(':')[1].split()[0]=='1' or item3.split(':')[1].split()[0]=='0' ]
    print(rejectedlist_1[0:5])
    rejectedlist = [item9 for item9 in rejectedlist_1 for item10 in genelist if item10.split() != [] if item9 == item10.split()[0] ]
    print("Secondlist",rejectedlist[0:5])
    trim = [rejectedlist_1.remove(item11) for item11 in rejectedlist]
    Rejlist = set(rejectedlist_1)
    Rejgenelist = list(Rejlist)    
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
    filtGCfile = [iteminGC for iteminGC in readGCfile for remgene in readRejgenefile if remgene.split()[0] in iteminGC ]
    removedGC1 = [item4 for item4 in filtGCfile for item in readRejgenefile if item.split()[0] == item4.split()[0]]
    removedGC2 = [item6 for item6 in filtGCfile for item7 in readRejgenefile if item7.split()[0] == item6.split()[1]]
    removedGC = removedGC1 + removedGC2
    remGCset = set(removedGC)
    remGClst = list(remGCset)
    dxy = [readGCfile.remove(item) for item in remGClst]
    print("\t\tTime: %.4F" % (time.time() - start))
    print('\tNew GC-write')
    start = time.time()
    writefilewospace("RevisedGC2_%s_GiantComponent"%filenametowrite,readGCfile)
    print("\t\tTime: %.4F" % (time.time() - start))

"""
processDegreefile_d and NewGiantComponent_d applicable for the case of single node removal except
MD genes

"""

def processDegreefile_d(filenametoread,filenametowrite,l=None):
    if l is None:
        linestoread = []
        rejectedlist = []
        Rejlist = []
        Rejgenelist =[]
        trim = []
        genelist = []
        rejectedlist_1 = []
    genelistMD = readfile('90pc_S1.txt')    
    genelistwSpace = set(genelistMD)
    genelst = list(genelistwSpace)
    genelist = [item10 for item10 in genelst if item10 != '\n']
    print (genelist[0:5])
    filename = '%s_degree.txt'%filenametoread
    linestoread = readfile(filename)
    rejectedlist_1 = [item3.split(':')[0] for item3 in linestoread if item3.split(':')[1].split()[0]=='1' or item3.split(':')[1].split()[0]=='0' ]
    print(rejectedlist_1[0:5])
#    rejlist = [item12 for item12 in genelist if item12 =! '']
    rejectedlist = [item9 for item9 in rejectedlist_1 for item10 in genelist if item10.split() != [] if item9 == item10.split()[0] ]
    print("Secondlist",rejectedlist[0:5])
    trim = [rejectedlist_1.remove(item11) for item11 in rejectedlist]
    Rejlist = set(rejectedlist_1)
    Rejgenelist = list(Rejlist)    
    Removegenes = "GeneswDeglessthan1_%s" %filenametowrite
    writefile(Removegenes,Rejgenelist)



def NewGiantComponent_d(filenametoread,filenametowrite, p2=None):
    if p2 is None:
        removedGC1 = []
        removedGC2 = []
        readGCfile = []
        dxy = []
        removedGC = []
        filtGCfile = []
        remGCset = []
    print('\tNew GC')
    start = time.time()  
    readGCfile = readfile("%s_processed_GC.txt" %filenametoread)
    readRejgenefile = readfile("GeneswDeglessthan1_%s.txt" %filenametoread) 
    filtGCfile = [iteminGC for iteminGC in readGCfile for remgene in readRejgenefile if remgene.split()[0] in iteminGC ]
    removedGC1 = [item4 for item4 in filtGCfile for item in readRejgenefile if item.split()[0] == item4.split()[0]]
    removedGC2 = [item6 for item6 in filtGCfile for item7 in readRejgenefile if item7.split()[0] == item6.split()[1]]
    removedGC = removedGC1 + removedGC2
    remGCset = set(removedGC)
    remGClst = list(remGCset)
    dxy = [readGCfile.remove(item) for item in remGClst]
    print("\t\tTime: %.4F" % (time.time() - start))
    print('\tNew GC-write')
    start = time.time()
    writefilewospace("RevisedGC2_%s_GiantComponent"%filenametowrite,readGCfile)
    print("\t\tTime: %.4F" % (time.time() - start))


    
def main1c():
    Diseasegenelist = readfile(RandomGL)
    PPIext = PPIextraction(Diseasegenelist,basefilename)  #Geneslist_NR2,Geneslist_NR2_16082018,1920_unique_genes 
    CC = AddConnectedComponents(PPIext,Diseasegenelist,PPIfile,basefilename)
    FirstGC = createGiantComponent(basefilename,basefilename)
    RemSingInt = RemoveSingleInteractions(basefilename)
    Firstdeg = calculateDegree(basefilename,basefilename)
    FilteredDeg = processDegreefile_c(Diseasegenelist,basefilename,basefilename)      
    SecondGC = NewGiantComponent_c(basefilename,basefilename)
    return Diseasegenelist

def main2b():    
    reffilename = "RevisedGC2_%s" %basefilename
    Bet_cent_parallel_2.runcent(reffilename,jobID)
    
#def main3b(ijk):    
#    tableinput = AppendTable(betweenName)
#    randlistinput= AppendRandListTable(ijk)
#    time.sleep(2)        
 
"""
processDegreefile_c and NewGiantComponent_c applicable for case of removing single degree
nodes except the ones in rand gene lists
"""    
    
dbname = "HPRD" #Option: HPRD, BioGRID,MINT    
Op = [3]
Setno ='Full'
#SetNo = '1'
#ONo = '01'
#betno = '01'
RunType = 'S%s_wCC'%Setno

for k in range(1):
    Option = Op[k]
    jobID = '0%s'%str(Option)
#jobID =sys.argv[1]
    print ("Saving to file " + jobID)
    basefilename = jobID + "CV_wCC_%s" %dbname      #Used to name all the files generated
    betweenName = jobID + "betweeness"      #Naming convention, required while running multiple jobs               
    degreeName = jobID + "degree"
    
    for c in range(1):
        timestamp = time.time()
        PPIfile = readfile('HPRD_GC_PPI.txt')
 
                    
        if Option == 3:
            for i in range(1):
                RandomGL = 'MDNodesMappingtoHPRD.txt'              
                if __name__ == "__main__":
                    f1= main1c()
                    main2b()
