# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 14:19:26 2019

@author: Admin
"""
import time
import csv
import pickle
import json
import networkx as nx

def readfile(filetoread, p=None):
    if p is None:
        listoflines = []
    openfile = open(filetoread)
    outputpar = openfile.readlines()
    listoflines = [linein for linein in outputpar if linein.split() != [] and linein.split() != '']
    openfile.close()
    listoflines[-1] = listoflines[-1] + '\n'
    return listoflines

#Writing files: without and with space
def writefilewospace(filenametowrite, datatowrite):
    with open("%s.txt" %filenametowrite, "w") as f: f.writelines("".join(datatowrite))

def writefile(filenametowrite, datatowrite):
    with open("%s.txt" %filenametowrite, "w") as f: f.writelines("\n".join(datatowrite))

def writefilewithspace(filenametowrite, datatowrite):
    datatowritewithspace = [item +'\n' for item in datatowrite]
    with open("%s.txt" %filenametowrite, "wb+") as f: f.writelines("".join(datatowritewithspace))

#PPI to Dictionary conversion


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
    G = nx.read_edgelist("%s.txt" %filenametoread)   # Using the complete PPI from above 
    Gcc=sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)
    G0=Gcc[0]
    nx.write_edgelist(G0, "%s_GiantComponent.txt" %filenametowrite)

#Extract all the nodes in the PPI
def NodeExtraction(PPIfname, Allnodesfname, ijk = None):
    if ijk is None:
        GL,G1,G2,G3,G4,G5 = ([] for i in range(6))
    Pfile = readfile('%s.txt' %PPIfname)
    GL = [item.split() for item in Pfile]
    G1 = [item2[0] for item2 in GL]
    G2 = [item3[1] for item3 in GL]
    G3 = G1 + G2

    G4 = set(G3)
    G5 = list(G4)

    writefile(Allnodesfname,G5)

#PPI to dict conversion 
def MDgenesMappingtoPPI(readMDgenesname,PPInodesfile,fnametowrite):
    MDgeneList = readfile(readMDgenesname)
    PPInodes = readfile(PPInodesfile)
    MD = [item23.split()[0] for item23 in MDgeneList]
    PPINode = [item24.split()[0] for item24 in PPInodes]
    MDnodesmapping = [item25 for item25 in MD if item25 in PPINode]
    Missing = [item26 for item26 in MD if item26 not in PPINode]
    Missingf = 'Missing_HPRD_nonmapping'
    writefile(fnametowrite,MDnodesmapping)
    writefile(Missingf, Missing)
    return MDnodesmapping
    
    
def save_dict_to_file(dic):
    f = open('dict.txt','w')
    f.write(str(dic))
    f.close()

def load_dict_from_file():
    f = open('dict.txt','r')
    data=f.read()
    f.close()
    return eval(data)
def writecentvalues(datafromappendtable):
    print('\tWriteCentVals')
    start = time.time()
    with open("Sorted_PPI.csv","w") as f:
        wr = csv.writer(f)
        wr.writerows(datafromappendtable)
    print("\t\tTime: %.4F" % (time.time() - start))

def PPItodict(Allnodesfname, PPIfname):
    res = dict()
    lst = readfile('%s.txt'%Allnodesfname)
    ppifile = readfile('%s.txt' %PPIfname)
    NRset = set(ppifile)
    NRPPI = list(NRset)
    #dictOfWords = dict.fromkeys(lst)
    Lst1 = [item.split() for item in lst]
    Lstgenes = [item.split() for item in lst]
    print('\tBuilding_set')
    start = time.time()
    createppilst =[item2.append(linein) for item2 in Lstgenes for linein in NRPPI if item2[0]==linein.split()[0] or item2[0]==linein.split()[1]]
    print("\t\tTime: %.4F" % (time.time() - start))
    
    writecentvalues(Lstgenes)
    with open('filepickle_NR_HPRD_GC_PPI.txt', 'wb') as handle:
        pickle.dump(Lstgenes, handle)
    with open('filepickle_NR_HPRD_GC_PPI.txt', 'rb') as handle:
        b = pickle.loads(handle.read())
    for item8 in b:
        res[item8[0]] = item8[1:]
    with open('PPI_dict_HPRD_GC_NR.txt', 'wb') as handle:
        pickle.dump(res, handle)

def PPIfileformat(fnametoread, fnametowrite):
    Pfiletoconvert = readfile(fnametoread)
    Datatowrite = [item27.split('{}')[0] for item27 in Pfiletoconvert]
    writefile(fnametowrite, Datatowrite)
    

PPIfile = 'PPItodict_Test'
createGiantComponent(PPIfile,'HPRD_')
NodeExtraction('HPRD__GiantComponent','All_Nodes_HPRD')
PPItodict('All_Nodes_HPRD','HPRD__GiantComponent' )
MDgenesMappingtoPPI('MD_HGNC_NR_Test.txt','All_Nodes_HPRD.txt','MDNodesMappingtoHPRD')   
PPIfileformat('HPRD__GiantComponent.txt','HPRD_GC_PPI')
