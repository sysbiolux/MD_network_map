# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 10:33:54 2018

@author: apurva.badkas

This file takes BioGRID MV raw datafile and converts it into a two column format. It selects human interactions.
"""

#Function definitions
def readfile(filetoread):    
    listoflines = []
    openfile = open(filetoread)
    outputpar = openfile.readlines()
    for linein in outputpar:
        listoflines.append(linein)
    openfile.close()    
    listoflines[-1] = listoflines[-1]+ '\n'
    return listoflines


def writefilewospace(filenametowrite, datatowrite):
    datatowritewithspace = [item +'\n' for item in datatowrite]
    with open("%s.txt" %filenametowrite, "w+") as f: f.writelines("".join(datatowritewithspace))
#End of definitions


#biogridfile = readfile('BIOGRID-MV-Physical-3.5.166.tab2.txt')

biogridfile = readfile('Test.txt')
firstline = biogridfile[0].split('\t')
colind1 = [ind1 for ind1,val1 in enumerate(firstline) if 'Type' in val1]
PhysicalInts2 = [item2 for item2 in biogridfile if 'hysical' in item2.split('\t')[colind1[0]]]
colind = [ind for ind,val in enumerate(firstline) if 'Organism' in val]
PhysicalInts = [item for item in PhysicalInts2 if '9606' in item.split('\t')[colind[0]] and '9606' in item.split('\t')[colind[1]]]
PhysicalInts_extracted =[item4.split('\t')[7]+'\t'+item4.split('\t')[8] for item4 in PhysicalInts]
setphys = set(PhysicalInts_extracted)
Listphys = list(setphys)


#writefilewospace('Test_Human_fromOrgansim_PhysInts_Test',PhysicalInts)  #Intermediate file
#writefilewospace('Test_ExtractedBioGRID_PPI_8112019_Test',PhysicalInts_extracted) #Intermediate file
writefilewospace('BioGRID_HGNC_PPI_NR_8112019_Test',Listphys)
