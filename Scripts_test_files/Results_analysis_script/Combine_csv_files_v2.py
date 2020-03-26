# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 15:53:48 2019
File used to combine Randlist output tables. Tested version
@author: apurva.badkas
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
import csv



def readcsvfile(filetoread):    
    listoflines = []
    with open(filetoread) as csv_file:
        openfile = csv.reader(csv_file)
        listoflines = [linein for linein in openfile if linein != []]
    return listoflines


def writecentvalues(datafromappendtable):    
    with open("Comb_RandGeneTable_S1_20.csv","w") as f: #Comb_File_Set4_12112018 Comb_RandGeneTable_S1_12
        wr = csv.writer(f)
        wr.writerows(datafromappendtable)
        
filename = [file for file in os.listdir('.')]          #Collects list of filenames to concatenate data from
        
        
F2convert = [item for item in filename if '3_RandListTable_' in item]   #Alternatively, 'BetweennessMeasures' 3_BetweenessMeasures_01_mod
F3convert = [item for item in F2convert if '.csv' in item]              #Only csv files
F4convert = sorted(F3convert)


firstfile = readcsvfile(F4convert[0]) #First file in list. Rest of the files are appended to this

for item in F4convert[1:]:
    f1 = readcsvfile(item)
    if len(firstfile) == len(f1):  
        lst = []
        d = []
        for x1,x2 in zip(firstfile,f1):
          temp = [x1, x2[1:]]   #Skip the gene name in the first columns, only append values
          lst.append(temp)   
    d = [list(sum(item1[:], [])) for item1 in lst]     
    firstfile = d                  #Keep appending to previous updated file

writecentvalues(firstfile)    