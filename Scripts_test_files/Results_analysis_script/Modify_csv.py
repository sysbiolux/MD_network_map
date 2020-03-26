# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 16:02:29 2019

This script takes as input the 2 .csv files, removes the single line spaces
between them and renames them as _mod.csv files. 

@author: apurva.badkas
"""
import os
import pandas as pd



filename = [file for file in os.listdir('.')]          #Collects list of filenames to concatenate data from
        
        
#F2convert = [item for item in filename if '3_RandListTable_' in item]   #Alternatively, 'BetweennessMeasures'
F3convert = [item00 for item00 in filename if '.csv' in item00]              #Only csv files
F4convert = sorted(F3convert)

for item01 in F4convert:
    DFrame = pd.read_csv(item01)
    DFrame = DFrame.drop(DFrame.columns[0], axis=1)
    Fname = item01[0:-4]
    ModFname = Fname + '_mod.csv'
    DFrame.to_csv(ModFname, index = False, header = False)
    
